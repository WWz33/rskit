from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from rskit.config import DESeq2Config
from rskit.core.salmon import SalmonExpressionExporter, merge_salmon_quant_tables
from rskit.input_contracts import (
    design_columns,
    ensure_genes_by_samples,
    load_coldata,
    read_table,
    require_no_missing_values,
    validate_sample_alignment,
)
from rskit.utils.logger import get_logger
from rskit.utils.manifest import write_manifest

logger = get_logger(__name__)


def parse_contrast(contrast_value: Optional[str], metadata_df: pd.DataFrame) -> Optional[List[str]]:
    """Parse and validate a DESeq2 contrast against sample metadata."""
    if not contrast_value:
        return None

    contrast = [part.strip() for part in contrast_value.split(",")]
    if len(contrast) != 3 or any(not part for part in contrast):
        raise ValueError("Contrast must be in format: factor,level1,level2 (e.g., condition,B,A)")

    factor, level1, level2 = contrast
    if factor not in metadata_df.columns:
        raise ValueError(
            f"Contrast factor '{factor}' is not present in coldata columns: "
            + ", ".join(metadata_df.columns)
        )

    levels = sorted(str(level) for level in metadata_df[factor].dropna().unique())
    missing_levels = [level for level in (level1, level2) if level not in levels]
    if missing_levels:
        raise ValueError(
            f"Contrast levels not found in coldata column '{factor}': "
            + ", ".join(missing_levels)
            + ". Available levels: "
            + ", ".join(levels)
        )

    return contrast


def parse_contrasts(contrast_values, metadata_df: pd.DataFrame) -> Optional[List[List[str]]]:
    """Parse and validate one or more DESeq2 contrasts against sample metadata.

    Accepts a single string or a sequence of strings (the CLI passes one entry
    per -c). Duplicates are rejected: they would produce colliding outputs for
    no benefit.
    """
    if not contrast_values:
        return None
    if isinstance(contrast_values, str):
        contrast_values = [contrast_values]

    contrasts: List[List[str]] = []
    for value in contrast_values:
        contrast = parse_contrast(value, metadata_df)
        if contrast is None:
            raise ValueError("Contrast must be in format: factor,level1,level2 (e.g., condition,B,A)")
        if contrast in contrasts:
            raise ValueError("Duplicate contrast: " + ",".join(contrast))
        contrasts.append(contrast)
    return contrasts


def _lfc_shrink_coefficient(contrast: List[str], lfc_columns) -> Optional[str]:
    """Resolve the pydeseq2 LFC coefficient column for a (factor, tested, reference) contrast.

    pydeseq2 names LFC columns after the design matrix (e.g. 'condition[T.treat]' via
    formulaic), so match the contrast against known naming schemes instead of guessing.
    """
    factor, tested, reference = contrast[0], contrast[1], contrast[2]
    candidates = (
        f"{factor}[T.{tested}]",
        f"{factor}_{tested}_vs_{reference}",
        f"{factor}[{tested}]",
        factor,
    )
    return next((name for name in candidates if name in lfc_columns), None)


class Deseq2Analyzer:
    def __init__(self, config: DESeq2Config):
        self.config = config
        self.logger = logger
        self.expression_exporter = SalmonExpressionExporter()
        self.dds = None
        self.inference = None
        self.stats_results = None
        self.stat_res = None
        self.contrast = None
        self.counts_df = None
        self.metadata_df = None
        
    def _create_tx2gene_from_gtf(self, gtf_file: str, output_dir: Optional[str] = None) -> pd.DataFrame:
        """Create transcript-to-gene mapping from GTF/GFF3 file."""
        tx2gene_df = self.expression_exporter._create_tx2gene_from_gtf(gtf_file, output_dir)
        self.logger.info(f"tx2gene preview:\n{tx2gene_df.head()}")
        return tx2gene_df
    
    def load_counts_from_file(self, counts_file: str, metadata_df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """Load count data from file.
        
        Args:
            counts_file: Path to counts matrix file (genes x samples)
            
        Returns:
            DataFrame with counts (samples x genes)
        """
        counts_df = read_table(counts_file, index_col=0)

        if metadata_df is not None:
            counts_df = ensure_genes_by_samples(counts_df, metadata_df, table_name="counts matrix")
        else:
            counts_df = counts_df.T

        if counts_df.isna().any().any():
            raise ValueError(
                f"Counts matrix {counts_file} contains NaN values; "
                "check for missing entries before running DESeq2"
            )

        # Ensure integer counts
        counts_df = counts_df.round().astype(int)
        
        self.counts_df = counts_df
        self.logger.info(f"Loaded counts for {counts_df.shape[0]} samples and {counts_df.shape[1]} genes")
        
        return counts_df

    def prefilter_counts(self, counts_df: pd.DataFrame) -> pd.DataFrame:
        """Keep genes whose total counts across samples meet the configured threshold."""
        min_count = self.config.prefilter_min_count
        if min_count <= 0:
            return counts_df

        genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= min_count]
        filtered_counts = counts_df.loc[:, genes_to_keep]
        removed_genes = counts_df.shape[1] - filtered_counts.shape[1]

        if filtered_counts.shape[1] == 0:
            raise ValueError(
                f"No genes remain after prefiltering with total count >= {min_count}"
            )

        if removed_genes:
            self.logger.info(
                f"Prefiltered {removed_genes} genes with total counts < {min_count}; "
                f"{filtered_counts.shape[1]} genes remain"
            )

        return filtered_counts
    
    def load_metadata(self, metadata_file: str, required_columns: Optional[List[str]] = None) -> pd.DataFrame:
        """Load metadata from CSV or TSV file.
        
        Expected format:
            sample,id,condition
            lhy-D-rep1,lhy_D,lhy-D
            lhy-D-rep2,lhy_D,lhy-D
            ...
        
        The 'sample' column is used as index to match with Salmon output directories.
        The 'id' column can be used for batch/group effects.
        The 'condition' column is used for differential expression analysis.
        
        Args:
            metadata_file: Path to metadata file (coldata)
            
        Returns:
            DataFrame with sample metadata (index=sample names)
        """
        metadata_df = load_coldata(metadata_file, required_columns=required_columns or [])
        require_no_missing_values(metadata_df, required_columns or [])
        
        self.metadata_df = metadata_df
        self.logger.info(f"Loaded metadata for {metadata_df.shape[0]} samples")
        self.logger.info(f"Columns: {list(metadata_df.columns)}")
        self.logger.info(f"Conditions: {dict(metadata_df['condition'].value_counts()) if 'condition' in metadata_df.columns else 'N/A'}")
        
        return metadata_df
    
    def fit(self, counts_df: Optional[pd.DataFrame] = None,
            metadata_df: Optional[pd.DataFrame] = None,
            contrast_factors: Optional[List[str]] = None) -> None:
        """Build and fit the DESeq2 model once.

        ``contrast_results()`` can then run any number of contrasts against
        this single fit, which is the expensive part of the analysis.

        Args:
            counts_df: Count matrix (samples x genes), or use self.counts_df
            metadata_df: Sample metadata, or use self.metadata_df
            contrast_factors: metadata columns that will be contrasted later;
                coerced to strings so formulaic treats them as categorical
        """
        try:
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.default_inference import DefaultInference
        except ImportError:
            raise ImportError("PyDESeq2 is not installed. Please install it with: pip install pydeseq2")

        # Use provided data or stored data
        if counts_df is not None:
            self.counts_df = counts_df
        if metadata_df is not None:
            self.metadata_df = metadata_df

        if self.counts_df is None:
            raise ValueError("No counts data provided. Load counts first.")
        if self.metadata_df is None:
            raise ValueError("No metadata provided. Load metadata first.")

        # Contrasts compare factor levels, so every contrasted factor must be
        # categorical; numeric coldata columns (condition = 0/1) fail inside
        # formulaic otherwise. Coerce before the design matrix is built.
        for factor in contrast_factors or ["condition"]:
            if factor in self.metadata_df.columns and self.metadata_df[factor].map(
                lambda value: not isinstance(value, str)
            ).any():
                self.logger.info(f"Coercing coldata column '{factor}' to string for contrast levels")
                self.metadata_df[factor] = self.metadata_df[factor].astype(str)

        # Ensure sample names match
        validate_sample_alignment(self.counts_df, self.metadata_df, table_name="counts matrix")
        self.counts_df = self.counts_df.loc[self.metadata_df.index]
        self.counts_df = self.prefilter_counts(self.counts_df)

        # Initialize inference
        self.inference = DefaultInference(n_cpus=self.config.n_cpus)

        # Create DeseqDataSet
        self.logger.info("Creating DeseqDataSet...")
        self.dds = DeseqDataSet(
            counts=self.counts_df,
            metadata=self.metadata_df,
            design=self.config.design,
            fit_type=self.config.fit_type,
            size_factors_fit_type=self.config.size_factors_fit_type,
            refit_cooks=self.config.refit_cooks,
            min_replicates=self.config.min_replicates,
            inference=self.inference,
            quiet=self.config.quiet
        )

        # Run DESeq2 pipeline
        self.logger.info("Running DESeq2 pipeline...")
        self.dds.deseq2()

    def contrast_results(self, contrast: Optional[List[str]] = None) -> pd.DataFrame:
        """Run one contrast (Wald test + LFC shrinkage) on the fitted model.

        Args:
            contrast: Contrast specification ['condition', 'B', 'A']; inferred
                from the design when omitted

        Returns:
            DataFrame with differential expression results for this contrast
        """
        try:
            from pydeseq2.ds import DeseqStats
        except ImportError:
            raise ImportError("PyDESeq2 is not installed. Please install it with: pip install pydeseq2")

        if self.dds is None:
            raise ValueError("No fitted model. Call fit() first.")

        # Set default contrast if not provided
        if contrast is None:
            contrast = self._infer_contrast()

        # Create stats object
        self.logger.info(f"Running statistical analysis with contrast: {contrast}")
        stat_res = DeseqStats(
            dds=self.dds,
            contrast=contrast,
            alpha=self.config.alpha,
            cooks_filter=self.config.cooks_filter,
            independent_filter=self.config.independent_filter,
            lfc_null=self.config.lfc_null,
            alt_hypothesis=self.config.alt_hypothesis,
            inference=self.inference,
            quiet=self.config.quiet
        )

        # Run Wald test
        stat_res.summary()

        # Store results
        self.stat_res = stat_res
        self.contrast = contrast
        self.stats_results = stat_res.results_df

        # Apply LFC shrinkage
        try:
            coeff = _lfc_shrink_coefficient(contrast, stat_res.LFC.columns)
            if coeff is None:
                self.logger.warning(
                    f"Could not apply LFC shrinkage: no coefficient for contrast {contrast} "
                    f"in LFC columns {list(stat_res.LFC.columns)}"
                )
            else:
                stat_res.lfc_shrink(coeff=coeff)
                self.logger.info(f"LFC shrinkage applied for coefficient: {coeff}")
        except Exception as e:
            self.logger.warning(f"Could not apply LFC shrinkage: {e}")

        return self.stats_results

    def analyze(self, counts_df: Optional[pd.DataFrame] = None,
                metadata_df: Optional[pd.DataFrame] = None,
                contrast: Optional[List[str]] = None) -> pd.DataFrame:
        """Fit the model and run a single contrast.

        Convenience wrapper around fit()/contrast_results() kept for API
        callers; multiple contrasts should call those two directly so the fit
        is shared.

        Args:
            counts_df: Count matrix (samples x genes), or use self.counts_df
            metadata_df: Sample metadata, or use self.metadata_df
            contrast: Contrast specification ['condition', 'B', 'A']

        Returns:
            DataFrame with differential expression results
        """
        self.fit(
            counts_df,
            metadata_df,
            contrast_factors=[contrast[0]] if contrast else ["condition"],
        )
        return self.contrast_results(contrast)
    
    def _infer_contrast(self) -> List[str]:
        """Infer contrast from design matrix and metadata.

        Prefers a 'condition' column; otherwise uses the first non-intercept
        design column that names a categorical metadata column.
        """
        # First, try to use 'condition' column if it exists
        if 'condition' in self.metadata_df.columns:
            unique_vals = self.metadata_df['condition'].unique()
            if len(unique_vals) >= 2:
                # Sort to ensure consistent ordering
                sorted_vals = sorted([str(v) for v in unique_vals])
                self.logger.info(f"Using 'condition' column for contrast: {sorted_vals[-1]} vs {sorted_vals[0]}")
                return ['condition', sorted_vals[-1], sorted_vals[0]]

        # Otherwise look for a usable design-matrix column
        design_cols = self.dds.obsm["design_matrix"].columns
        if len(design_cols) <= 1:
            raise ValueError("Design matrix has only intercept. Please check your design formula.")
        for col in design_cols[1:]:
            # formulaic treatment coding, e.g. "condition[T.B]"
            if '[T.' in col:
                factor_name = col.split('[T.')[0]
                level = col.split('[T.')[1].rstrip(']')
                if factor_name in self.metadata_df.columns:
                    unique_vals = sorted([str(v) for v in self.metadata_df[factor_name].unique()])
                    if len(unique_vals) >= 2:
                        ref_level = unique_vals[0]
                        if level != ref_level:
                            return [factor_name, level, ref_level]
            # a column that names a metadata column outright (e.g. continuous coding)
            elif col in self.metadata_df.columns:
                unique_vals = sorted([str(v) for v in self.metadata_df[col].unique()])
                if len(unique_vals) >= 2:
                    return [col, unique_vals[-1], unique_vals[0]]

        raise ValueError("Cannot determine contrast from design matrix. Please specify --contrast manually.")
    
    def save_results(self, output_dir: str, prefix: str = "deseq2") -> Dict[str, str]:
        """Save analysis results to files.
        
        Args:
            output_dir: Output directory
            prefix: File prefix
            
        Returns:
            Dictionary of saved file paths
        """
        if self.stats_results is None:
            raise ValueError("No results to save. Run analyze() first.")
        
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # name the index column so the CSV header reads gene_id,... like
        # every other exported table
        results = self.stats_results.copy()
        results.index.name = "gene_id"

        saved_files = {}

        # Save full results
        results_file = output_path / f"{prefix}_results.csv"
        results.to_csv(results_file)
        saved_files['results'] = str(results_file)

        # Save significant genes (padj < alpha AND abs(log2FoldChange) > lfc_threshold)
        sig_genes = results[
            (results['padj'] < self.config.alpha) &
            (results['log2FoldChange'].abs() > self.config.lfc_threshold)
        ]
        sig_file = output_path / f"{prefix}_significant.csv"
        sig_genes.to_csv(sig_file)
        saved_files['significant'] = str(sig_file)

        # Save up-regulated genes
        up_genes = sig_genes[sig_genes['log2FoldChange'] > 0]
        up_file = output_path / f"{prefix}_upregulated.csv"
        up_genes.to_csv(up_file)
        saved_files['upregulated'] = str(up_file)

        # Save down-regulated genes
        down_genes = sig_genes[sig_genes['log2FoldChange'] < 0]
        down_file = output_path / f"{prefix}_downregulated.csv"
        down_genes.to_csv(down_file)
        saved_files['downregulated'] = str(down_file)
        
        self.logger.info(f"Results saved to {output_dir}")
        return saved_files
    
    def plot_ma(self, save_path: Optional[str] = None) -> None:
        """Create MA plot from the fitted DeseqStats object.

        Args:
            save_path: Path to save the plot
        """
        if self.stat_res is None:
            raise ValueError("No results to plot. Run analyze() first.")

        try:
            self.stat_res.plot_MA(save_path=save_path)
        except Exception as e:
            self.logger.error(f"Error creating MA plot: {e}")
    
    def plot_volcano(self, save_path: Optional[str] = None) -> None:
        """Create volcano plot.
        
        Args:
            save_path: Path to save the plot
        """
        if self.stats_results is None:
            raise ValueError("No results to plot. Run analyze() first.")
        
        try:
            import matplotlib.pyplot as plt

            # DESeq2 can emit pvalue == 0.0 for extreme effects; -log10(0) = inf
            # breaks the axis and the whole plot is lost to the broad except below
            min_p = np.finfo(float).tiny

            # Create volcano plot
            fig, ax = plt.subplots(figsize=(10, 8))
            try:
                # Plot non-significant genes
                non_sig = self.stats_results[self.stats_results['padj'] >= self.config.alpha]
                ax.scatter(non_sig['log2FoldChange'], -np.log10(non_sig['pvalue'].clip(lower=min_p)),
                          alpha=0.5, label='Non-significant', color='gray', s=10)

                # Plot significant up-regulated genes (same criteria as save_results)
                sig_up = self.stats_results[
                    (self.stats_results['padj'] < self.config.alpha) &
                    (self.stats_results['log2FoldChange'] > self.config.lfc_threshold)
                ]
                ax.scatter(sig_up['log2FoldChange'], -np.log10(sig_up['pvalue'].clip(lower=min_p)),
                          alpha=0.7, label='Up-regulated', color='red', s=20)

                # Plot significant down-regulated genes
                sig_down = self.stats_results[
                    (self.stats_results['padj'] < self.config.alpha) &
                    (self.stats_results['log2FoldChange'] < -self.config.lfc_threshold)
                ]
                ax.scatter(sig_down['log2FoldChange'], -np.log10(sig_down['pvalue'].clip(lower=min_p)),
                          alpha=0.7, label='Down-regulated', color='blue', s=20)

                # Add labels and title
                ax.set_xlabel('log2 Fold Change', fontsize=12)
                ax.set_ylabel('-log10(p-value)', fontsize=12)
                ax.set_title('Volcano Plot', fontsize=14)
                ax.legend(loc='upper right')
                ax.grid(True, alpha=0.3)

                ax.axvline(x=0, color='black', linestyle='-', alpha=0.3)

                plt.tight_layout()

                if save_path:
                    plt.savefig(save_path, dpi=300, bbox_inches='tight')
                    self.logger.info(f"Volcano plot saved to {save_path}")
            finally:
                plt.close(fig)

        except ImportError:
            self.logger.error("Matplotlib is required for plotting. Install with: pip install matplotlib")
        except Exception as e:
            self.logger.error(f"Error creating volcano plot: {e}")
    
    def plot_pca(self, save_path: Optional[str] = None, n_top_genes: int = 500) -> None:
        """Create PCA plot from normalized counts.
        
        Args:
            save_path: Path to save the plot
            n_top_genes: Number of top variable genes to use for PCA
        """
        if self.dds is None:
            raise ValueError("No DeseqDataSet available. Run analyze() first.")
        
        try:
            import matplotlib.pyplot as plt
            from sklearn.decomposition import PCA
            from sklearn.preprocessing import StandardScaler
            
            # Get normalized counts
            if "normed_counts" in self.dds.layers:
                normed_counts = self.dds.layers["normed_counts"]
            else:
                # Use size factor normalized counts
                normed_counts = self.dds.X / self.dds.obs["size_factors"].values[:, None]
            
            # Log transform
            log_counts = np.log1p(normed_counts)
            
            # Select top variable genes
            gene_var = np.var(log_counts, axis=0)
            top_gene_idx = np.argsort(gene_var)[-n_top_genes:]
            log_counts_top = log_counts[:, top_gene_idx]
            
            # Standardize
            scaler = StandardScaler()
            log_counts_scaled = scaler.fit_transform(log_counts_top)
            
            # Run PCA
            pca = PCA(n_components=2)
            pca_result = pca.fit_transform(log_counts_scaled)
            
            # Get condition labels: prefer the last design factor (the factor
            # of interest by convention), then 'condition', else unlabeled
            design_factor = design_columns(self.config.design)
            color_column = design_factor[-1] if design_factor else 'condition'
            if color_column not in self.dds.obs.columns:
                color_column = 'condition' if 'condition' in self.dds.obs.columns else None
            if color_column:
                conditions = self.dds.obs[color_column].values
            else:
                conditions = ['Unknown'] * len(pca_result)
            
            # Create PCA plot
            fig, ax = plt.subplots(figsize=(10, 8))
            try:
                # Plot each condition with different color
                unique_conditions = np.unique(conditions)
                colors = plt.cm.tab10(np.linspace(0, 1, len(unique_conditions)))

                for i, condition in enumerate(unique_conditions):
                    mask = conditions == condition
                    ax.scatter(pca_result[mask, 0], pca_result[mask, 1],
                              label=condition, color=colors[i], s=100, alpha=0.7)

                # Add sample labels
                for i, sample_name in enumerate(self.dds.obs_names):
                    ax.annotate(sample_name, (pca_result[i, 0], pca_result[i, 1]),
                               fontsize=8, alpha=0.7)

                # Add labels and title
                ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=12)
                ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=12)
                ax.set_title('PCA Plot (Top Variable Genes)', fontsize=14)
                ax.legend(loc='best')
                ax.grid(True, alpha=0.3)

                plt.tight_layout()

                if save_path:
                    plt.savefig(save_path, dpi=300, bbox_inches='tight')
                    self.logger.info(f"PCA plot saved to {save_path}")
            finally:
                plt.close(fig)

        except ImportError as e:
            self.logger.error(f"Required package not installed: {e}. Install with: pip install matplotlib scikit-learn")
        except Exception as e:
            self.logger.error(f"Error creating PCA plot: {e}")
    
    def get_summary(self) -> Dict:
        """Get summary statistics of the analysis.
        
        Returns:
            Dictionary with summary statistics
        """
        if self.stats_results is None:
            raise ValueError("No results available. Run analyze() first.")
        
        total_genes = len(self.stats_results)
        
        # Significant genes: padj < alpha AND abs(log2FoldChange) > lfc_threshold
        sig_mask = (
            (self.stats_results['padj'] < self.config.alpha) & 
            (self.stats_results['log2FoldChange'].abs() > self.config.lfc_threshold)
        )
        sig_genes = len(self.stats_results[sig_mask])
        up_genes = len(self.stats_results[sig_mask & (self.stats_results['log2FoldChange'] > 0)])
        down_genes = len(self.stats_results[sig_mask & (self.stats_results['log2FoldChange'] < 0)])
        
        return {
            'total_genes': total_genes,
            'significant_genes': sig_genes,
            'upregulated_genes': up_genes,
            'downregulated_genes': down_genes,
            'alpha': self.config.alpha,
            'lfc_threshold': self.config.lfc_threshold,
            'prefilter_min_count': self.config.prefilter_min_count,
            'design': self.config.design
        }


def run_deseq2_cli(args):
    """Run DESeq2 analysis from CLI arguments.
    
    Args:
        args: Parsed argument namespace
    """
    from rskit.config import DESeq2Config
    
    # Create config
    config = DESeq2Config(
        design=args.design,
        alpha=args.alpha,
        lfc_threshold=args.lfc_threshold,
        prefilter_min_count=args.prefilter_min_count,
        n_cpus=args.threads
    )
    
    # Create analyzer
    analyzer = Deseq2Analyzer(config)
    
    # Setup output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load metadata (coldata)
    logger.info(f"Loading metadata from {args.coldata}")
    metadata_df = analyzer.load_metadata(args.coldata, required_columns=design_columns(args.design))
    contrasts = parse_contrasts(args.contrast, metadata_df)
    
    # Load counts
    if args.salmon_dir:
        existing_counts = SalmonExpressionExporter.find_existing_gene_counts(args.salmon_dir)
        if existing_counts is not None:
            if args.gtf or args.tx2gene:
                logger.warning(
                    f"Reusing precomputed gene counts from {existing_counts}; the given "
                    "--gtf/--tx2gene are ignored. Delete that file to force a re-export."
                )
            logger.info(f"Using precomputed gene counts from {existing_counts}")
            counts_file = str(existing_counts)
        else:
            logger.info(f"Exporting gene-level quantification tables before DESeq2: {args.salmon_dir}")
            expression_outputs = merge_salmon_quant_tables(
                salmon_dir=args.salmon_dir,
                output_dir=args.salmon_dir,
                gtf_file=args.gtf,
                tx2gene=args.tx2gene,
            )
            counts_file = expression_outputs["gene_counts"]
        
    elif args.gene_counts:
        # Load from gene counts file
        logger.info(f"Loading gene counts from {args.gene_counts}")
        counts_file = args.gene_counts
    else:
        raise ValueError("Either --salmon-dir or --gene-counts must be provided")

    counts_df = analyzer.load_counts_from_file(counts_file, metadata_df=metadata_df)

    # Fit the model once; every contrast reuses this fit (the expensive step)
    logger.info("Running DESeq2 analysis...")
    analyzer.fit(
        counts_df=counts_df,
        contrast_factors=[c[0] for c in contrasts] if contrasts else ["condition"],
    )
    if not contrasts:
        contrasts = [analyzer._infer_contrast()]

    multi = len(contrasts) > 1
    if multi:
        logger.info(f"Running {len(contrasts)} contrasts against a single model fit")

    contrast_entries = []
    combined_significant = []
    pca_plot = output_dir / "deseq2_pca_plot.pdf"
    for contrast in contrasts:
        label = f"{contrast[0]}_{contrast[1]}_vs_{contrast[2]}"
        target_dir = output_dir / label if multi else output_dir
        target_dir.mkdir(parents=True, exist_ok=True)

        analyzer.contrast_results(contrast)
        saved_files = analyzer.save_results(str(target_dir))
        logger.info(f"[{label}] Saved results: {saved_files}")

        plot_paths = {
            "volcano_plot": target_dir / "deseq2_volcano_plot.pdf",
            "ma_plot": target_dir / "deseq2_ma_plot.pdf",
        }
        analyzer.plot_volcano(str(plot_paths["volcano_plot"]))
        analyzer.plot_ma(str(plot_paths["ma_plot"]))
        if not pca_plot.exists():
            # PCA is design-level: identical for every contrast, so plot once
            analyzer.plot_pca(str(pca_plot))
        # plot helpers log-and-continue on failure; only record files that exist
        written_plots = {name: str(path) for name, path in plot_paths.items() if path.exists()}
        if pca_plot.exists():
            written_plots["pca_plot"] = str(pca_plot)

        summary = analyzer.get_summary()
        logger.info("\n" + "="*50)
        logger.info(f"DESeq2 Analysis Summary [{label}]")
        logger.info("="*50)
        logger.info(f"Total genes: {summary['total_genes']}")
        logger.info(
            f"Significant genes (padj < {summary['alpha']} and |log2FC| > {summary['lfc_threshold']}): "
            f"{summary['significant_genes']}"
        )
        logger.info(f"  - Up-regulated: {summary['upregulated_genes']}")
        logger.info(f"  - Down-regulated: {summary['downregulated_genes']}")
        logger.info("="*50)

        contrast_entries.append({
            "contrast": contrast,
            "label": label,
            "summary": summary,
            "outputs": {**saved_files, **written_plots},
        })

        if multi:
            significant = pd.read_csv(saved_files["significant"], index_col=0)
            significant.insert(0, "contrast", label)
            combined_significant.append(significant)

    if multi:
        combined_path = output_dir / "deseq2_significant_all.csv"
        pd.concat(combined_significant).to_csv(combined_path)
        logger.info(f"Combined significant genes written to {combined_path}")
        manifest_outputs = {
            "significant_all": str(combined_path),
            "pca_plot": str(pca_plot),
            "contrasts": [entry["label"] for entry in contrast_entries],
        }
    else:
        manifest_outputs = contrast_entries[0]["outputs"]

    manifest = {
        "command": "deseq2",
        "inputs": {
            "coldata": args.coldata,
            "gene_counts": args.gene_counts,
            "salmon_dir": args.salmon_dir,
            "gtf": args.gtf,
            "tx2gene": args.tx2gene,
            "counts_source": "salmon_dir" if args.salmon_dir else "gene_counts",
            "counts_file": counts_file,
        },
        "samples": list(metadata_df.index),
        "design": args.design,
        "contrasts": contrast_entries,
        "outputs": manifest_outputs,
    }
    if not multi:
        # keep the single-contrast manifest shape stable for existing consumers
        manifest["contrast"] = contrasts[0]
        manifest["summary"] = contrast_entries[0]["summary"]

    manifest_path = write_manifest(output_dir, manifest)
    logger.info(f"Saved manifest: {manifest_path}")
    
    return analyzer
