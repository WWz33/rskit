from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import pandas as pd
import numpy as np
from rskit.config import DESeq2Config
from rskit.utils.logger import get_logger

logger = get_logger(__name__)

class Deseq2Analyzer:
    def __init__(self, config: DESeq2Config):
        self.config = config
        self.logger = logger
        self.dds = None
        self.stats_results = None
        self.counts_df = None
        self.metadata_df = None
        
    def _create_tx2gene_from_gtf(self, gtf_file: str, output_dir: Optional[str] = None) -> pd.DataFrame:
        """Create transcript-to-gene mapping from GTF/GFF3 file.
        
        Args:
            gtf_file: Path to GTF/GFF3 annotation file
            output_dir: If provided, save tx2gene.tsv to this directory
            
        Returns:
            DataFrame with transcript_id and gene_id columns
        """
        from rskit.utils.gtf import open as gtf_open
        
        self.logger.info(f"Parsing GTF/GFF3 file: {gtf_file}")
        
        # Use dict to automatically deduplicate by transcript_id
        tx2gene_map = {}
        num_records = 0
        
        with open(gtf_file, 'r', encoding='utf-8', errors='ignore') as reader:
            for rec in gtf_open(reader, 'ensembl'):
                num_records += 1
                if rec.feature == 'transcript' and rec.transcript_id and rec.gene_id:
                    if rec.transcript_id not in tx2gene_map:
                        tx2gene_map[rec.transcript_id] = rec.gene_id
        
        self.logger.info(f"Scanned {num_records} GTF/GFF3 lines")
        self.logger.info(f"Extracted {len(tx2gene_map)} unique transcript-to-gene mappings")
        
        # Convert to DataFrame
        tx2gene_df = pd.DataFrame(
            [(tx, gene) for tx, gene in tx2gene_map.items()],
            columns=['transcript_id', 'gene_id']
        )
        
        # Save to file for debugging
        if output_dir:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            tx2gene_file = output_path / "tx2gene.tsv"
            tx2gene_df.to_csv(tx2gene_file, sep='\t', index=False)
            self.logger.info(f"Saved tx2gene mapping to: {tx2gene_file}")
        
        self.logger.info(f"tx2gene preview:\n{tx2gene_df.head()}")
        return tx2gene_df
    
    def load_counts_from_salmon(self, salmon_dir: str, coldata: pd.DataFrame, 
                                 gtf_file: Optional[str] = None,
                                 tx2gene: Optional[str] = None,
                                 output_dir: Optional[str] = None) -> pd.DataFrame:
        """Load count data from Salmon quantification results using pytximport.
        
        Args:
            salmon_dir: Directory containing Salmon quantification results
            coldata: DataFrame with sample information (index=sample names)
            gtf_file: Path to GTF annotation file (to create tx2gene map)
            tx2gene: Path to tx2gene mapping file (CSV or TSV)
            output_dir: Directory to save tx2gene.tsv for debugging
            
        Returns:
            DataFrame with gene-level counts (samples x genes)
        """
        try:
            from pytximport import tximport
        except ImportError:
            raise ImportError("pytximport is not installed. Please install it with: pip install pytximport")
        
        salmon_path = Path(salmon_dir)
        
        # Get list of quant.sf files
        file_paths = []
        sample_names = []
        for sample_name in coldata.index:
            quant_file = salmon_path / sample_name / "quant.sf"
            if quant_file.exists():
                file_paths.append(str(quant_file))
                sample_names.append(sample_name)
            else:
                self.logger.warning(f"Quantification file not found for sample {sample_name}")
        
        if not file_paths:
            raise FileNotFoundError(f"No quant.sf files found in {salmon_dir}")
        
        # Get transcript-to-gene mapping
        if tx2gene is not None:
            sep = '\t' if tx2gene.endswith('.tsv') else ','
            tx2gene_map = pd.read_csv(tx2gene, sep=sep)
            self.logger.info(f"Loaded tx2gene map from {tx2gene}")
        elif gtf_file is not None:
            # Use local GTF parsing instead of BioMart
            tx2gene_map = self._create_tx2gene_from_gtf(gtf_file, output_dir)
        else:
            raise ValueError("Either gtf_file or tx2gene must be provided")
        
        # Ensure correct column names and save if loaded from file
        if 'transcript_id' not in tx2gene_map.columns or 'gene_id' not in tx2gene_map.columns:
            if len(tx2gene_map.columns) >= 2:
                tx2gene_map = tx2gene_map.iloc[:, :2]
                tx2gene_map.columns = ['transcript_id', 'gene_id']
                self.logger.warning("Renamed tx2gene columns to transcript_id and gene_id")
            else:
                raise ValueError("tx2gene map must have at least 2 columns (transcript_id, gene_id)")
        
        # Save tx2gene for debugging if loaded from file (GTF path already saves)
        if tx2gene is not None and output_dir:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            tx2gene_file = output_path / "tx2gene.tsv"
            tx2gene_map.to_csv(tx2gene_file, sep='\t', index=False)
            self.logger.info(f"Saved tx2gene mapping to {tx2gene_file}")
        
        # Final stats before tximport
        self.logger.info(f"Ready for tximport: {len(tx2gene_map)} transcripts, {tx2gene_map['gene_id'].nunique()} genes")
        
        # Run tximport
        self.logger.info(f"Running pytximport on {len(file_paths)} samples...")
        try:
            results = tximport(
                file_paths=file_paths,
                data_type="salmon",
                transcript_gene_map=tx2gene_map,
                counts_from_abundance="length_scaled_tpm",
                ignore_transcript_version=False,
                ignore_after_bar=False,
                output_type="dict"
            )
        except AssertionError as e:
            self.logger.error(f"pytximport failed: {e}")
            self.logger.error(f"tx2gene map has {len(tx2gene_map)} entries")
            # Check for duplicates in tx2gene
            dup_count = tx2gene_map['transcript_id'].duplicated().sum()
            self.logger.error(f"Duplicate transcript_ids in tx2gene: {dup_count}")
            if output_dir:
                # Save duplicate info
                dups = tx2gene_map[tx2gene_map['transcript_id'].duplicated(keep=False)]
                dup_file = Path(output_dir) / "tx2gene_duplicates.tsv"
                dups.to_csv(dup_file, sep='\t', index=False)
                self.logger.error(f"Saved duplicates to {dup_file}")
            raise
        
        # Extract counts matrix - pytximport returns xarray DataArray
        counts_matrix = results["counts"]
        
        # Convert xarray DataArray to pandas DataFrame
        # DataArray has dims (sample, transcript/gene), need to transpose
        counts_df = counts_matrix.to_pandas().T
        counts_df.index.name = None
        counts_df.columns.name = None
        
        # Ensure sample names match
        if list(counts_df.columns) != sample_names:
            counts_df = counts_df[sample_names]
        
        counts_df = counts_df.round().astype(int)
        
        self.counts_df = counts_df
        self.logger.info(f"Loaded counts for {counts_df.shape[0]} samples and {counts_df.shape[1]} genes")
        
        return counts_df
    
    def load_counts_from_file(self, counts_file: str) -> pd.DataFrame:
        """Load count data from file.
        
        Args:
            counts_file: Path to counts matrix file (genes x samples or samples x genes)
            
        Returns:
            DataFrame with counts (samples x genes)
        """
        counts_df = pd.read_csv(counts_file, index_col=0)
        
        # Check if we need to transpose (genes x samples -> samples x genes)
        # Heuristic: if columns are mostly numeric-looking sample names, transpose
        if counts_df.shape[0] > counts_df.shape[1] * 2:
            self.logger.info("Transposing counts matrix to samples x genes format")
            counts_df = counts_df.T
        
        # Ensure integer counts
        counts_df = counts_df.round().astype(int)
        
        self.counts_df = counts_df
        self.logger.info(f"Loaded counts for {counts_df.shape[0]} samples and {counts_df.shape[1]} genes")
        
        return counts_df
    
    def load_metadata(self, metadata_file: str) -> pd.DataFrame:
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
        sep = '\t' if metadata_file.endswith('.tsv') else ','
        
        # Read the file
        metadata_df = pd.read_csv(metadata_file, sep=sep)
        
        # Check if 'sample' column exists, if so use it as index
        if 'sample' in metadata_df.columns:
            metadata_df = metadata_df.set_index('sample')
        elif metadata_df.columns[0] != metadata_df.index.name:
            # Use first column as index if no 'sample' column
            metadata_df = metadata_df.set_index(metadata_df.columns[0])
        
        self.metadata_df = metadata_df
        self.logger.info(f"Loaded metadata for {metadata_df.shape[0]} samples")
        self.logger.info(f"Columns: {list(metadata_df.columns)}")
        self.logger.info(f"Conditions: {dict(metadata_df['condition'].value_counts()) if 'condition' in metadata_df.columns else 'N/A'}")
        
        return metadata_df
    
    def analyze(self, counts_df: Optional[pd.DataFrame] = None, 
                metadata_df: Optional[pd.DataFrame] = None,
                contrast: Optional[List[str]] = None) -> pd.DataFrame:
        """Perform differential expression analysis using PyDESeq2.
        
        Args:
            counts_df: Count matrix (samples x genes), or use self.counts_df
            metadata_df: Sample metadata, or use self.metadata_df
            contrast: Contrast specification ['condition', 'B', 'A']
            
        Returns:
            DataFrame with differential expression results
        """
        try:
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats
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
        
        # Ensure sample names match
        common_samples = self.counts_df.index.intersection(self.metadata_df.index)
        if len(common_samples) == 0:
            raise ValueError("No common samples between counts and metadata")
        if len(common_samples) < len(self.counts_df.index):
            self.logger.warning(f"Only {len(common_samples)} samples match between counts and metadata")
            self.counts_df = self.counts_df.loc[common_samples]
            self.metadata_df = self.metadata_df.loc[common_samples]
        
        # Initialize inference
        inference = DefaultInference(n_cpus=self.config.n_cpus)
        
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
            inference=inference,
            quiet=self.config.quiet
        )
        
        # Run DESeq2 pipeline
        self.logger.info("Running DESeq2 pipeline...")
        self.dds.deseq2()
        
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
            inference=inference,
            quiet=self.config.quiet
        )
        
        # Run Wald test
        stat_res.summary()
        
        # Store results
        self.stats_results = stat_res.results_df
        
        # Apply LFC shrinkage
        try:
            coeff = f"{contrast[0]}_{contrast[1]}_vs_{contrast[2]}"
            stat_res.lfc_shrink(coeff=coeff)
            self.logger.info(f"LFC shrinkage applied for coefficient: {coeff}")
        except Exception as e:
            self.logger.warning(f"Could not apply LFC shrinkage: {e}")
        
        return self.stats_results
    
    def _infer_contrast(self) -> List[str]:
        """Infer contrast from design matrix and metadata.
        
        For multi-factor designs, tries to find the main factor of interest.
        Typically looks for 'condition' column or the column with most levels.
        """
        # First, try to use 'condition' column if it exists
        if 'condition' in self.metadata_df.columns:
            unique_vals = self.metadata_df['condition'].unique()
            if len(unique_vals) >= 2:
                # Sort to ensure consistent ordering
                sorted_vals = sorted([str(v) for v in unique_vals])
                self.logger.info(f"Using 'condition' column for contrast: {sorted_vals[-1]} vs {sorted_vals[0]}")
                return ['condition', sorted_vals[-1], sorted_vals[0]]
        
        # If no 'condition' column, try to infer from design matrix
        design_cols = self.dds.obsm["design_matrix"].columns
        if len(design_cols) > 1:
            # Get the first non-intercept column
            for col in design_cols[1:]:
                # Try different parsing strategies for column names
                # Strategy 1: column format like "condition[T.B]" or "condition[T.lhy-D]"
                if '[T.' in col:
                    factor_name = col.split('[T.')[0]
                    # Extract level from "factor[T.level]"
                    level = col.split('[T.')[1].rstrip(']')
                    if factor_name in self.metadata_df.columns:
                        unique_vals = sorted([str(v) for v in self.metadata_df[factor_name].unique()])
                        if len(unique_vals) >= 2:
                            ref_level = unique_vals[0]
                            if level != ref_level:
                                return [factor_name, level, ref_level]
                
                # Strategy 2: column format like "factor_level"
                parts = col.split('_')
                if len(parts) >= 2:
                    factor_name = parts[0]
                    if factor_name in self.metadata_df.columns:
                        unique_vals = sorted([str(v) for v in self.metadata_df[factor_name].unique()])
                        if len(unique_vals) >= 2:
                            return [factor_name, unique_vals[-1], unique_vals[0]]
            
            raise ValueError("Cannot determine contrast from design matrix. Please specify --contrast manually.")
        else:
            raise ValueError("Design matrix has only intercept. Please check your design formula.")
    
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
        
        saved_files = {}
        
        # Save full results
        results_file = output_path / f"{prefix}_results.csv"
        self.stats_results.to_csv(results_file)
        saved_files['results'] = str(results_file)
        
        # Save significant genes
        sig_genes = self.stats_results[self.stats_results['padj'] < self.config.alpha]
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
    
    def save_gene_counts(self, output_dir: str, prefix: str = "gene_counts") -> str:
        """Save gene counts matrix.
        
        Args:
            output_dir: Output directory
            prefix: File prefix
            
        Returns:
            Path to saved file
        """
        if self.counts_df is None:
            raise ValueError("No counts data to save.")
        
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        counts_file = output_path / f"{prefix}.csv"
        self.counts_df.to_csv(counts_file)
        
        self.logger.info(f"Gene counts saved to {counts_file}")
        return str(counts_file)
    
    def plot_ma(self, save_path: Optional[str] = None) -> None:
        """Create MA plot.
        
        Args:
            save_path: Path to save the plot
        """
        if self.dds is None or self.stats_results is None:
            raise ValueError("No results to plot. Run analyze() first.")
        
        try:
            from pydeseq2.ds import DeseqStats
            stat_res = DeseqStats(
                dds=self.dds,
                contrast=self.stats_results.columns[0],
                alpha=self.config.alpha
            )
            stat_res.results_df = self.stats_results
            stat_res.plot_MA(save_path=save_path)
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
            
            # Create volcano plot
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # Plot non-significant genes
            non_sig = self.stats_results[self.stats_results['padj'] >= self.config.alpha]
            ax.scatter(non_sig['log2FoldChange'], -np.log10(non_sig['pvalue']), 
                      alpha=0.5, label='Non-significant', color='gray', s=10)
            
            # Plot significant up-regulated genes
            sig_up = self.stats_results[
                (self.stats_results['padj'] < self.config.alpha) & 
                (self.stats_results['log2FoldChange'] > 0)
            ]
            ax.scatter(sig_up['log2FoldChange'], -np.log10(sig_up['pvalue']), 
                      alpha=0.7, label='Up-regulated', color='red', s=20)
            
            # Plot significant down-regulated genes
            sig_down = self.stats_results[
                (self.stats_results['padj'] < self.config.alpha) & 
                (self.stats_results['log2FoldChange'] < 0)
            ]
            ax.scatter(sig_down['log2FoldChange'], -np.log10(sig_down['pvalue']), 
                      alpha=0.7, label='Down-regulated', color='blue', s=20)
            
            # Add labels and title
            ax.set_xlabel('log2 Fold Change', fontsize=12)
            ax.set_ylabel('-log10(p-value)', fontsize=12)
            ax.set_title('Volcano Plot', fontsize=14)
            ax.legend(loc='upper right')
            ax.grid(True, alpha=0.3)
            
            # Add threshold lines
            ax.axhline(y=-np.log10(self.config.alpha), color='blue', linestyle='--', alpha=0.5, label=f'p={self.config.alpha}')
            ax.axvline(x=0, color='black', linestyle='-', alpha=0.3)
            
            plt.tight_layout()
            
            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches='tight')
                self.logger.info(f"Volcano plot saved to {save_path}")
            
            plt.close()
            
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
            
            # Get condition labels
            if hasattr(self.dds, 'obs') and 'condition' in self.dds.obs.columns:
                conditions = self.dds.obs['condition'].values
            else:
                conditions = ['Unknown'] * len(pca_result)
            
            # Create PCA plot
            fig, ax = plt.subplots(figsize=(10, 8))
            
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
            
            plt.close()
            
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
        sig_genes = len(self.stats_results[self.stats_results['padj'] < self.config.alpha])
        up_genes = len(self.stats_results[
            (self.stats_results['padj'] < self.config.alpha) & 
            (self.stats_results['log2FoldChange'] > 0)
        ])
        down_genes = len(self.stats_results[
            (self.stats_results['padj'] < self.config.alpha) & 
            (self.stats_results['log2FoldChange'] < 0)
        ])
        
        return {
            'total_genes': total_genes,
            'significant_genes': sig_genes,
            'upregulated_genes': up_genes,
            'downregulated_genes': down_genes,
            'alpha': self.config.alpha,
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
        n_cpus=args.threads
    )
    
    # Create analyzer
    analyzer = Deseq2Analyzer(config)
    
    # Setup output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load metadata (coldata)
    logger.info(f"Loading metadata from {args.coldata}")
    metadata_df = analyzer.load_metadata(args.coldata)
    
    # Load counts
    if args.salmon_dir:
        # Load from Salmon using pytximport
        logger.info(f"Loading counts from Salmon directory: {args.salmon_dir}")
        counts_df = analyzer.load_counts_from_salmon(
            salmon_dir=args.salmon_dir,
            coldata=metadata_df,
            gtf_file=args.gtf,
            tx2gene=args.tx2gene,
            output_dir=str(output_dir)
        )
        
        # Save tximport-processed gene counts
        counts_file = analyzer.save_gene_counts(str(output_dir), "gene_counts_tximport")
        logger.info(f"Saved tximport gene counts to {counts_file}")
        
    elif args.gene_counts:
        # Load from gene counts file
        logger.info(f"Loading gene counts from {args.gene_counts}")
        counts_df = analyzer.load_counts_from_file(args.gene_counts)
    else:
        raise ValueError("Either --salmon-dir or --gene-counts must be provided")
    
    # Parse contrast if provided
    contrast = None
    if args.contrast:
        contrast = args.contrast.split(',')
        if len(contrast) != 3:
            raise ValueError("Contrast must be in format: factor,level1,level2 (e.g., condition,B,A)")
    
    # Run analysis
    logger.info("Running DESeq2 analysis...")
    results = analyzer.analyze(contrast=contrast)
    
    # Save results
    saved_files = analyzer.save_results(str(output_dir))
    logger.info(f"Saved results: {saved_files}")
    
    # Generate plots
    logger.info("Generating plots...")
    
    # Volcano plot
    volcano_path = output_dir / "volcano_plot.pdf"
    analyzer.plot_volcano(str(volcano_path))
    
    # PCA plot
    pca_path = output_dir / "pca_plot.pdf"
    analyzer.plot_pca(str(pca_path))
    
    # Print summary
    summary = analyzer.get_summary()
    logger.info("\n" + "="*50)
    logger.info("DESeq2 Analysis Summary")
    logger.info("="*50)
    logger.info(f"Total genes: {summary['total_genes']}")
    logger.info(f"Significant genes (padj < {summary['alpha']}): {summary['significant_genes']}")
    logger.info(f"  - Up-regulated: {summary['upregulated_genes']}")
    logger.info(f"  - Down-regulated: {summary['downregulated_genes']}")
    logger.info("="*50)
    
    return analyzer
