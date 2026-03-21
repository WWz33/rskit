from pathlib import Path
from typing import Dict, List, Optional
from rskit.config import PipelineConfig
from rskit.core.star import StarIndexer, StarAligner
from rskit.core.salmon import SalmonQuantifier
from rskit.core.deseq2 import Deseq2Analyzer
from rskit.utils.logger import get_logger
from rskit.utils.validators import check_star_index

logger = get_logger(__name__)

class RNAseqPipeline:
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.indexer = StarIndexer(config.star)
        self.aligner = StarAligner(config.star)
        self.quantifier = SalmonQuantifier(config.salmon)
        self.deseq2_analyzer = Deseq2Analyzer(config.deseq2)
        self.logger = logger

    def run(self, samples: Dict[str, Dict], genome_fasta: str, gtf_file: str,
            transcript_fasta: str, index_dir: str, output_dir: str,
            quant_output_dir: str, force_index: bool = False) -> Dict:
        results = {}
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        quant_path = Path(quant_output_dir)
        quant_path.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"Checking STAR index at {index_dir}")
        if not Path(index_dir).exists() or not check_star_index(index_dir):
            self.logger.info("Index not found, building...")
            self.indexer.build_index(genome_fasta, gtf_file, index_dir, force=force_index)
        else:
            self.logger.info("Index found, skipping build")

        for sample_name, sample_data in samples.items():
            self.logger.info(f"Processing sample: {sample_name}")
            sample_output = output_path / sample_name
            sample_output.mkdir(parents=True, exist_ok=True)

            align_prefix = str(sample_output / f"{sample_name}_")
            align_results = self.aligner.align(index_dir, sample_data["fq1"], sample_data["fq2"],
                                              align_prefix, sample_name=sample_name)

            salmon_output = quant_path / sample_name
            quant_results = self.quantifier.quantify(transcript_fasta, align_results["transcriptome_bam"],
                                                     str(salmon_output), sample_name=sample_name)

            results[sample_name] = {"alignment": align_results, "quantification": quant_results}
            self.logger.info(f"Completed sample: {sample_name}")

        return results
    
    def run_with_deseq2(self, samples: Dict[str, Dict], genome_fasta: str, gtf_file: str,
                       transcript_fasta: str, index_dir: str, output_dir: str,
                       quant_output_dir: str, metadata: Dict[str, str],
                       contrast: Optional[List[str]] = None,
                       force_index: bool = False) -> Dict:
        """Run pipeline with DESeq2 differential expression analysis.
        
        Args:
            samples: Dictionary of sample data with fq1 and fq2 paths
            genome_fasta: Path to genome FASTA file
            gtf_file: Path to GTF annotation file
            transcript_fasta: Path to transcript FASTA file
            index_dir: Directory for STAR index
            output_dir: Directory for alignment output
            quant_output_dir: Directory for quantification output
            metadata: Dictionary mapping sample names to conditions
            contrast: Contrast for DESeq2 analysis ['condition', 'B', 'A']
            force_index: Whether to force rebuild STAR index
            
        Returns:
            Dictionary with pipeline results and DESeq2 analysis
        """
        # Run standard pipeline
        results = self.run(samples, genome_fasta, gtf_file, transcript_fasta,
                          index_dir, output_dir, quant_output_dir, force_index)
        
        # Prepare data for DESeq2 analysis
        sample_names = list(samples.keys())
        
        # Create counts DataFrame from Salmon results
        import pandas as pd
        counts_data = {}
        for sample_name in sample_names:
            quant_file = Path(quant_output_dir) / sample_name / "quant.sf"
            if quant_file.exists():
                quant_df = pd.read_csv(quant_file, sep='\t')
                gene_counts = quant_df.groupby('Name')['NumReads'].sum()
                counts_data[sample_name] = gene_counts
        
        counts_df = pd.DataFrame(counts_data).T.fillna(0).astype(int)
        
        # Create metadata DataFrame
        metadata_df = pd.DataFrame({
            'sample': sample_names,
            'condition': [metadata.get(sample_name, 'unknown') for sample_name in sample_names]
        })
        metadata_df = metadata_df.set_index('sample')
        
        # Run DESeq2 analysis
        deseq2_output_dir = Path(output_dir) / "deseq2"
        deseq2_output_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            de_results = self.deseq2_analyzer.analyze(counts_df, metadata_df, contrast)
            
            # Save results
            saved_files = self.deseq2_analyzer.save_results(str(deseq2_output_dir))
            
            # Get summary
            summary = self.deseq2_analyzer.get_summary()
            
            results['deseq2'] = {
                'results': de_results,
                'saved_files': saved_files,
                'summary': summary
            }
            
            self.logger.info(f"DESeq2 analysis completed. {summary['significant_genes']} significant genes found.")
            
        except Exception as e:
            self.logger.error(f"DESeq2 analysis failed: {e}")
            results['deseq2'] = {'error': str(e)}
        
        return results
