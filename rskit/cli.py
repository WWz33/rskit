import argparse
import sys
import os
import pandas as pd
import subprocess
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, List, Dict, Tuple
from rskit.config import StarConfig, SalmonConfig, PipelineConfig, DESeq2Config
from rskit.core.pipeline import RNAseqPipeline
from rskit.core.deseq2 import Deseq2Analyzer, run_deseq2_cli
from rskit.core.wgcna import run_wgcna_cli
from rskit.utils.logger import get_logger
from rskit.utils.validators import check_and_prepare_index
from rskit.utils.parallel import calculate_threads_per_sample, run_samples_parallel
from rskit.core.star import StarIndexer
from rskit.core.salmon import SalmonExpressionExporter, SalmonQuantifier

logger = get_logger(__name__)


@dataclass
class Deseq2Args:
    """Arguments for DESeq2 analysis"""
    salmon_dir: Optional[str]
    coldata: str
    gtf: Optional[str]
    output_dir: str
    gene_counts: Optional[str] = None
    tx2gene: Optional[str] = None
    work_dir: str = "."
    design: str = "~condition"
    contrast: Optional[List[str]] = None
    alpha: float = 0.05
    lfc_threshold: float = 2.0
    threads: Optional[int] = None


def setup_workdir(output_dir: str) -> Dict[str, Path]:
    """Setup work directory structure"""
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    
    dirs = {
        'index': output_path / '00_index',
        'clean_data': output_path / '01_clean_data',
        'clean_data_html': output_path / '01_clean_data' / 'html',
        'clean_data_json': output_path / '01_clean_data' / 'json',
        'bam': output_path / '02_bam',
        'quant': output_path / '03_quant',
        'deseq2': output_path / '04_deseq2'
    }
    
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    
    return dirs


def trim_reads(read1: Path, read2: Path, sample: str, workdirs: Dict[str, Path], threads: int = 8) -> Tuple[str, str]:
    """Run fastp trimming"""
    prefix = workdirs['clean_data'] / sample
    cmd = [
        'fastp',
        '-i', str(read1),
        '-I', str(read2),
        '-o', str(prefix) + '_1.fq',
        '-O', str(prefix) + '_2.fq',
        '--thread', str(threads),
        '-j', str(workdirs['clean_data_json'] / sample) + '.json',
        '-h', str(workdirs['clean_data_html'] / sample) + '.html'
    ]
    logger.info(f"[{sample}] Trimming...")
    subprocess.run(cmd, check=True, capture_output=True)
    logger.info(f"[{sample}] Trimming completed")
    return str(prefix) + '_1.fq', str(prefix) + '_2.fq'


def parse_samples_from_coldata(coldata: str):
    """Parse samples from coldata file"""
    sep = '\t' if coldata.endswith('.tsv') else ','
    samples_df = pd.read_csv(coldata, sep=sep)
    required_cols = {'sample', 'r1', 'r2'}
    if not required_cols.issubset(samples_df.columns):
        raise ValueError(f"Coldata file must contain columns: {required_cols}")
    return [(row['sample'], Path(row['r1']).resolve(), Path(row['r2']).resolve()) 
            for _, row in samples_df.iterrows()]


def trim_sample_wrapper(args):
    """Wrapper for parallel trimming"""
    sample_name, r1_path, r2_path, workdirs, threads = args
    r1_clean, r2_clean = trim_reads(r1_path, r2_path, sample_name, workdirs, threads)
    return sample_name, r1_clean, r2_clean


def prepare_samples(samples_list, workdirs: Dict[str, Path], 
                    trim: bool, threads: int, parallel: bool = False) -> Dict[str, Dict]:
    """Prepare samples dict with optional trimming"""
    from concurrent.futures import ThreadPoolExecutor, as_completed
    
    if trim:
        if parallel and len(samples_list) > 1:
            # Parallel trimming
            num_samples = len(samples_list)
            logger.info(f"Parallel trimming: {num_samples} samples")
            
            trim_args = [(name, r1, r2, workdirs, threads) for name, r1, r2 in samples_list]
            
            samples = {}
            with ThreadPoolExecutor(max_workers=num_samples) as executor:
                futures = {executor.submit(trim_sample_wrapper, args): args[0] for args in trim_args}
                for i, future in enumerate(as_completed(futures), 1):
                    sample_name, r1_clean, r2_clean = future.result()
                    samples[sample_name] = {'fq1': r1_clean, 'fq2': r2_clean}
                    logger.info(f"Trim progress: {i}/{num_samples} completed")
        else:
            # Sequential trimming
            samples = {}
            for sample_name, r1_path, r2_path in samples_list:
                r1_clean, r2_clean = trim_reads(r1_path, r2_path, sample_name, workdirs, threads)
                samples[sample_name] = {'fq1': r1_clean, 'fq2': r2_clean}
    else:
        samples = {name: {'fq1': str(r1), 'fq2': str(r2)} for name, r1, r2 in samples_list}
    
    return samples


def build_index_if_needed(index_dir: Path, genome_fasta: str, gtf_file: str, 
                          threads: int, force_index: bool) -> bool:
    """Check and build index if needed, returns True if index was built"""
    index_dir, needs_build = check_and_prepare_index(str(index_dir), force_index)
    if needs_build:
        indexer = StarIndexer(StarConfig(threads=threads))
        indexer.build_index(genome_fasta, gtf_file, str(index_dir), force=force_index)
    return needs_build


def run_quantification(samples: Dict[str, Dict], genome_fasta: str, gtf_file: str,
                       transcript_fasta: str, index_dir: Path, workdirs: Dict[str, Path],
                       threads_per_sample: int, parallel: bool, skip_existing: bool = False) -> Dict:
    """Run quantification pipeline (parallel or sequential)"""
    num_samples = len(samples)
    
    if parallel and num_samples > 1:
        return run_samples_parallel(samples, str(index_dir), transcript_fasta, 
                                    workdirs, threads_per_sample, skip_existing)
    else:
        config = PipelineConfig(
            star=StarConfig(threads=threads_per_sample),
            salmon=SalmonConfig(threads=threads_per_sample),
            output_dir=str(workdirs['bam'])
        )
        pipeline = RNAseqPipeline(config)
        return pipeline.run(
            samples=samples,
            genome_fasta=genome_fasta,
            gtf_file=gtf_file,
            transcript_fasta=transcript_fasta,
            index_dir=str(index_dir),
            output_dir=str(workdirs['bam']),
            quant_output_dir=str(workdirs['quant']),
            force_index=False
        )


def export_quant_expression_tables(
    quant_dir: Path,
    gtf_file: Optional[str],
    tx2gene: Optional[str] = None,
    sample_names: Optional[List[str]] = None,
) -> Dict[str, str]:
    """Export gene-level expression tables from Salmon quantification output."""
    exporter = SalmonExpressionExporter()
    outputs = exporter.export_gene_tables(
        salmon_dir=str(quant_dir),
        output_dir=str(quant_dir),
        gtf_file=gtf_file,
        tx2gene=tx2gene,
        sample_names=sample_names,
    )
    logger.info(f"Expression tables written to {quant_dir}")
    return outputs


def main_quant(args):
    """Run quantification pipeline"""
    # Convert paths to absolute
    r1 = Path(args.r1).resolve() if args.r1 else None
    r2 = Path(args.r2).resolve() if args.r2 else None
    genome_fasta = str(Path(args.genome_fasta).resolve())
    gtf_file = str(Path(args.gtf_file).resolve())
    transcript_fasta = str(Path(args.transcript_fasta).resolve())
    
    # Setup work directory
    workdirs = setup_workdir(args.output_dir)
    
    # Determine and check index directory
    index_dir = Path(args.index_dir).resolve() if args.index_dir else workdirs['index']
    build_index_if_needed(index_dir, genome_fasta, gtf_file, args.threads, args.force_index)
    
    # Parse samples
    if args.coldata:
        samples_list = parse_samples_from_coldata(args.coldata)
    else:
        if not all([args.sample, r1, r2]):
            raise ValueError("Must provide --sample, --r1, --r2 or use --coldata")
        samples_list = [(args.sample, r1, r2)]
    
    # Calculate threads per sample
    num_samples = len(samples_list)
    use_parallel = args.parallel is not None and num_samples > 1
    if use_parallel:
        threads_per_sample = calculate_threads_per_sample(args.parallel, num_samples)
        logger.info(f"Parallel: {args.parallel} cores / {num_samples} samples = {threads_per_sample} threads/sample")
    else:
        threads_per_sample = args.threads
    
    # Prepare and run samples
    samples = prepare_samples(samples_list, workdirs, args.trim, threads_per_sample, parallel=use_parallel)
    results = run_quantification(samples, genome_fasta, gtf_file, transcript_fasta,
                                 index_dir, workdirs, threads_per_sample, use_parallel, args.skip_existing)
    expression_outputs = export_quant_expression_tables(
        quant_dir=workdirs['quant'],
        gtf_file=gtf_file,
        tx2gene=getattr(args, "tx2gene", None),
        sample_names=[sample_name for sample_name, _, _ in samples_list],
    )
    
    logger.info(f"Pipeline completed. Processed {len(results)} samples.")
    logger.info(f"Gene-level outputs: {expression_outputs}")


def main_deseq2(args):
    """Run DESeq2 differential expression analysis"""
    logger.info("="*60)
    logger.info("DESeq2 Differential Expression Analysis")
    logger.info("="*60)
    
    # Validate input arguments
    if not args.salmon_dir and not args.gene_counts:
        logger.error("Either --salmon-dir or --gene-counts must be provided")
        sys.exit(1)
    
    if args.salmon_dir and args.gene_counts:
        logger.warning("Both --salmon-dir and --gene-counts provided. Using --salmon-dir with pytximport")
    
    # Setup work directory structure
    workdirs = setup_workdir(args.work_dir)
    
    # Use default 04_deseq2 directory if output_dir not specified
    output_dir = Path(args.output_dir) if args.output_dir else workdirs['deseq2']
    output_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir = str(output_dir)
    
    # Run DESeq2 analysis
    try:
        analyzer = run_deseq2_cli(args)
        logger.info("DESeq2 analysis completed successfully!")
    except Exception as e:
        logger.error(f"DESeq2 analysis failed: {e}")
        raise


def main_wgcna(args):
    """Run WGCNA co-expression network analysis"""
    logger.info("="*60)
    logger.info("WGCNA Co-expression Network Analysis")
    logger.info("="*60)
    
    if not args.expression:
        logger.error("Expression file must be provided")
        sys.exit(1)
    
    try:
        analyzer = run_wgcna_cli(args)
        logger.info("WGCNA analysis completed successfully!")
    except Exception as e:
        logger.error(f"WGCNA analysis failed: {e}")
        raise


def main_all(args):
    """Run complete pipeline: quant -> deseq2"""
    logger.info("="*60)
    logger.info("Complete Pipeline: Quantification + DESeq2")
    logger.info("="*60)
    
    # Convert paths to absolute
    genome_fasta = str(Path(args.genome_fasta).resolve())
    gtf_file = str(Path(args.gtf_file).resolve())
    transcript_fasta = str(Path(args.transcript_fasta).resolve())
    coldata_file = str(Path(args.coldata).resolve())
    
    # Setup work directory
    workdirs = setup_workdir(args.output_dir)
    
    # Determine and check index directory
    index_dir = Path(args.index_dir).resolve() if args.index_dir else workdirs['index']
    build_index_if_needed(index_dir, genome_fasta, gtf_file, args.threads, args.force_index)
    
    # Parse samples from coldata
    samples_list = parse_samples_from_coldata(coldata_file)
    
    # Calculate threads per sample
    num_samples = len(samples_list)
    use_parallel = args.parallel is not None and num_samples > 1
    if use_parallel:
        threads_per_sample = calculate_threads_per_sample(args.parallel, num_samples)
        logger.info(f"Parallel: {args.parallel} cores / {num_samples} samples = {threads_per_sample} threads/sample")
    else:
        threads_per_sample = args.threads
    
    # Step 1: Run quantification
    logger.info("="*60)
    logger.info("Step 1: Quantification Pipeline")
    logger.info("="*60)
    
    samples = prepare_samples(samples_list, workdirs, args.trim, threads_per_sample, parallel=use_parallel)
    results = run_quantification(samples, genome_fasta, gtf_file, transcript_fasta,
                                 index_dir, workdirs, threads_per_sample, use_parallel, args.skip_existing)
    logger.info(f"Quantification completed. Processed {len(results)} samples.")
    expression_outputs = export_quant_expression_tables(
        quant_dir=workdirs['quant'],
        gtf_file=gtf_file,
        tx2gene=args.tx2gene,
        sample_names=[sample_name for sample_name, _, _ in samples_list],
    )
    
    # Step 2: Run DESeq2
    logger.info("="*60)
    logger.info("Step 2: DESeq2 Differential Expression Analysis")
    logger.info("="*60)
    
    deseq2_args = Deseq2Args(
        salmon_dir=None,
        coldata=coldata_file,
        gtf=None,
        output_dir=str(workdirs['deseq2']),
        gene_counts=expression_outputs["gene_counts"],
        tx2gene=args.tx2gene,
        work_dir=args.output_dir,
        design=args.design,
        contrast=args.contrast,
        alpha=args.alpha,
        lfc_threshold=getattr(args, 'lfc_threshold', 2.0),
        threads=args.threads
    )
    
    try:
        analyzer = run_deseq2_cli(deseq2_args)
        logger.info("DESeq2 analysis completed successfully!")
    except Exception as e:
        logger.error(f"DESeq2 analysis failed: {e}")
        raise
    
    logger.info("="*60)
    logger.info("Complete pipeline finished successfully!")
    logger.info(f"Results saved to: {args.output_dir}")
    logger.info("="*60)


def main():
    parser = argparse.ArgumentParser(
        description="RNA-seq analysis toolkit",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparsers = parser.add_subparsers(help="Task to perform")
    
    # quant command
    parser_quant = subparsers.add_parser("quant", help="Complete quantification pipeline (index -> align -> quant)")
    
    parser_quant.add_argument("-s", "--sample", help="Sample name")
    parser_quant.add_argument("-S", "--coldata", help="Sample file (CSV/TSV) with columns: sample,r1,r2")
    parser_quant.add_argument("-1", "--r1", help="First read file")
    parser_quant.add_argument("-2", "--r2", help="Second read file")
    parser_quant.add_argument("-g", "--genome-fasta", dest="genome_fasta", required=True, help="Genome FASTA file")
    parser_quant.add_argument("-gtf", "--gtf-file", dest="gtf_file", required=True, help="GTF annotation file")
    parser_quant.add_argument("-gf", "--transcript-fasta", dest="transcript_fasta", required=True, help="Transcript FASTA file")
    parser_quant.add_argument("-o", "--output-dir", dest="output_dir", required=True, help="Output directory (work directory)")
    parser_quant.add_argument("-idx", "--index-dir", dest="index_dir", help="STAR index directory (default: <output_dir>/00_index)")
    parser_quant.add_argument("-t2g", "--tx2gene", dest="tx2gene",
        help="Path to transcript-to-gene mapping file for gene-level expression export")
    parser_quant.add_argument("-t", "--threads", type=int, default=8, help="Number of threads per sample")
    parser_quant.add_argument("-p", "--parallel", type=int, help="Total cores for parallel processing")
    parser_quant.add_argument("--trim", action="store_true", help="Trim reads with fastp")
    parser_quant.add_argument("--force-index", action="store_true", help="Force rebuild index")
    parser_quant.add_argument("--skip-existing", action="store_true", help="Skip samples if output already exists")
    parser_quant.set_defaults(func=main_quant)
    
    # deseq2 command
    parser_deseq2 = subparsers.add_parser("deseq2", 
        help="DESeq2 differential expression analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
Coldata format (CSV):
    sample,id,condition
    sample1,ctrl,control
    sample2,treat,treatment

Examples:
    rskit deseq2 -sd ./03_quant -S coldata.csv -gtf annotation.gtf
    rskit deseq2 -gc counts.csv -S coldata.csv --design "~batch + condition"
        """)
    
    input_group = parser_deseq2.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-sd", "--salmon-dir", dest="salmon_dir", 
        help="Directory containing Salmon quant subfolders")
    input_group.add_argument("-gc", "--gene-counts", dest="gene_counts",
        help="Path to gene counts matrix file (CSV/TSV)")
    
    parser_deseq2.add_argument("-S", "--coldata", required=True,
        help="Path to coldata/metadata file with columns: sample,id,condition")
    parser_deseq2.add_argument("-gtf", "--gtf", 
        help="GTF annotation file (required when using --salmon-dir without --tx2gene)")
    parser_deseq2.add_argument("-t2g", "--tx2gene", dest="tx2gene",
        help="Path to transcript-to-gene mapping file")
    parser_deseq2.add_argument("-w", "--work-dir", dest="work_dir", default=".",
        help="Work directory")
    parser_deseq2.add_argument("-o", "--output-dir", dest="output_dir", default=None,
        help="Custom output directory for DESeq2 results")
    parser_deseq2.add_argument("--design", default="~condition",
        help="Design formula (e.g., '~condition', '~batch + condition')")
    parser_deseq2.add_argument("--contrast",
        help="Contrast specification (e.g., 'condition,treatment,control')")
    parser_deseq2.add_argument("--alpha", type=float, default=0.05,
        help="Significance threshold for adjusted p-values")
    parser_deseq2.add_argument("--lfc", dest="lfc_threshold", type=float, default=2.0,
        help="Log2 fold change threshold for significant genes")
    parser_deseq2.add_argument("-t", "--threads", type=int, default=None,
        help="Number of threads for parallel processing")
    
    parser_deseq2.set_defaults(func=main_deseq2)
    
    # wgcna command
    parser_wgcna = subparsers.add_parser("wgcna", 
        help="WGCNA co-expression network analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
Examples:
    rskit wgcna -e expression.csv -o ./wgcna_results
    rskit wgcna -e expression.csv -S coldata.csv -G gene_info.csv -o ./wgcna_results
        """)
    
    parser_wgcna.add_argument("-e", "--expression", required=True,
        help="Path to gene expression matrix file (CSV/TSV)")
    parser_wgcna.add_argument("-o", "--output-dir", dest="output_dir", required=True,
        help="Output directory for WGCNA results")
    parser_wgcna.add_argument("-S", "--coldata", dest="coldata",
        help="Path to sample metadata file (CSV/TSV)")
    parser_wgcna.add_argument("-G", "--gene-info", dest="gene_info",
        help="Path to gene metadata file (CSV/TSV)")
    parser_wgcna.add_argument("-sep", "--sep", default=",",
        help="Separator for input files (comma or tab)")
    parser_wgcna.add_argument("-n", "--name", default="WGCNA",
        help="Name for the WGCNA analysis")
    parser_wgcna.add_argument("-s", "--species", 
        help="Species for enrichment analysis (Human, Mouse, Yeast, Fly, Fish, Worm)")
    parser_wgcna.add_argument("-l", "--level", default="gene", choices=["gene", "transcript"],
        help="Data level (gene or transcript)")
    parser_wgcna.add_argument("-nt", "--network-type", dest="network_type", default="signed hybrid",
        choices=["unsigned", "signed", "signed hybrid"], help="Type of network")
    parser_wgcna.add_argument("-tom", "--tom-type", dest="tom_type", default="signed",
        choices=["unsigned", "signed"], help="Type of TOM")
    parser_wgcna.add_argument("-min", "--min-module-size", dest="min_module_size", type=int, default=50,
        help="Minimum module size")
    parser_wgcna.add_argument("-p", "--power", type=int,
        help="Soft thresholding power (auto-detected if not specified)")
    parser_wgcna.add_argument("-rsquared", "--rsquared-cut", dest="rsquared_cut", type=float, default=0.9,
        help="R-squared cutoff")
    parser_wgcna.add_argument("-mean", "--mean-cut", dest="mean_cut", type=int, default=100,
        help="Mean connectivity cutoff")
    parser_wgcna.add_argument("-mediss", "--mediss-thresh", dest="mediss_thresh", type=float, default=0.2,
        help="Module merging threshold")
    parser_wgcna.add_argument("-tpm", "--tpm-cutoff", dest="tpm_cutoff", type=int, default=1,
        help="TPM cutoff for filtering")
    
    parser_wgcna.set_defaults(func=main_wgcna)
    
    # all command (quant -> deseq2)
    parser_all = subparsers.add_parser("all", 
        help="Complete pipeline: quantification + DESeq2 analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
Coldata format (CSV):
    sample,id,condition,r1,r2
    sample1,ctrl,control,sample1_R1.fq.gz,sample1_R2.fq.gz
    sample2,treat,treatment,sample2_R1.fq.gz,sample2_R2.fq.gz

Examples:
    rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/
    rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ -p 24
        """)
    
    parser_all.add_argument("-S", "--coldata", required=True,
        help="Coldata file (CSV/TSV) with columns: sample,id,condition,r1,r2")
    parser_all.add_argument("-g", "--genome-fasta", dest="genome_fasta", required=True,
        help="Genome FASTA file")
    parser_all.add_argument("-gtf", "--gtf-file", dest="gtf_file", required=True,
        help="GTF annotation file")
    parser_all.add_argument("-gf", "--transcript-fasta", dest="transcript_fasta", required=True,
        help="Transcript FASTA file")
    parser_all.add_argument("-o", "--output-dir", dest="output_dir", required=True,
        help="Output directory (work directory)")
    parser_all.add_argument("-idx", "--index-dir", dest="index_dir",
        help="STAR index directory (default: <output_dir>/00_index)")
    parser_all.add_argument("-t2g", "--tx2gene", dest="tx2gene",
        help="Path to transcript-to-gene mapping file")
    parser_all.add_argument("-t", "--threads", type=int, default=8,
        help="Number of threads per sample")
    parser_all.add_argument("-p", "--parallel", type=int,
        help="Total cores for parallel processing")
    parser_all.add_argument("--trim", action="store_true",
        help="Trim reads with fastp")
    parser_all.add_argument("--force-index", action="store_true",
        help="Force rebuild index")
    parser_all.add_argument("--skip-existing", action="store_true",
        help="Skip samples if output already exists")
    parser_all.add_argument("--design", default="~condition",
        help="Design formula (e.g., '~condition', '~batch + condition')")
    parser_all.add_argument("--contrast",
        help="Contrast specification (e.g., 'condition,treatment,control')")
    parser_all.add_argument("--alpha", type=float, default=0.05,
        help="Significance threshold for adjusted p-values")
    parser_all.add_argument("--lfc", dest="lfc_threshold", type=float, default=2.0,
        help="Log2 fold change threshold for significant genes")
    
    parser_all.set_defaults(func=main_all)
    
    args = parser.parse_args()
    if vars(args) == {}:
        parser.print_help(sys.stderr)
    else:
        try:
            args.func(args)
        except Exception as e:
            logger.error(f"Error: {e}")
            raise


if __name__ == "__main__":
    main()
