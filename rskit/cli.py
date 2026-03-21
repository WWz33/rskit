import argparse
import sys
import os
import pandas as pd
import subprocess
from pathlib import Path
from rskit.config import StarConfig, SalmonConfig, PipelineConfig, DESeq2Config
from rskit.core.pipeline import RNAseqPipeline
from rskit.core.deseq2 import Deseq2Analyzer, run_deseq2_cli
from rskit.core.wgcna import run_wgcna_cli
from rskit.utils.logger import get_logger

logger = get_logger(__name__)

def setup_workdir(output_dir):
    """Setup work directory structure"""
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    
    dirs = {
        'index': output_dir / '00_index',
        'clean_data': output_dir / '01_clean_data',
        'clean_data_html': output_dir / '01_clean_data' / 'html',
        'clean_data_json': output_dir / '01_clean_data' / 'json',
        'bam': output_dir / '02_bam',
        'quant': output_dir / '03_quant',
        'deseq2': output_dir / '04_deseq2'
    }
    
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    
    return dirs

def trim_reads(read1, read2, sample, clean_data_dir, html_dir, json_dir, threads=8):
    """Run fastp trimming"""
    prefix = Path(clean_data_dir) / sample
    cmd = [
        'fastp',
        '-i', str(read1),
        '-I', str(read2),
        '-o', str(prefix) + '_1.fq',
        '-O', str(prefix) + '_2.fq',
        '--thread', str(threads),
        '-j', str(Path(json_dir) / sample) + '.json',
        '-h', str(Path(html_dir) / sample) + '.html'
    ]
    logger.info(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    return str(prefix) + '_1.fq', str(prefix) + '_2.fq'

def main_quant(args):
    """Run quantification pipeline"""
    from rskit.utils.validators import check_and_prepare_index
    
    # Convert paths to absolute before changing directory
    r1 = Path(args.r1).resolve() if args.r1 else None
    r2 = Path(args.r2).resolve() if args.r2 else None
    genome_fasta = Path(args.genome_fasta).resolve()
    gtf_file = Path(args.gtf_file).resolve()
    transcript_fasta = Path(args.transcript_fasta).resolve()
    
    # Setup work directory
    workdirs = setup_workdir(args.output_dir)
    
    # Determine and check index directory
    if args.index_dir:
        index_dir, needs_build = check_and_prepare_index(args.index_dir, args.force_index)
    else:
        index_dir, needs_build = check_and_prepare_index(workdirs['index'], args.force_index)
    
    # Parse samples
    if args.coldata:
        # Auto-detect separator based on file extension
        sep = '\t' if args.coldata.endswith('.tsv') else ','
        samples_df = pd.read_csv(args.coldata, sep=sep)
        required_cols = {'sample', 'r1', 'r2'}
        if not required_cols.issubset(samples_df.columns):
            raise ValueError(f"Sample file must contain columns: {required_cols}")
        samples_list = [(row['sample'], Path(row['r1']).resolve(), Path(row['r2']).resolve()) 
                        for _, row in samples_df.iterrows()]
    else:
        if not all([args.sample, r1, r2]):
            raise ValueError("Must provide --sample, --r1, --r2 or use --coldata")
        samples_list = [(args.sample, r1, r2)]
    
    # Process samples
    samples = {}
    for sample_name, r1_path, r2_path in samples_list:
        if args.trim:
            logger.info(f"Trimming {sample_name}...")
            r1_clean, r2_clean = trim_reads(r1_path, r2_path, sample_name, workdirs['clean_data'], 
                                           workdirs['clean_data_html'], workdirs['clean_data_json'], args.threads)
            samples[sample_name] = {'fq1': r1_clean, 'fq2': r2_clean}
        else:
            samples[sample_name] = {'fq1': str(r1_path), 'fq2': str(r2_path)}
    
    config = PipelineConfig(
        star=StarConfig(threads=args.threads),
        salmon=SalmonConfig(threads=args.threads),
        output_dir=str(workdirs['bam'])
    )
    pipeline = RNAseqPipeline(config)
    results = pipeline.run(
        samples=samples,
        genome_fasta=str(genome_fasta),
        gtf_file=str(gtf_file),
        transcript_fasta=str(transcript_fasta),
        index_dir=str(index_dir),
        output_dir=str(workdirs['bam']),
        quant_output_dir=str(workdirs['quant']),
        force_index=args.force_index
    )
    logger.info(f"Pipeline completed. Results: {results}")

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
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = workdirs['deseq2']
    
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
    
    # Validate input arguments
    if not args.expression:
        logger.error("Expression file must be provided")
        sys.exit(1)
    
    # Run WGCNA analysis
    try:
        analyzer = run_wgcna_cli(args)
        logger.info("WGCNA analysis completed successfully!")
    except Exception as e:
        logger.error(f"WGCNA analysis failed: {e}")
        raise

def main_all(args):
    """Run complete pipeline: quant -> deseq2"""
    from rskit.utils.validators import check_and_prepare_index
    
    logger.info("="*60)
    logger.info("Complete Pipeline: Quantification + DESeq2")
    logger.info("="*60)
    
    # Convert paths to absolute
    genome_fasta = Path(args.genome_fasta).resolve()
    gtf_file = Path(args.gtf_file).resolve()
    transcript_fasta = Path(args.transcript_fasta).resolve()
    coldata_file = Path(args.coldata).resolve()
    
    # Setup work directory
    workdirs = setup_workdir(args.output_dir)
    
    # Determine and check index directory
    if args.index_dir:
        index_dir, needs_build = check_and_prepare_index(args.index_dir, args.force_index)
    else:
        index_dir, needs_build = check_and_prepare_index(workdirs['index'], args.force_index)
    
    # Parse samples from coldata
    sep = '\t' if str(coldata_file).endswith('.tsv') else ','
    samples_df = pd.read_csv(coldata_file, sep=sep)
    required_cols = {'sample', 'r1', 'r2'}
    if not required_cols.issubset(samples_df.columns):
        raise ValueError(f"Coldata file must contain columns: {required_cols}")
    
    samples_list = [(row['sample'], Path(row['r1']).resolve(), Path(row['r2']).resolve()) 
                    for _, row in samples_df.iterrows()]
    
    # Step 1: Run quantification
    logger.info("="*60)
    logger.info("Step 1: Quantification Pipeline")
    logger.info("="*60)
    
    samples = {}
    for sample_name, r1_path, r2_path in samples_list:
        if args.trim:
            logger.info(f"Trimming {sample_name}...")
            r1_clean, r2_clean = trim_reads(r1_path, r2_path, sample_name, workdirs['clean_data'], 
                                           workdirs['clean_data_html'], workdirs['clean_data_json'], args.threads)
            samples[sample_name] = {'fq1': r1_clean, 'fq2': r2_clean}
        else:
            samples[sample_name] = {'fq1': str(r1_path), 'fq2': str(r2_path)}
    
    config = PipelineConfig(
        star=StarConfig(threads=args.threads),
        salmon=SalmonConfig(threads=args.threads),
        output_dir=str(workdirs['bam'])
    )
    pipeline = RNAseqPipeline(config)
    results = pipeline.run(
        samples=samples,
        genome_fasta=str(genome_fasta),
        gtf_file=str(gtf_file),
        transcript_fasta=str(transcript_fasta),
        index_dir=str(index_dir),
        output_dir=str(workdirs['bam']),
        quant_output_dir=str(workdirs['quant']),
        force_index=args.force_index
    )
    logger.info(f"Quantification completed. Results: {results}")
    
    # Step 2: Run DESeq2
    logger.info("="*60)
    logger.info("Step 2: DESeq2 Differential Expression Analysis")
    logger.info("="*60)
    
    # Create args namespace for deseq2
    class Deseq2Args:
        def __init__(self):
            self.salmon_dir = str(workdirs['quant'])
            self.gene_counts = None
            self.coldata = str(coldata_file)
            self.gtf = str(gtf_file)
            self.tx2gene = args.tx2gene
            self.work_dir = args.output_dir
            self.output_dir = str(workdirs['deseq2'])
            self.design = args.design
            self.contrast = args.contrast
            self.alpha = args.alpha
            self.threads = args.threads
    
    deseq2_args = Deseq2Args()
    
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
    parser_quant.add_argument("-S", "--coldata", help="Sample file (CSV/TSV) with columns: sample, r1, r2 (and optionally id, condition). Auto-detects separator based on file extension (.tsv for tab, .csv for comma)")
    parser_quant.add_argument("-1", "--r1", help="First read file")
    parser_quant.add_argument("-2", "--r2", help="Second read file")
    parser_quant.add_argument("-g", "--genome-fasta", dest="genome_fasta", required=True, help="Genome FASTA file")
    parser_quant.add_argument("-gtf", "--gtf-file", dest="gtf_file", required=True, help="GTF annotation file")
    parser_quant.add_argument("-gf", "--transcript-fasta", dest="transcript_fasta", required=True, help="Transcript FASTA file")
    parser_quant.add_argument("-o", "--output-dir", dest="output_dir", required=True, help="Output directory (work directory)")
    parser_quant.add_argument("-idx", "--index-dir", dest="index_dir", help="STAR index directory (default: <output_dir>/00_index)")
    parser_quant.add_argument("-t", "--threads", type=int, default=8, help="Number of threads")
    parser_quant.add_argument("--trim", action="store_true", help="Trim reads with fastp")
    parser_quant.add_argument("--force-index", action="store_true", help="Force rebuild index")
    parser_quant.set_defaults(func=main_quant)
    
    # deseq2 command
    parser_deseq2 = subparsers.add_parser("deseq2", 
        help="DESeq2 differential expression analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
Coldata format example (CSV):
    sample,id,condition
    lhy-D-rep1,lhy_D,lhy-D
    lhy-D-rep2,lhy_D,lhy-D
    lhy-rep1,lhy,lhy
    lhy-rep2,lhy,lhy
    WT-D-rep1,WT_D,WT-D
    WT-D-rep2,WT_D,WT-D
    WT-rep1,WT,WT
    WT-rep2,WT,WT

Examples:
    # Basic analysis with default design (~condition)
    rskit deseq2 --salmon-dir ./03_quant --coldata coldata.csv --gtf annotation.gtf
    
    # Multi-factor design with batch effect
    rskit deseq2 --salmon-dir ./03_quant --coldata coldata.csv --gtf annotation.gtf --design "~id + condition"
    
    # Specify contrast explicitly
    rskit deseq2 --salmon-dir ./03_quant --coldata coldata.csv --gtf annotation.gtf --contrast "condition,lhy-D,WT"
        """)
    
    # Input options (mutually exclusive)
    input_group = parser_deseq2.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-sd", "--salmon-dir", dest="salmon_dir", 
        help="Directory containing Salmon quant subfolders (uses pytximport for gene count estimation)")
    input_group.add_argument("-gc", "--gene-counts", dest="gene_counts",
        help="Path to gene counts matrix file (CSV/TSV, samples x genes)")
    
    # Required arguments
    parser_deseq2.add_argument("-S", "--coldata", required=True,
        help="Path to coldata/metadata file with columns: sample,id,condition (sample=folder name, id=group, condition=DE factor)")
    parser_deseq2.add_argument("-gtf", "--gtf", 
        help="GTF annotation file (required when using --salmon-dir without --tx2gene)")
    parser_deseq2.add_argument("-t2g", "--tx2gene", dest="tx2gene",
        help="Path to transcript-to-gene mapping file (CSV/TSV with transcript_id,gene_id columns)")
    
    # Output options
    parser_deseq2.add_argument("-w", "--work-dir", dest="work_dir", default=".",
        help="Work directory (will create 00_index, 01_clean_data, 02_bam, 03_quant, 04_deseq2)")
    parser_deseq2.add_argument("-o", "--output-dir", dest="output_dir", default=None,
        help="Custom output directory for DESeq2 results (default: work_dir/04_deseq2)")
    
    # Analysis options
    parser_deseq2.add_argument("--design", default="~condition",
        help="Design formula (e.g., '~condition' for single factor, '~id + condition' for multi-factor)")
    parser_deseq2.add_argument("--contrast", 
        help="Contrast specification: 'factor,level1,level2' (e.g., 'condition,lhy-D,WT'). If not specified, will auto-detect from 'condition' column.")
    parser_deseq2.add_argument("--alpha", type=float, default=0.05,
        help="Significance threshold for adjusted p-values")
    parser_deseq2.add_argument("-t", "--threads", type=int, default=None,
        help="Number of threads for parallel processing")
    
    parser_deseq2.set_defaults(func=main_deseq2)
    
    # wgcna command
    parser_wgcna = subparsers.add_parser("wgcna", 
        help="WGCNA co-expression network analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
Examples:
    # Basic analysis with expression data
    rskit wgcna --expression expression.csv --output-dir ./wgcna_results
    
    # Analysis with sample and gene metadata
    rskit wgcna --expression expression.csv --sample-info sample_info.csv --gene-info gene_info.csv --output-dir ./wgcna_results
    
    # Custom network parameters
    rskit wgcna --expression expression.csv --output-dir ./wgcna_results --network-type signed --min-module-size 30
        """)
    
    # Required arguments
    parser_wgcna.add_argument("-e", "--expression", required=True,
        help="Path to gene expression matrix file (CSV/TSV, samples x genes)")
    parser_wgcna.add_argument("-o", "--output-dir", dest="output_dir", required=True,
        help="Output directory for WGCNA results")
    
    # Optional metadata files
    parser_wgcna.add_argument("-S", "--coldata", dest="coldata",
        help="Path to sample metadata file (CSV/TSV)")
    parser_wgcna.add_argument("-G", "--gene-info", dest="gene_info",
        help="Path to gene metadata file (CSV/TSV)")
    
    # Data options
    parser_wgcna.add_argument("-sep", "--sep", default=",",
        help="Separator for input files (comma or tab)")
    parser_wgcna.add_argument("-n", "--name", default="WGCNA",
        help="Name for the WGCNA analysis")
    parser_wgcna.add_argument("-s", "--species", 
        help="Species for enrichment analysis (Human, Mouse, Yeast, Fly, Fish, Worm)")
    parser_wgcna.add_argument("-l", "--level", default="gene", choices=["gene", "transcript"],
        help="Data level (gene or transcript)")
    
    # WGCNA parameters
    parser_wgcna.add_argument("-nt", "--network-type", dest="network_type", default="signed hybrid",
        choices=["unsigned", "signed", "signed hybrid"],
        help="Type of network to generate")
    parser_wgcna.add_argument("-tom", "--tom-type", dest="tom_type", default="signed",
        choices=["unsigned", "signed"],
        help="Type of topological overlap matrix")
    parser_wgcna.add_argument("-min", "--min-module-size", dest="min_module_size", type=int, default=50,
        help="Minimum module size")
    parser_wgcna.add_argument("-p", "--power", type=int,
        help="Soft thresholding power (auto-detected if not specified)")
    parser_wgcna.add_argument("-rsquared", "--rsquared-cut", dest="rsquared_cut", type=float, default=0.9,
        help="R-squared cutoff for scale-free topology")
    parser_wgcna.add_argument("-mean", "--mean-cut", dest="mean_cut", type=int, default=100,
        help="Mean connectivity cutoff")
    parser_wgcna.add_argument("-mediss", "--mediss-thresh", dest="mediss_thresh", type=float, default=0.2,
        help="Module eigengene dissimilarity threshold for merging")
    parser_wgcna.add_argument("-tpm", "--tpm-cutoff", dest="tpm_cutoff", type=int, default=1,
        help="TPM cutoff for filtering low-expression genes")
    
    parser_wgcna.set_defaults(func=main_wgcna)
    
    # all command (quant -> deseq2)
    parser_all = subparsers.add_parser("all", 
        help="Complete pipeline: quantification + DESeq2 analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
Coldata format (CSV):
    sample,id,condition,r1,r2
    sample1,ctrl,control,sample1_R1.fq.gz,sample1_R2.fq.gz
    sample2,ctrl,control,sample2_R1.fq.gz,sample2_R2.fq.gz
    sample3,treat,treatment,sample3_R1.fq.gz,sample3_R2.fq.gz
    sample4,treat,treatment,sample4_R1.fq.gz,sample4_R2.fq.gz

Examples:
    # Basic analysis with default design (~condition)
    rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/
    
    # Multi-factor design with batch effect
    rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ --design "~id + condition"
    
    # Specify contrast explicitly
    rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ --contrast "condition,treatment,control"
        """)
    
    # Required arguments
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
    
    # Optional arguments
    parser_all.add_argument("-idx", "--index-dir", dest="index_dir",
        help="STAR index directory (default: <output_dir>/00_index)")
    parser_all.add_argument("-t2g", "--tx2gene", dest="tx2gene",
        help="Path to transcript-to-gene mapping file (CSV/TSV with transcript_id,gene_id columns)")
    parser_all.add_argument("-t", "--threads", type=int, default=8,
        help="Number of threads")
    parser_all.add_argument("--trim", action="store_true",
        help="Trim reads with fastp")
    parser_all.add_argument("--force-index", action="store_true",
        help="Force rebuild index")
    
    # DESeq2 options
    parser_all.add_argument("--design", default="~condition",
        help="Design formula (e.g., '~condition' for single factor, '~id + condition' for multi-factor)")
    parser_all.add_argument("--contrast",
        help="Contrast specification: 'factor,level1,level2' (e.g., 'condition,treatment,control'). If not specified, will auto-detect from 'condition' column.")
    parser_all.add_argument("--alpha", type=float, default=0.05,
        help="Significance threshold for adjusted p-values")
    
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
