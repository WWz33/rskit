import argparse
import sys
import os
import pandas as pd
import subprocess
from pathlib import Path
from rskit.config import StarConfig, SalmonConfig, PipelineConfig
from rskit.core.pipeline import RNAseqPipeline
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
        'quant': output_dir / '03_quant'
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
    # Convert paths to absolute before changing directory
    r1 = Path(args.r1).resolve()
    r2 = Path(args.r2).resolve()
    genome_fasta = Path(args.genome_fasta).resolve()
    gtf_file = Path(args.gtf_file).resolve()
    transcript_fasta = Path(args.transcript_fasta).resolve()
    index_dir = Path(args.index_dir).resolve()
    
    # Setup work directory
    workdirs = setup_workdir(args.output_dir)
    
    # Parse samples
    if args.sample_tsv:
        samples_df = pd.read_csv(args.sample_tsv, sep='\t')
        required_cols = {'sample', 'r1', 'r2'}
        if not required_cols.issubset(samples_df.columns):
            raise ValueError(f"TSV must contain columns: {required_cols}")
        samples_list = [(row['sample'], Path(row['r1']).resolve(), Path(row['r2']).resolve()) 
                        for _, row in samples_df.iterrows()]
    else:
        if not all([args.sample, r1, r2]):
            raise ValueError("Must provide --sample, --r1, --r2 or use --sample-tsv")
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
    logger.info(f"DESeq2 analysis for: {args.quant_dir}")

def main():
    parser = argparse.ArgumentParser(
        description="RNA-seq analysis toolkit",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparsers = parser.add_subparsers(help="Task to perform")
    
    # quant command
    parser_quant = subparsers.add_parser("quant", help="Complete quantification pipeline (index -> align -> quant)")
    
    parser_quant.add_argument("-s", "--sample", help="Sample name")
    parser_quant.add_argument("-S", "--sample-tsv", help="TSV file with columns: sample, r1, r2")
    parser_quant.add_argument("-1", "--r1", required=True, help="First read file")
    parser_quant.add_argument("-2", "--r2", required=True, help="Second read file")
    parser_quant.add_argument("-g", "--genome-fasta", dest="genome_fasta", required=True, help="Genome FASTA file")
    parser_quant.add_argument("-gtf", "--gtf-file", dest="gtf_file", required=True, help="GTF annotation file")
    parser_quant.add_argument("-gf", "--transcript-fasta", dest="transcript_fasta", required=True, help="Transcript FASTA file")
    parser_quant.add_argument("-o", "--output-dir", dest="output_dir", required=True, help="Output directory (work directory)")
    parser_quant.add_argument("--index-dir", dest="index_dir", default="STAR_index", help="STAR index directory")
    parser_quant.add_argument("-t", "--threads", type=int, default=8, help="Number of threads")
    parser_quant.add_argument("--trim", action="store_true", help="Trim reads with fastp")
    parser_quant.add_argument("--force-index", action="store_true", help="Force rebuild index")
    parser_quant.set_defaults(func=main_quant)
    
    # deseq2 command
    parser_deseq2 = subparsers.add_parser("deseq2", help="DESeq2 differential expression analysis")
    parser_deseq2.add_argument("quant_dir", help="Quantification results directory")
    parser_deseq2.set_defaults(func=main_deseq2)
    
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
