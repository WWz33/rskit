import argparse
import sys
import pandas as pd
from rnaseq_tools.config import StarConfig, SalmonConfig, PipelineConfig
from rnaseq_tools.core.pipeline import RNAseqPipeline
from rnaseq_tools.utils.logger import get_logger

logger = get_logger(__name__)

def main_quant(args):
    # Parse samples
    if args.sample_tsv:
        samples_df = pd.read_csv(args.sample_tsv, sep='\t')
        required_cols = {'sample', 'r1', 'r2'}
        if not required_cols.issubset(samples_df.columns):
            raise ValueError(f"TSV must contain columns: {required_cols}")
        samples = {row['sample']: {'fq1': row['r1'], 'fq2': row['r2']} 
                   for _, row in samples_df.iterrows()}
    else:
        if not all([args.sample, args.r1, args.r2]):
            raise ValueError("Must provide --sample, --r1, --r2 or use --sample-tsv")
        samples = {args.sample: {'fq1': args.r1, 'fq2': args.r2}}
    
    config = PipelineConfig(
        star=StarConfig(threads=args.threads),
        salmon=SalmonConfig(threads=args.threads),
        output_dir=args.output_dir
    )
    pipeline = RNAseqPipeline(config)
    results = pipeline.run(
        samples=samples,
        genome_fasta=args.genome_fasta,
        gtf_file=args.gtf_file,
        transcript_fasta=args.transcript_fasta,
        index_dir=args.index_dir,
        output_dir=args.output_dir,
        force_index=args.force_index
    )
    logger.info(f"Pipeline completed. Results: {results}")

def main_deseq2(args):
    logger.info(f"DESeq2 analysis for: {args.quant_dir}")

def main():
    parser = argparse.ArgumentParser(
        description="RNA-seq analysis tools",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparsers = parser.add_subparsers(help="Task to perform")
    
    # quant command
    parser_quant = subparsers.add_parser("quant", help="Complete quantification pipeline (index -> align -> quant)")
    parser_quant.add_argument("-g", "--genome_fasta", help="Genome FASTA file")
    parser_quant.add_argument("-gtf, ""gtf_file", help="GTF annotation file")
    parser_quant.add_argument("-tf", "transcript_fasta", help="Transcript FASTA file")
    parser_quant.add_argument("output_dir", help="Output directory")
    
    sample_group = parser_quant.add_mutually_exclusive_group(required=True)
    sample_group.add_argument("-s", "--sample", help="Sample name")
    sample_group.add_argument("-S", "--sample-tsv", help="TSV file with columns: sample, r1, r2")
    
    parser_quant.add_argument("-1", "--r1", help="First read file (required with -s)")
    parser_quant.add_argument("-2", "--r2", help="Second read file (required with -s)")
    parser_quant.add_argument("--index-dir", default="STAR_index", help="STAR index directory")
    parser_quant.add_argument("-t", "--threads", type=int, default=56, help="Number of threads")
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
