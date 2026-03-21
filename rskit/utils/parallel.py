from pathlib import Path
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from rskit.utils.logger import get_logger

logger = get_logger(__name__)

def count_samples_from_coldata(coldata_file):
    """Count number of samples from coldata file."""
    coldata_path = Path(coldata_file)
    sep = '\t' if coldata_path.suffix == '.tsv' else ','
    df = pd.read_csv(coldata_path, sep=sep)
    return len(df)

def calculate_threads_per_sample(parallel_cores, num_samples):
    """Calculate threads per sample based on total parallel cores and number of samples."""
    if num_samples <= 0:
        return parallel_cores
    return max(1, parallel_cores // num_samples)


def process_single_sample(args):
    """Process a single sample: alignment + quantification."""
    sample_name, sample_data, index_dir, transcript_fasta, workdirs, threads, skip_existing = args
    
    from rskit.core.star import StarAligner
    from rskit.core.salmon import SalmonQuantifier
    from rskit.config import StarConfig, SalmonConfig
    
    aligner = StarAligner(StarConfig(threads=threads))
    quantifier = SalmonQuantifier(SalmonConfig(threads=threads))
    
    # Setup output directories
    sample_bam_dir = Path(workdirs['bam']) / sample_name
    sample_bam_dir.mkdir(parents=True, exist_ok=True)
    
    sample_quant_dir = Path(workdirs['quant']) / sample_name
    sample_quant_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if output exists and skip if requested
    quant_file = sample_quant_dir / "quant.sf"
    if skip_existing and quant_file.exists():
        logger.info(f"[{sample_name}] Output exists, skipping")
        return sample_name, {"quantification": {"quant": str(quant_file)}}
    
    # Run alignment
    logger.info(f"[{sample_name}] Aligning with {threads} threads...")
    align_prefix = str(sample_bam_dir / f"{sample_name}_")
    align_results = aligner.align(index_dir, sample_data["fq1"], sample_data["fq2"],
                                  align_prefix, sample_name=sample_name)
    
    # Run quantification
    logger.info(f"[{sample_name}] Quantifying...")
    quant_results = quantifier.quantify(transcript_fasta, align_results["transcriptome_bam"],
                                        str(sample_quant_dir), sample_name=sample_name,
                                        skip_if_exists=skip_existing)
    
    logger.info(f"[{sample_name}] Completed")
    return sample_name, {"alignment": align_results, "quantification": quant_results}


def run_samples_parallel(samples, index_dir, transcript_fasta, workdirs, threads_per_sample, skip_existing=False):
    """Run alignment and quantification for multiple samples in parallel."""
    num_samples = len(samples)
    max_workers = num_samples
    
    logger.info(f"Parallel processing: {num_samples} samples, {threads_per_sample} threads each")
    
    sample_args = [
        (name, data, index_dir, transcript_fasta, workdirs, threads_per_sample, skip_existing)
        for name, data in samples.items()
    ]
    
    results = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_single_sample, args): args[0] for args in sample_args}
        
        for i, future in enumerate(as_completed(futures), 1):
            sample_name, sample_results = future.result()
            results[sample_name] = sample_results
            logger.info(f"Progress: {i}/{num_samples} completed")
    
    logger.info(f"All {num_samples} samples completed!")
    return results
