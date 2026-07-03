from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import List
from rskit.utils.logger import get_logger

logger = get_logger(__name__)


@dataclass(frozen=True)
class SamplePlan:
    """Resolved per-sample scheduling plan."""

    total_threads: int
    requested_jobs: int
    active_jobs: int
    threads_per_sample: int
    cap_reasons: List[str]


def calculate_sample_plan(total_threads: int, requested_jobs: int, num_samples: int) -> SamplePlan:
    """Resolve total thread budget and requested jobs into a bounded sample plan."""
    if total_threads <= 0:
        raise ValueError("--threads must be a positive integer")
    if requested_jobs <= 0:
        raise ValueError("--jobs must be a positive integer")
    if num_samples <= 0:
        raise ValueError("At least one sample is required")

    active_jobs = min(requested_jobs, num_samples, total_threads)
    cap_reasons = []
    if active_jobs < requested_jobs:
        if active_jobs == num_samples or num_samples <= total_threads:
            cap_reasons.append("sample count")
        if active_jobs == total_threads or total_threads < num_samples:
            cap_reasons.append("total threads")

    return SamplePlan(
        total_threads=total_threads,
        requested_jobs=requested_jobs,
        active_jobs=active_jobs,
        threads_per_sample=max(1, total_threads // active_jobs),
        cap_reasons=cap_reasons,
    )


def process_single_sample(args):
    """Process a single sample: alignment + quantification."""
    sample_name, sample_data, index_dir, transcript_fasta, workdirs, threads, skip_existing, star_args, salmon_args = args
    
    from rskit.core.star import StarAligner
    from rskit.core.salmon import SalmonQuantifier
    from rskit.config import StarConfig, SalmonConfig
    
    aligner = StarAligner(StarConfig(threads=threads, extra_args=star_args))
    quantifier = SalmonQuantifier(SalmonConfig(threads=threads, extra_args=salmon_args))
    
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


def run_samples_parallel(
    samples,
    index_dir,
    transcript_fasta,
    workdirs,
    threads_per_sample,
    jobs,
    skip_existing=False,
    star_args="",
    salmon_args="",
):
    """Run alignment and quantification for multiple samples in parallel."""
    num_samples = len(samples)
    max_workers = min(jobs, num_samples)
    
    logger.info(
        f"Parallel processing: {num_samples} samples, {max_workers} active jobs, "
        f"{threads_per_sample} threads/sample"
    )
    
    sample_args = [
        (name, data, index_dir, transcript_fasta, workdirs, threads_per_sample, skip_existing, star_args, salmon_args)
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
