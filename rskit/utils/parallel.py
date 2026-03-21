from pathlib import Path
import pandas as pd
from rskit.utils.logger import get_logger

logger = get_logger(__name__)

def count_samples_from_coldata(coldata_file):
    """Count number of samples from coldata file.
    
    Args:
        coldata_file: Path to coldata file (CSV/TSV)
        
    Returns:
        int: Number of samples (rows excluding header)
    """
    coldata_path = Path(coldata_file)
    
    # Auto-detect separator
    sep = '\t' if coldata_path.suffix == '.tsv' else ','
    
    df = pd.read_csv(coldata_path, sep=sep)
    return len(df)

def calculate_threads_per_sample(parallel_cores, num_samples):
    """Calculate threads per sample based on total parallel cores and number of samples.
    
    Args:
        parallel_cores: Total number of cores for parallel processing
        num_samples: Number of samples to process
        
    Returns:
        int: Threads per sample (minimum 1)
    """
    if num_samples <= 0:
        return parallel_cores
    
    threads_per_sample = max(1, parallel_cores // num_samples)
    return threads_per_sample

def get_optimal_threads(parallel_cores, coldata_file=None, num_samples=None):
    """Get optimal threads per sample.
    
    Args:
        parallel_cores: Total number of cores specified by --parallel
        coldata_file: Path to coldata file (optional)
        num_samples: Number of samples (optional, will be calculated from coldata if not provided)
        
    Returns:
        tuple: (threads_per_sample, num_samples)
    """
    if num_samples is None:
        if coldata_file:
            num_samples = count_samples_from_coldata(coldata_file)
        else:
            num_samples = 1
    
    threads_per_sample = calculate_threads_per_sample(parallel_cores, num_samples)
    
    logger.info(f"Parallel processing: {parallel_cores} cores / {num_samples} samples = {threads_per_sample} threads per sample")
    
    return threads_per_sample, num_samples
