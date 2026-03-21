from pathlib import Path
from rskit.utils.logger import get_logger

logger = get_logger(__name__)

def validate_file(file_path: str) -> bool:
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    if not path.is_file():
        raise ValueError(f"Not a file: {file_path}")
    return True

def validate_dir(dir_path: str) -> bool:
    path = Path(dir_path)
    if not path.exists():
        raise FileNotFoundError(f"Directory not found: {dir_path}")
    if not path.is_dir():
        raise ValueError(f"Not a directory: {dir_path}")
    return True

def check_star_index(index_dir: str) -> bool:
    index_path = Path(index_dir)
    required_files = ["SA", "SAindex", "Genome", "genomeParameters.txt"]
    return all((index_path / f).exists() for f in required_files)

def check_salmon_index(index_dir: str) -> bool:
    index_path = Path(index_dir)
    required_files = ["sa.bin", "txpInfo.bin", "refInfo.json", "versionInfo.json"]
    return all((index_path / f).exists() for f in required_files)

def check_and_prepare_index(index_dir, force_index=False):
    """Check if STAR index exists and is complete.
    
    Args:
        index_dir: Path to index directory
        force_index: If True, will rebuild even if index exists
        
    Returns:
        tuple: (index_path, needs_build)
            - index_path: Path object of index directory
            - needs_build: True if index needs to be built
    """
    index_path = Path(index_dir).resolve()
    
    if index_path.exists() and check_star_index(str(index_path)):
        if force_index:
            logger.info(f"Valid STAR index found at {index_path}, but --force-index specified, will rebuild")
            return index_path, True
        else:
            logger.info(f"Valid STAR index found at {index_path}, skipping index building")
            return index_path, False
    else:
        if force_index:
            logger.info(f"Will rebuild index at {index_path}")
        else:
            logger.info(f"No valid index found at {index_path}, will build new index")
        return index_path, True
