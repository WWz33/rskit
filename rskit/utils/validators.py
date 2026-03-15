from pathlib import Path

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
