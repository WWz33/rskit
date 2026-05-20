from pathlib import Path
from typing import Optional, Sequence

import pandas as pd


def detect_separator(path: str, sep: Optional[str] = None) -> str:
    """Return the table separator for a path, unless explicitly provided."""
    if sep is not None:
        return sep
    return "\t" if Path(path).suffix.lower() in {".tsv", ".txt"} else ","


def read_table(path: str, sep: Optional[str] = None, index_col=None) -> pd.DataFrame:
    """Read a CSV/TSV-style table using the repository's separator contract."""
    return pd.read_csv(path, sep=detect_separator(path, sep), index_col=index_col)


def load_coldata(path: str, required_columns: Sequence[str] = ()) -> pd.DataFrame:
    """Load sample metadata with a required sample identifier column."""
    metadata = read_table(path)
    if "sample" not in metadata.columns:
        raise ValueError("Coldata file must contain a 'sample' column")

    missing_columns = [column for column in required_columns if column not in metadata.columns]
    if missing_columns:
        raise ValueError(
            "Coldata file is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    return metadata.set_index("sample")


def validate_sample_alignment(
    table: pd.DataFrame,
    metadata: pd.DataFrame,
    table_name: str = "input table",
) -> None:
    """Require table rows to match metadata sample IDs exactly."""
    table_samples = {str(sample_id) for sample_id in table.index}
    metadata_samples = {str(sample_id) for sample_id in metadata.index}

    missing_from_table = sorted(metadata_samples - table_samples)
    extra_in_table = sorted(table_samples - metadata_samples)

    problems = []
    if missing_from_table:
        problems.append(
            "missing from "
            + table_name
            + ": "
            + _preview_values(missing_from_table)
        )
    if extra_in_table:
        problems.append(
            "not present in coldata: "
            + _preview_values(extra_in_table)
        )

    if problems:
        raise ValueError("Sample IDs do not match (" + "; ".join(problems) + ")")


def _preview_values(values: Sequence[str], limit: int = 5) -> str:
    preview = ", ".join(values[:limit])
    return preview + ("..." if len(values) > limit else "")
