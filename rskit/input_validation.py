from typing import List, Optional, Sequence, Tuple

import pandas as pd

from rskit.input_contracts import (
    design_columns,
    load_coldata,
    read_table,
    resolve_path_from_table,
    validate_sample_alignment,
)
from rskit.utils.logger import get_logger


logger = get_logger(__name__)


def validate_input_files(
    coldata: str,
    design: str = "~condition",
    check_reads: bool = False,
    gene_counts: Optional[str] = None,
    expression: Optional[str] = None,
) -> List[str]:
    """Validate user-provided input files without running analysis tools."""
    required_columns = list(design_columns(design))
    if check_reads:
        required_columns.extend(["r1", "r2"])

    metadata = load_coldata(coldata, required_columns=_unique(required_columns))
    messages = [f"coldata: {len(metadata.index)} samples"]

    if check_reads:
        _validate_read_paths(coldata, metadata)
        messages.append("reads: r1/r2 paths exist")

    if gene_counts:
        counts = read_table(gene_counts, index_col=0)
        samples, genes = _validate_gene_by_sample_table(counts, metadata, table_name="gene counts")
        messages.append(
            "gene counts: "
            + f"{samples} samples x {genes} genes"
        )

    if expression:
        expression_table = read_table(expression, index_col=0)
        samples, genes = _validate_gene_by_sample_table(expression_table, metadata, table_name="expression matrix")
        messages.append(
            "expression: "
            + f"{samples} samples x {genes} genes"
        )

    return messages


def _validate_read_paths(coldata: str, metadata: pd.DataFrame) -> None:
    missing_paths = []
    for sample_name, row in metadata.iterrows():
        for column in ("r1", "r2"):
            read_path = resolve_path_from_table(row[column], coldata)
            if not read_path.exists():
                missing_paths.append(f"{sample_name}:{column}={read_path}")

    if missing_paths:
        raise ValueError("Read files do not exist: " + _preview(missing_paths))


def _validate_gene_by_sample_table(
    table: pd.DataFrame,
    metadata: pd.DataFrame,
    table_name: str,
) -> Tuple[int, int]:
    metadata_samples = {str(sample_id) for sample_id in metadata.index}
    table_rows = {str(sample_id) for sample_id in table.index}
    table_columns = {str(column) for column in table.columns}

    appears_samples_by_genes = metadata_samples.issubset(table_rows) and not metadata_samples.issubset(table_columns)
    if appears_samples_by_genes:
        logger.warning(
            f"{table_name} appears to be samples x genes, but rskit expects "
            "genes x samples for user input; please transpose the file."
        )
        validate_sample_alignment(table, metadata, table_name=table_name)
        return table.shape

    oriented = table.T
    validate_sample_alignment(oriented, metadata, table_name=table_name)

    return oriented.shape


def _unique(values: Sequence[str]) -> List[str]:
    result = []
    for value in values:
        if value not in result:
            result.append(value)
    return result


def _preview(values: Sequence[str], limit: int = 5) -> str:
    preview = ", ".join(values[:limit])
    return preview + ("..." if len(values) > limit else "")
