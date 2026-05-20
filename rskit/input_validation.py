from typing import List, Optional, Sequence

import pandas as pd

from rskit.input_contracts import (
    design_columns,
    load_coldata,
    orient_sample_table,
    read_table,
    resolve_path_from_table,
    validate_sample_alignment,
)


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
        oriented_counts = orient_sample_table(counts, metadata, table_name="gene counts")
        messages.append(
            "gene counts: "
            + f"{oriented_counts.shape[0]} samples x {oriented_counts.shape[1]} genes"
        )

    if expression:
        expression_table = read_table(expression, index_col=0)
        _validate_expression_table(expression_table, metadata)
        messages.append(
            "expression: "
            + f"{expression_table.shape[0]} samples x {expression_table.shape[1]} genes"
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


def _validate_expression_table(expression: pd.DataFrame, metadata: pd.DataFrame) -> None:
    metadata_samples = {str(sample_id) for sample_id in metadata.index}
    expression_rows = {str(sample_id) for sample_id in expression.index}
    expression_columns = {str(column) for column in expression.columns}

    if metadata_samples.issubset(expression_columns) and not metadata_samples.issubset(expression_rows):
        raise ValueError(
            "Expression matrix appears transposed: sample IDs match columns, "
            "but rows must be samples and columns must be genes"
        )

    validate_sample_alignment(expression, metadata, table_name="expression matrix")


def _unique(values: Sequence[str]) -> List[str]:
    result = []
    for value in values:
        if value not in result:
            result.append(value)
    return result


def _preview(values: Sequence[str], limit: int = 5) -> str:
    preview = ", ".join(values[:limit])
    return preview + ("..." if len(values) > limit else "")
