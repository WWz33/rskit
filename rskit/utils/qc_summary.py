"""Aggregate per-sample QC metrics into one summary table.

Collects fastp json, STAR Log.final.out and salmon meta_info.json for each
sample and writes 00_summary/summary.csv plus a human-readable copy, so a
finished run has a single place to check trimming, alignment and mapping
rates without walking per-sample directories.
"""

import json
from pathlib import Path
from typing import Dict, Optional

from rskit.utils.logger import get_logger

logger = get_logger(__name__)

SUMMARY_COLUMNS = [
    "sample",
    "input_reads",
    "clean_reads",
    "clean_bases",
    "q30_rate",
    "star_total_reads",
    "star_unique_mapped",
    "star_unique_rate",
    "star_multimapped",
    "salmon_mapped_reads",
    "salmon_mapping_rate",
]


def _fastp_metrics(json_path: Path) -> Dict[str, object]:
    report = json.loads(json_path.read_text(encoding="utf-8"))
    summary = report.get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})
    clean_reads = after.get("total_reads", 0)
    return {
        "input_reads": before.get("total_reads", 0),
        "clean_reads": clean_reads,
        "clean_bases": after.get("total_bases", 0),
        "q30_rate": after.get("q30_rate", ""),
    }


def _star_metrics(log_path: Path) -> Dict[str, object]:
    metrics = {}
    for line in log_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if "|" not in line:
            continue
        key, _, value = line.partition("|")
        key = key.strip()
        if key in (
            "Number of input reads",
            "Uniquely mapped reads number",
            "Uniquely mapped reads %",
            "Number of reads mapped to multiple loci",
        ):
            metrics[key] = value.strip()
    unique_rate = metrics.get("Uniquely mapped reads %", "")
    if unique_rate:
        unique = float(unique_rate.rstrip("%") or 0)
        unique_rate = f"{unique / 100:.4f}"
    return {
        "star_total_reads": metrics.get("Number of input reads", ""),
        "star_unique_mapped": metrics.get("Uniquely mapped reads number", ""),
        "star_unique_rate": unique_rate,
        "star_multimapped": metrics.get("Number of reads mapped to multiple loci", ""),
    }


def _salmon_metrics(meta_path: Path) -> Dict[str, object]:
    info = json.loads(meta_path.read_text(encoding="utf-8"))
    observed = info.get("num_processed", 0)
    mapped = info.get("num_mapped", 0)
    return {
        "salmon_mapped_reads": mapped,
        "salmon_mapping_rate": f"{mapped / observed:.4f}" if observed else "",
    }


def write_qc_summary(workdirs: Dict[str, Path], skip_trimming: bool = False) -> Optional[Path]:
    """Write 00_summary/summary.csv from per-sample QC files. Returns its path."""
    rows = []
    quant_dir = workdirs["quant"]
    for sample_dir in sorted(p for p in quant_dir.iterdir() if p.is_dir()):
        row: Dict[str, object] = {"sample": sample_dir.name}
        fastp_json = workdirs["clean_data_json"] / f"{sample_dir.name}.json"
        if fastp_json.exists():
            row.update(_fastp_metrics(fastp_json))
        star_log = workdirs["bam"] / sample_dir.name / f"{sample_dir.name}_Log.final.out"
        if star_log.exists():
            row.update(_star_metrics(star_log))
        salmon_meta = sample_dir / "aux_info" / "meta_info.json"
        if salmon_meta.exists():
            row.update(_salmon_metrics(salmon_meta))
        rows.append(row)

    if not rows:
        return None

    summary_dir = workdirs.get("summary")
    if summary_dir is None:
        return None
    summary_dir.mkdir(parents=True, exist_ok=True)

    import pandas as pd

    table = pd.DataFrame(rows, columns=SUMMARY_COLUMNS)
    summary_path = summary_dir / "summary.csv"
    table.to_csv(summary_path, index=False)
    logger.info(f"QC summary written to {summary_path} ({len(rows)} samples)")
    return summary_path
