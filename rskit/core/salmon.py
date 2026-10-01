import gzip
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from rskit.cli_args import merge_extra_args
from rskit.config import SalmonConfig
from rskit.core.base import Tool, is_dry_run
from rskit.utils.gtf import iter_gtf
from rskit.utils.logger import get_logger
from rskit.utils.validators import validate_file

logger = get_logger(__name__)


def _gff3_value(value: Optional[str], prefix: str) -> Optional[str]:
    """Extract a GFF3 ID/Parent value, stripping the 'transcript:'/'gene:' prefix."""
    if not value:
        return None
    value = value.split(",")[0].strip()
    tag = prefix + ":"
    if value.startswith(tag):
        value = value[len(tag):]
    return value or None


# GTF annotation writes "transcript" records; GFF3 (Ensembl/NCBI) calls the
# same feature "mRNA", so both must be accepted when building tx2gene
TRANSCRIPT_FEATURES = frozenset({"transcript", "mRNA"})


def _columns_look_like_identifiers(columns) -> bool:
    """Header rows are words ('transcript_id'); transcript/gene IDs contain digits."""
    values = [str(column) for column in columns]
    return bool(values) and all(any(char.isdigit() for char in value) for value in values)

SALMON_QUANT_PROTECTED_OPTIONS = {
    "-t",
    "--targets",
    "-a",
    "--alignments",
    "-o",
    "--output",
    "-p",
    "--threads",
    "-l",
    "--libType",
    # read inputs (mapping-based / direct mode)
    "-r",
    "--unmatedReads",
    "-1",
    "--mates1",
    "-2",
    "--mates2",
}

class SalmonQuantifier:
    def __init__(self, config: SalmonConfig):
        self.config = config
        self.tool = Tool("salmon")
        self.logger = self.tool.logger
    
    def quantify(self, transcript_fasta: str, bam_file: str, output_dir: str,
                 sample_name: Optional[str] = None, skip_if_exists: bool = True) -> dict:
        if not is_dry_run():
            validate_file(transcript_fasta)
            validate_file(bam_file)
        
        output_path = Path(output_dir)
        quant_file = output_path / "quant.sf"
        
        # only trust a non-empty quant.sf (a crashed run can leave a truncated file)
        if skip_if_exists and quant_file.exists() and quant_file.stat().st_size > 0:
            self.logger.info(f"Quantification output already exists at {output_dir}, skipping")
            existing = {"quant": str(quant_file)}
            lib_counts = output_path / "lib_format_counts.json"
            if lib_counts.exists():
                existing["lib_format_counts"] = str(lib_counts)
            return existing
        
        output_path.mkdir(parents=True, exist_ok=True)
        
        cmd = ["salmon", "quant", "-t", transcript_fasta, "-l", self.config.lib_type,
               "-a", bam_file, "-o", output_dir, "-p", str(self.config.threads)]
        
        if self.config.seq_bias:
            cmd.append("--seqBias")
        if self.config.gc_bias:
            cmd.append("--gcBias")
        if self.config.pos_bias:
            cmd.append("--posBias")
        cmd = merge_extra_args(cmd, self.config.extra_args, SALMON_QUANT_PROTECTED_OPTIONS)
        
        self.logger.info(f"Quantifying {sample_name or 'sample'} with Salmon")
        self.tool._run_command(cmd)
        
        return {"quant": str(quant_file), "lib_format_counts": str(output_path / "lib_format_counts.json")}
    
    def quantify_from_reads(self, transcript_fasta: str, read1: str, read2: str,
                            output_dir: str, sample_name: Optional[str] = None,
                            skip_if_exists: bool = True) -> dict:
        """Quantify directly from reads (selective alignment), without STAR.

        Only the transcriptome is needed, so users without a genome build can
        still get gene-level counts/TPM.
        """
        if not is_dry_run():
            validate_file(transcript_fasta)
            validate_file(read1)
            validate_file(read2)

        output_path = Path(output_dir)
        quant_file = output_path / "quant.sf"
        if skip_if_exists and quant_file.exists() and quant_file.stat().st_size > 0:
            self.logger.info(f"Quantification output already exists at {output_dir}, skipping")
            existing = {"quant": str(quant_file)}
            lib_counts = output_path / "lib_format_counts.json"
            if lib_counts.exists():
                existing["lib_format_counts"] = str(lib_counts)
            return existing

        output_path.mkdir(parents=True, exist_ok=True)

        cmd = ["salmon", "quant", "-t", transcript_fasta, "-l", self.config.lib_type,
               "-1", read1, "-2", read2, "-o", output_dir, "-p", str(self.config.threads)]
        if self.config.seq_bias:
            cmd.append("--seqBias")
        if self.config.gc_bias:
            cmd.append("--gcBias")
        if self.config.pos_bias:
            cmd.append("--posBias")
        cmd = merge_extra_args(cmd, self.config.extra_args, SALMON_QUANT_PROTECTED_OPTIONS)

        self.logger.info(f"Quantifying {sample_name or 'sample'} from reads with Salmon")
        self.tool._run_command(cmd)

        return {"quant": str(quant_file), "lib_format_counts": str(output_path / "lib_format_counts.json")}

    def validate_inputs(self) -> bool:
        return self.tool._check_tool_installed()


class SalmonExpressionExporter:
    def __init__(self):
        self.logger = logger

    def _create_tx2gene_from_gtf(self, gtf_file: str, output_dir: Optional[str] = None) -> pd.DataFrame:
        self.logger.info(f"Parsing GTF/GFF3 file: {gtf_file}")

        tx2gene_map = {}
        num_records = 0

        opener = gzip.open if str(gtf_file).endswith(".gz") else open
        with opener(gtf_file, "rt", encoding="utf-8", errors="ignore") as reader:
            for rec in iter_gtf(reader, "ensembl"):
                num_records += 1
                if rec.feature not in TRANSCRIPT_FEATURES:
                    continue
                meta = rec.meta or {}
                transcript_id = meta.get("transcript_id") or _gff3_value(meta.get("ID"), "transcript")
                gene_id = meta.get("gene_id") or _gff3_value(meta.get("Parent"), "gene")
                if transcript_id and gene_id:
                    tx2gene_map.setdefault(transcript_id, gene_id)

        if not tx2gene_map:
            raise ValueError(
                f"No transcript-to-gene mappings extracted from {gtf_file}; "
                "check that it is a GTF with transcript records or a GFF3 with mRNA records"
            )

        self.logger.info(f"Scanned {num_records} GTF/GFF3 lines")
        self.logger.info(f"Extracted {len(tx2gene_map)} unique transcript-to-gene mappings")

        tx2gene_df = pd.DataFrame(
            [(tx, gene) for tx, gene in tx2gene_map.items()],
            columns=["transcript_id", "gene_id"],
        )

        if output_dir:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            tx2gene_file = output_path / "tx2gene.tsv"
            tx2gene_df.to_csv(tx2gene_file, sep="\t", index=False)
            self.logger.info(f"Saved tx2gene mapping to: {tx2gene_file}")

        return tx2gene_df

    def _load_tx2gene_map(
        self,
        gtf_file: Optional[str] = None,
        tx2gene: Optional[str] = None,
        output_dir: Optional[str] = None,
    ) -> Tuple[pd.DataFrame, Optional[Path]]:
        tx2gene_file = None

        if tx2gene is not None:
            tx2gene_path = Path(tx2gene)
            separator = "\t" if tx2gene_path.suffix.lower() in {".tsv", ".txt"} else ","
            tx2gene_map = pd.read_csv(tx2gene_path, sep=separator)
            self.logger.info(f"Loaded tx2gene map from {tx2gene_path}")

            if "transcript_id" not in tx2gene_map.columns or "gene_id" not in tx2gene_map.columns:
                # headerless files are common; pandas treats their first row as column names.
                # Ensembl-style IDs (ENST/ENSG) are recognized directly; for other
                # annotations, identifiers almost always contain digits while header
                # words like "transcript_id" do not. Without this check the first
                # mapping row would silently become the header.
                columns_are_ids = tx2gene_map.columns.astype(str).str.match(
                    r"^\w*(ENST|ENSG)\d"
                ).all()
                if columns_are_ids or _columns_look_like_identifiers(tx2gene_map.columns):
                    tx2gene_map = pd.read_csv(tx2gene_path, sep=separator, header=None)
                    self.logger.warning("tx2gene file has no header row; reading as headerless")
                if len(tx2gene_map.columns) < 2:
                    raise ValueError("tx2gene map must have at least 2 columns (transcript_id, gene_id)")
                tx2gene_map = tx2gene_map.iloc[:, :2].copy()
                # recover gene-first ordering (e.g. ENSG,ENST)
                first_column = tx2gene_map.iloc[:, 0].astype(str)
                second_column = tx2gene_map.iloc[:, 1].astype(str)
                if (
                    first_column.str.match(r"^ENSG\d").mean() > 0.5
                    and second_column.str.match(r"^\w*ENST").mean() > 0.5
                ):
                    tx2gene_map = tx2gene_map.iloc[:, [1, 0]]
                tx2gene_map.columns = ["transcript_id", "gene_id"]
                self.logger.warning("Renamed tx2gene columns to transcript_id and gene_id")

            if output_dir:
                tx2gene_file = Path(output_dir) / "tx2gene.tsv"
                tx2gene_file.parent.mkdir(parents=True, exist_ok=True)
                tx2gene_map.to_csv(tx2gene_file, sep="\t", index=False)
                self.logger.info(f"Saved tx2gene mapping to {tx2gene_file}")
        elif gtf_file is not None:
            tx2gene_map = self._create_tx2gene_from_gtf(gtf_file, output_dir)
            if output_dir:
                tx2gene_file = Path(output_dir) / "tx2gene.tsv"
        else:
            raise ValueError("Either gtf_file or tx2gene must be provided")

        return tx2gene_map, tx2gene_file

    def _find_quant_files(
        self,
        salmon_dir: str,
        sample_names: Optional[Sequence[str]] = None,
        sample_pattern: str = "quant.sf",
    ) -> Tuple[List[str], List[str]]:
        salmon_path = Path(salmon_dir)

        if sample_names is None:
            quant_files = sorted(salmon_path.rglob(sample_pattern))
            resolved_sample_names = [path.parent.name for path in quant_files]
            file_paths = [str(path) for path in quant_files]
        else:
            resolved_sample_names = []
            file_paths = []
            for sample_name in sample_names:
                quant_file = salmon_path / sample_name / sample_pattern
                if quant_file.exists():
                    resolved_sample_names.append(sample_name)
                    file_paths.append(str(quant_file))
                else:
                    self.logger.warning(f"Quantification file not found for sample {sample_name}")

        if not file_paths:
            raise FileNotFoundError(f"No files named {sample_pattern!r} found under {salmon_dir}")

        duplicated = sorted({n for n in resolved_sample_names if resolved_sample_names.count(n) > 1})
        if duplicated:
            raise ValueError(
                "Duplicate sample names from nested quant directories: " + ", ".join(duplicated)
            )

        return file_paths, resolved_sample_names

    def _to_dataframe(self, dataset, field: str, sample_names: Sequence[str]) -> pd.DataFrame:
        table = dataset[field].to_pandas()
        table.columns = list(sample_names)
        table.index.name = "gene_id"
        return table

    def build_gene_tables(
        self,
        salmon_dir: str,
        gtf_file: Optional[str] = None,
        tx2gene: Optional[str] = None,
        output_dir: Optional[str] = None,
        sample_names: Optional[Sequence[str]] = None,
        sample_pattern: str = "quant.sf",
        # Ensembl quant.sf carries versioned transcript IDs while GTF-derived tx2gene
        # maps do not; stripping versions is the only combination that matches out of the box
        ignore_transcript_version: bool = True,
    ) -> Dict[str, pd.DataFrame]:
        try:
            from pytximport import tximport
        except ImportError:
            raise ImportError("pytximport is not installed. Please install it with: pip install pytximport")

        file_paths, resolved_sample_names = self._find_quant_files(
            salmon_dir=salmon_dir,
            sample_names=sample_names,
            sample_pattern=sample_pattern,
        )
        tx2gene_map, _ = self._load_tx2gene_map(
            gtf_file=gtf_file,
            tx2gene=tx2gene,
            output_dir=output_dir,
        )

        self.logger.info(
            f"Ready for tximport: {len(tx2gene_map)} transcripts, {tx2gene_map['gene_id'].nunique()} genes"
        )
        self.logger.info(f"Running pytximport on {len(file_paths)} samples...")

        dataset = tximport(
            file_paths=file_paths,
            data_type="salmon",
            transcript_gene_map=tx2gene_map,
            counts_from_abundance="length_scaled_tpm",
            ignore_transcript_version=ignore_transcript_version,
            ignore_after_bar=False,
            output_type="xarray",
        )

        counts = self._to_dataframe(dataset, "counts", resolved_sample_names)
        counts.index.name = "gene_id"
        tpm = self._to_dataframe(dataset, "abundance", resolved_sample_names)
        tpm.index.name = "gene_id"
        log_tpm = np.log2(tpm + 1.0)
        log_tpm.index.name = "gene_id"

        return {
            "counts": counts,
            "tpm": tpm,
            "log2_tpm": log_tpm,
        }

    def export_gene_tables(
        self,
        salmon_dir: str,
        output_dir: str,
        gtf_file: Optional[str] = None,
        tx2gene: Optional[str] = None,
        sample_names: Optional[Sequence[str]] = None,
        sample_pattern: str = "quant.sf",
        ignore_transcript_version: bool = True,
    ) -> Dict[str, str]:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        tables = self.build_gene_tables(
            salmon_dir=salmon_dir,
            gtf_file=gtf_file,
            tx2gene=tx2gene,
            output_dir=output_dir,
            sample_names=sample_names,
            sample_pattern=sample_pattern,
            ignore_transcript_version=ignore_transcript_version,
        )

        outputs = {
            "gene_counts": output_path / "gene_counts.csv",
            "gene_tpm": output_path / "gene_tpm.csv",
            "gene_log2_tpm": output_path / "gene_log2_tpm.csv",
            "tx2gene": output_path / "tx2gene.tsv",
        }

        tables["counts"].to_csv(outputs["gene_counts"])
        tables["tpm"].to_csv(outputs["gene_tpm"])
        tables["log2_tpm"].to_csv(outputs["gene_log2_tpm"])

        return {name: str(path) for name, path in outputs.items()}

    @staticmethod
    def find_existing_gene_counts(salmon_dir: str) -> Optional[Path]:
        """Return a reusable gene counts file, or None if missing or stale.

        Stale means older than any quant.sf under salmon_dir: re-running quant
        must invalidate previously exported gene-level tables.
        """
        salmon_path = Path(salmon_dir)
        for candidate in ("gene_counts.csv", "gene_counts.tsv"):
            counts_path = salmon_path / candidate
            if not counts_path.exists():
                continue
            newest_quant = max(
                (q.stat().st_mtime for q in salmon_path.rglob("quant.sf")),
                default=0.0,
            )
            if counts_path.stat().st_mtime < newest_quant:
                logger.warning(
                    f"{counts_path} is older than the newest quant.sf in {salmon_dir}; "
                    "re-exporting gene-level tables"
                )
                return None
            return counts_path
        return None

    @staticmethod
    def find_reusable_gene_tables(
        salmon_dir: str,
        sample_names: Sequence[str],
    ) -> Optional[Dict[str, str]]:
        """Return existing gene-level tables when they are fresh for these samples.

        Re-importing every quant.sf is the only real work left on a fully
        resumed run, so the tables are reused when none of them predates the
        current samples' quant.sf files and the counts matrix covers exactly
        the requested samples.
        """
        salmon_path = Path(salmon_dir)
        outputs = {
            "gene_counts": salmon_path / "gene_counts.csv",
            "gene_tpm": salmon_path / "gene_tpm.csv",
            "gene_log2_tpm": salmon_path / "gene_log2_tpm.csv",
        }
        if not all(path.exists() for path in outputs.values()):
            return None

        quant_files = [salmon_path / name / "quant.sf" for name in sample_names]
        quant_files = [path for path in quant_files if path.exists()]
        if not quant_files:
            return None
        newest_quant = max(path.stat().st_mtime for path in quant_files)
        if any(path.stat().st_mtime < newest_quant for path in outputs.values()):
            return None

        try:
            header = pd.read_csv(outputs["gene_counts"], nrows=0)
        except (OSError, pd.errors.EmptyDataError):
            return None
        if set(str(column) for column in header.columns[1:]) != {
            str(name) for name in sample_names
        }:
            return None

        tx2gene = salmon_path / "tx2gene.tsv"
        if tx2gene.exists():
            outputs["tx2gene"] = tx2gene
        return {name: str(path) for name, path in outputs.items()}


def merge_salmon_quant_tables(
    salmon_dir: str,
    output_dir: str,
    gtf_file: Optional[str] = None,
    tx2gene: Optional[str] = None,
    sample_names: Optional[Sequence[str]] = None,
) -> Dict[str, str]:
    """Merge Salmon quant.sf files under a directory into gene-level tables.

    sample_names restricts the merge to those samples; None merges every
    quant.sf found (which picks up stale sample directories from older runs).
    """
    return SalmonExpressionExporter().export_gene_tables(
        salmon_dir=salmon_dir,
        output_dir=output_dir,
        gtf_file=gtf_file,
        tx2gene=tx2gene,
        sample_names=sample_names,
    )
