from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from rskit.core.base import Tool
from rskit.config import SalmonConfig
from rskit.utils.gtf import open as gtf_open
from rskit.utils.logger import get_logger
from rskit.utils.validators import validate_file

logger = get_logger(__name__)

class SalmonQuantifier:
    def __init__(self, config: SalmonConfig):
        self.config = config
        self.tool = Tool("salmon")
        self.logger = self.tool.logger
    
    def quantify(self, transcript_fasta: str, bam_file: str, output_dir: str,
                 sample_name: Optional[str] = None, skip_if_exists: bool = True) -> dict:
        validate_file(transcript_fasta)
        validate_file(bam_file)
        
        output_path = Path(output_dir)
        quant_file = output_path / "quant.sf"
        
        if skip_if_exists and quant_file.exists():
            self.logger.info(f"Quantification output already exists at {output_dir}, skipping")
            return {"quant": str(quant_file), "lib_format_counts": str(output_path / "lib_format_counts.json")}
        
        output_path.mkdir(parents=True, exist_ok=True)
        
        cmd = ["salmon", "quant", "-t", transcript_fasta, "-l", self.config.lib_type,
               "-a", bam_file, "-o", output_dir, "-p", str(self.config.threads)]
        
        if self.config.seq_bias:
            cmd.append("--seqBias")
        if self.config.gc_bias:
            cmd.append("--gcBias")
        if self.config.pos_bias:
            cmd.append("--posBias")
        if self.config.validate_mappings:
            cmd.append("--validateMappings")
        
        self.logger.info(f"Quantifying {sample_name or 'sample'} with Salmon")
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

        with open(gtf_file, "r", encoding="utf-8", errors="ignore") as reader:
            for rec in gtf_open(reader, "ensembl"):
                num_records += 1
                if rec.feature == "transcript" and rec.transcript_id and rec.gene_id:
                    tx2gene_map.setdefault(rec.transcript_id, rec.gene_id)

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
                if len(tx2gene_map.columns) < 2:
                    raise ValueError("tx2gene map must have at least 2 columns (transcript_id, gene_id)")
                tx2gene_map = tx2gene_map.iloc[:, :2].copy()
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
        ignore_transcript_version: bool = False,
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

        counts = self._to_dataframe(dataset, "counts", resolved_sample_names).T
        counts.index.name = "sample"
        tpm = self._to_dataframe(dataset, "abundance", resolved_sample_names).T
        tpm.index.name = "sample"
        log_tpm = np.log2(tpm + 1.0)
        log_tpm.index.name = "sample"

        return {
            "counts": counts,
            "tpm": tpm,
            "log2_tpm_plus1": log_tpm,
        }

    def export_gene_tables(
        self,
        salmon_dir: str,
        output_dir: str,
        gtf_file: Optional[str] = None,
        tx2gene: Optional[str] = None,
        sample_names: Optional[Sequence[str]] = None,
        sample_pattern: str = "quant.sf",
        ignore_transcript_version: bool = False,
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
            "gene_log2_tpm_plus1": output_path / "gene_log2_tpm_plus1.csv",
            "tx2gene": output_path / "tx2gene.tsv",
        }

        tables["counts"].to_csv(outputs["gene_counts"])
        tables["tpm"].to_csv(outputs["gene_tpm"])
        tables["log2_tpm_plus1"].to_csv(outputs["gene_log2_tpm_plus1"])

        return {name: str(path) for name, path in outputs.items()}

    @staticmethod
    def find_existing_gene_counts(salmon_dir: str) -> Optional[Path]:
        salmon_path = Path(salmon_dir)
        for candidate in ("gene_counts.csv", "gene_counts.tsv"):
            counts_path = salmon_path / candidate
            if counts_path.exists():
                return counts_path
        return None
