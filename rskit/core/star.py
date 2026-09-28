import json
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional
from rskit.cli_args import merge_extra_args
from rskit.core.base import ToolBase, Tool
from rskit.config import StarConfig
from rskit.utils.validators import validate_file, check_star_index

STAR_INDEX_PROTECTED_OPTIONS = {
    "--runThreadN",
    "--runMode",
    "--genomeDir",
    "--genomeFastaFiles",
    "--sjdbGTFfile",
    # not managed by rskit, but kept consistent with the alignment guard: the
    # index build writes its Log files through this prefix
    "--outFileNamePrefix",
}

STAR_ALIGN_PROTECTED_OPTIONS = {
    "--runMode",
    "--runThreadN",
    "--genomeDir",
    "--readFilesIn",
    "--readFilesCommand",
    "--outFileNamePrefix",
    "--outSAMtype",
    "--quantMode",
}

# fingerprint of the inputs a STAR index was built from; lets rskit warn when
# the genome FASTA or GTF changed after the index was built (STAR itself would
# silently align against the stale index)
INDEX_FINGERPRINT_NAME = ".rskit_index.json"


def _file_fingerprint(path) -> Dict[str, object]:
    file_path = Path(path)
    stat = file_path.stat()
    return {"path": str(file_path.resolve()), "size": stat.st_size, "mtime": stat.st_mtime}


def _star_version() -> Optional[str]:
    try:
        result = subprocess.run(["STAR", "--version"], capture_output=True, text=True, check=True)
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip() or None


def read_index_fingerprint(index_dir) -> Optional[dict]:
    fingerprint_path = Path(index_dir) / INDEX_FINGERPRINT_NAME
    if not fingerprint_path.exists():
        return None
    try:
        return json.loads(fingerprint_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def index_input_changes(index_dir, genome_fasta, gtf_file) -> Optional[List[str]]:
    """Compare current inputs against the fingerprint recorded at index build.

    None: the index carries no fingerprint (built by hand or by an older
    rskit) and cannot be verified. []: inputs unchanged. Otherwise the names
    of the inputs that changed since the build.
    """
    fingerprint = read_index_fingerprint(index_dir)
    if fingerprint is None:
        return None
    changed = []
    for key, path in (("genome_fasta", genome_fasta), ("gtf_file", gtf_file)):
        try:
            current = _file_fingerprint(path)
        except OSError:
            # input vanished; alignment only needs the index, so surface it
            # as a change rather than failing here
            changed.append(key)
            continue
        if fingerprint.get(key) != current:
            changed.append(key)
    return changed

class StarIndexer:
    def __init__(self, config: StarConfig):
        self.config = config
        self.tool = Tool("STAR")
        self.logger = self.tool.logger
    
    def build_index(self, genome_fasta: str, gtf_file: str, index_dir: str, force: bool = False) -> bool:
        validate_file(genome_fasta)
        validate_file(gtf_file)
        index_path = Path(index_dir)
        
        if index_path.exists() and check_star_index(index_dir) and not force:
            self.logger.info(f"STAR index already exists at {index_dir}, skipping")
            return True

        # STAR genomeGenerate refuses a non-empty --genomeDir, so clear stale/partial indexes first
        if index_path.exists():
            self.logger.info(f"Clearing existing index directory {index_dir}")
            shutil.rmtree(index_path)

        index_path.mkdir(parents=True, exist_ok=True)
        cmd = ["STAR", "--runThreadN", str(self.config.threads), "--runMode", "genomeGenerate",
               "--genomeDir", str(index_path), "--genomeFastaFiles", genome_fasta,
               "--sjdbGTFfile", gtf_file, "--sjdbOverhang", str(self.config.sjdb_overhang),
               "--outFileNamePrefix", str(index_path) + "/"]
        cmd = merge_extra_args(cmd, self.config.extra_args, STAR_INDEX_PROTECTED_OPTIONS)
        
        self.logger.info(f"Building STAR index in {index_dir}")
        built = self.tool._run_command(cmd)
        self._write_fingerprint(index_path, genome_fasta, gtf_file)
        return built

    def _write_fingerprint(self, index_path: Path, genome_fasta: str, gtf_file: str) -> None:
        """Record the build inputs so later runs can detect annotation changes."""
        fingerprint = {
            "genome_fasta": _file_fingerprint(genome_fasta),
            "gtf_file": _file_fingerprint(gtf_file),
            "sjdb_overhang": self.config.sjdb_overhang,
            "star_version": _star_version(),
        }
        fingerprint_path = index_path / INDEX_FINGERPRINT_NAME
        fingerprint_path.write_text(json.dumps(fingerprint, indent=2), encoding="utf-8")
        self.logger.info(f"Recorded index build inputs in {fingerprint_path}")

class StarAligner:
    def __init__(self, config: StarConfig):
        self.config = config
        self.tool = Tool("STAR")
        self.logger = self.tool.logger
    
    def align(self, index_dir: str, fq1: str, fq2: str, output_prefix: str, 
              sample_name: Optional[str] = None, auto_index: bool = False,
              genome_fasta: Optional[str] = None, gtf_file: Optional[str] = None) -> dict:
        validate_file(fq1)
        validate_file(fq2)
        
        if not Path(index_dir).exists() or not check_star_index(index_dir):
            if auto_index and genome_fasta and gtf_file:
                self.logger.info(f"Index not found, auto-creating at {index_dir}")
                indexer = StarIndexer(self.config)
                indexer.build_index(genome_fasta, gtf_file, index_dir)
            else:
                raise FileNotFoundError(f"STAR index not found at {index_dir}")
        
        output_path = Path(output_prefix).parent
        output_path.mkdir(parents=True, exist_ok=True)
        
        # 检测输入文件格式；两段 read 压缩格式必须一致，否则 STAR 会对纯文本跑解压命令
        if fq1.endswith(".gz") != fq2.endswith(".gz"):
            raise ValueError(f"Mates must both be gzipped or both plain: {fq1}, {fq2}")
        # gzip -dc 比 zcat 更可移植（部分系统只有 gzip 没有 zcat）
        read_cmd = 'gzip -dc' if fq1.endswith('.gz') else 'cat'
        
        cmd = ["STAR", "--runThreadN", str(self.config.threads), "--genomeDir", index_dir,
               "--readFilesIn", fq1, fq2, "--readFilesCommand", read_cmd,
               "--outFileNamePrefix", output_prefix, "--outSAMtype", "BAM", "Unsorted",
               "--quantMode", "TranscriptomeSAM", "--twopassMode", self.config.two_pass_mode,
               "--outSAMunmapped", self.config.out_sam_unmapped, "--outFilterType", self.config.out_filter_type,
               "--quantTranscriptomeSAMoutput", self.config.quant_transcriptome_sam_output,
               "--outSAMattributes", "NH", "HI", "AS", "nM", "NM", "MD", "jM", "jI",
               "--sjdbOverhang", str(self.config.sjdb_overhang),
               "--alignIntronMin", str(self.config.align_intron_min),
               "--alignIntronMax", str(self.config.align_intron_max),
               "--alignMatesGapMax", str(self.config.align_mates_gap_max),
               "--alignSJoverhangMin", str(self.config.align_sj_overhang_min),
               "--outFilterMismatchNoverReadLmax", str(self.config.out_filter_mismatch_n_over_read_lmax),
               "--outFilterMismatchNmax", str(self.config.out_filter_mismatch_nmax),
               "--outFilterMultimapNmax", str(self.config.out_filter_multimap_nmax),
               "--alignSJDBoverhangMin", str(self.config.align_sjdb_overhang_min)]
        cmd = merge_extra_args(cmd, self.config.extra_args, STAR_ALIGN_PROTECTED_OPTIONS)
        
        self.logger.info(f"Aligning {sample_name or 'sample'} with STAR")
        self.tool._run_command(cmd)
        
        return {
            "bam": f"{output_prefix}Aligned.out.bam",
            "transcriptome_bam": f"{output_prefix}Aligned.toTranscriptome.out.bam",
            "log": f"{output_prefix}Log.final.out"
        }
    
    def validate_inputs(self) -> bool:
        return self.tool._check_tool_installed()
