from pathlib import Path
from typing import Optional
from rskit.core.base import ToolBase, Tool
from rskit.config import StarConfig
from rskit.utils.validators import validate_file, check_star_index

class StarIndexer:
    def __init__(self, config: StarConfig):
        self.config = config
        self.tool = Tool("STAR")
        self.logger = self.tool.logger
    
    def build_index(self, genome_fasta: str, gtf_file: str, index_dir: str, force: bool = False) -> bool:
        validate_file(genome_fasta)
        validate_file(gtf_file)
        index_path = Path(index_dir)
        
        if index_path.exists() and check_star_index(index_dir):
            if not force:
                self.logger.info(f"STAR index already exists at {index_dir}, skipping")
                return True
            else:
                self.logger.info(f"Force rebuilding STAR index at {index_dir}")
        
        index_path.mkdir(parents=True, exist_ok=True)
        cmd = ["STAR", "--runThreadN", str(self.config.threads), "--runMode", "genomeGenerate",
               "--genomeDir", str(index_path), "--genomeFastaFiles", genome_fasta,
               "--sjdbGTFfile", gtf_file, "--sjdbOverhang", str(self.config.sjdb_overhang)]
        
        self.logger.info(f"Building STAR index in {index_dir}")
        return self.tool._run_command(cmd)

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
        
        # 检测输入文件格式
        read_cmd = 'zcat' if fq1.endswith('.gz') else 'cat'
        
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
        
        self.logger.info(f"Aligning {sample_name or 'sample'} with STAR")
        self.tool._run_command(cmd)
        
        return {
            "bam": f"{output_prefix}Aligned.out.bam",
            "transcriptome_bam": f"{output_prefix}Aligned.toTranscriptome.out.bam",
            "log": f"{output_prefix}Log.final.out"
        }
    
    def validate_inputs(self) -> bool:
        return self.tool._check_tool_installed()
