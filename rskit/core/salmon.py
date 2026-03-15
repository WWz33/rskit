from rskit.core.base import Tool
from pathlib import Path
from typing import Optional
from rskit.core.base import ToolBase
from rskit.config import SalmonConfig
from rskit.utils.validators import validate_file

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
