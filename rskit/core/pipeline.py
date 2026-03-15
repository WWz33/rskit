from pathlib import Path
from typing import Dict
from rskit.config import PipelineConfig
from rskit.core.star import StarIndexer, StarAligner
from rskit.core.salmon import SalmonQuantifier
from rskit.utils.logger import get_logger
from rskit.utils.validators import check_star_index

logger = get_logger(__name__)

class RNAseqPipeline:
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.indexer = StarIndexer(config.star)
        self.aligner = StarAligner(config.star)
        self.quantifier = SalmonQuantifier(config.salmon)
        self.logger = logger

    def run(self, samples: Dict[str, Dict], genome_fasta: str, gtf_file: str,
            transcript_fasta: str, index_dir: str, output_dir: str,
            quant_output_dir: str, force_index: bool = False) -> Dict:
        results = {}
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        quant_path = Path(quant_output_dir)
        quant_path.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"Checking STAR index at {index_dir}")
        if not Path(index_dir).exists() or not check_star_index(index_dir):
            self.logger.info("Index not found, building...")
            self.indexer.build_index(genome_fasta, gtf_file, index_dir, force=force_index)
        else:
            self.logger.info("Index found, skipping build")

        for sample_name, sample_data in samples.items():
            self.logger.info(f"Processing sample: {sample_name}")
            sample_output = output_path / sample_name
            sample_output.mkdir(parents=True, exist_ok=True)

            align_prefix = str(sample_output / f"{sample_name}_")
            align_results = self.aligner.align(index_dir, sample_data["fq1"], sample_data["fq2"],
                                              align_prefix, sample_name=sample_name)

            salmon_output = quant_path / sample_name
            quant_results = self.quantifier.quantify(transcript_fasta, align_results["transcriptome_bam"],
                                                     str(salmon_output), sample_name=sample_name)

            results[sample_name] = {"alignment": align_results, "quantification": quant_results}
            self.logger.info(f"Completed sample: {sample_name}")

        return results
