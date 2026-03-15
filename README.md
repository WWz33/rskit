# RNA-seq Tools

Python interface for STAR alignment and Salmon quantification.

## Installation

pip install -e .

## Usage

### Command Line

rnaseq-tools star-index genome.fa annotation.gtf index_dir --threads 56
rnaseq-tools star-align index_dir reads_1.fq reads_2.fq output_prefix --threads 56
rnaseq-tools salmon-quant transcripts.fa aligned.bam output_dir --threads 56

### Python API

from rnaseq_tools import RNAseqPipeline, PipelineConfig

config = PipelineConfig()
pipeline = RNAseqPipeline(config)

samples = {
    "sample1": {
        "fq1": "data/sample1_1.fq",
        "fq2": "data/sample1_2.fq"
    }
}

results = pipeline.run(
    samples=samples,
    genome_fasta="genome.fa",
    gtf_file="annotation.gtf",
    transcript_fasta="transcripts.fa",
    index_dir="STAR_index",
    output_dir="results"
)
