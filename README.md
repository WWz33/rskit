## Installation
```
pip install -e .
```
## Usage
```
rskit --help
```
### Command Line
```
rskit quant -s ${i} -1 ../${i}_1 -2 ../${i}_2 --index-dir ../STAR_Gmax_index/  -t 96 -g ../Gmax_508_Wm82.a4.v1.fa -gtf ../Gmax_508_Wm82.a4.v1.gene_exons.gff3.gtf -gf ../gffread_transcript.fa  -o data/
rskit quant -S sample.tsv 
|  sample | r1          | r2         |
| ------: | ----------- | ---------- |
| sample1 | s11_1.fq    | s1_1.fq    |
| sample2 | s2_r1.fq.gz | 2_r1.fq.gz |
```
### Python API
```
from rskit import RNAseqPipeline, PipelineConfig

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
```