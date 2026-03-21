# rskit - RNA-seq Analysis Toolkit

A Python toolkit for RNA-seq analysis, providing command-line interface and Python API for common bioinformatics tasks including read alignment, quantification, differential expression analysis, and co-expression network analysis.

## Features

- **Quantification Pipeline**: STAR alignment + Salmon quantification
- **Differential Expression Analysis**: DESeq2-based analysis with automatic data processing
- **Co-expression Network Analysis**: WGCNA for gene module detection
- **Flexible Input Formats**: Automatic detection of CSV/TSV files
- **Command-line Interface**: Easy-to-use CLI for all analyses
- **Python API**: Programmatic access for custom workflows

## Installation

### Install from source
```bash
git clone https://github.com/your-repo/rskit.git
cd rskit
pip install -e .
```

### Dependencies
- Python >= 3.8
- pandas
- numpy
- pydeseq2
- PyWGCNA
- STAR (for alignment)
- Salmon (for quantification)
- fastp (optional, for trimming)

## Quick Start

### 1. Quantification Pipeline

Run complete quantification from raw reads to gene counts:

```bash
# Single sample
rskit quant \
  -s sample1 \
  -1 sample1_R1.fq.gz \
  -2 sample1_R2.fq.gz \
  -g genome.fa \
  -gtf annotation.gtf \
  -gf transcripts.fa \
  -o results/

# Multiple samples using a sample sheet
rskit quant -S samples.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/
```

**Sample sheet format** (CSV or TSV):

Basic format (quant only):
```csv
sample,r1,r2
sample1,sample1_R1.fq.gz,sample1_R2.fq.gz
sample2,sample2_R1.fq.gz,sample2_R2.fq.gz
```

Extended format (compatible with both quant and deseq2):
```csv
sample,id,condition,r1,r2
sample1,ctrl,control,sample1_R1.fq.gz,sample1_R2.fq.gz
sample2,ctrl,control,sample2_R1.fq.gz,sample2_R2.fq.gz
sample3,treat,treatment,sample3_R1.fq.gz,sample3_R2.fq.gz
sample4,treat,treatment,sample4_R1.fq.gz,sample4_R2.fq.gz
```

> **Note**: The extended format with `id` and `condition` columns can be used directly as coldata for deseq2 analysis after quantification.

### 2. Differential Expression Analysis

Perform DESeq2 analysis on quantified data:

```bash
# From Salmon quantification directory
rskit deseq2 --salmon-dir ./03_quant --coldata coldata.csv --gtf annotation.gtf

# From gene counts matrix
rskit deseq2 --gene-counts counts.csv --coldata coldata.csv
```

**Coldata format** (CSV):

Basic format (deseq2 only):
```csv
sample,id,condition
sample1,ctrl,control
sample2,ctrl,control
sample3,treat,treatment
sample4,treat,treatment
```

> **Tip**: Use the same CSV file for both quant and deseq2 by including `r1` and `r2` columns:
> ```csv
> sample,id,condition,r1,r2
> sample1,ctrl,control,sample1_R1.fq.gz,sample1_R2.fq.gz
> sample2,ctrl,control,sample2_R1.fq.gz,sample2_R2.fq.gz
> sample3,treat,treatment,sample3_R1.fq.gz,sample3_R2.fq.gz
> sample4,treat,treatment,sample4_R1.fq.gz,sample4_R2.fq.gz
> ```

**Advanced options:**
```bash
# Multi-factor design
rskit deseq2 --salmon-dir ./03_quant --coldata coldata.csv --gtf annotation.gtf --design "~batch + condition"

# Specify contrast
rskit deseq2 --salmon-dir ./03_quant --coldata coldata.csv --gtf annotation.gtf --contrast "condition,treatment,control"
```

### 3. WGCNA Co-expression Network Analysis

Perform weighted gene co-expression network analysis:

```bash
# Basic analysis
rskit wgcna --expression expression.csv --output-dir ./wgcna_results

# With sample and gene metadata
rskit wgcna \
  --expression expression.csv \
  --sample-info sample_info.csv \
  --gene-info gene_info.csv \
  --output-dir ./wgcna_results

# Custom network parameters
rskit wgcna \
  --expression expression.csv \
  --output-dir ./wgcna_results \
  --network-type signed \
  --min-module-size 30 \
  --power 6
```

**Input file formats:**

Expression matrix (CSV, samples x genes):
```csv
sample,GENE1,GENE2,GENE3
sample1,12.5,8.3,15.2
sample2,11.8,9.1,14.8
```

Sample metadata (CSV):

Basic format:
```csv
sample,condition
sample1,control
sample2,treatment
```

Extended format (compatible with deseq2):
```csv
sample,id,condition
sample1,ctrl,control
sample2,treat,treatment
```

> **Note**: The sample metadata format is compatible with deseq2 coldata. You can use the same file for both analyses.

Gene metadata (CSV):
```csv
gene_id,gene_name,gene_type
GENE1,BRAF,protein_coding
GENE2,TP53,protein_coding
```

## Command Reference

### rskit quant

Complete quantification pipeline (index → align → quant).

| Option | Description |
|--------|-------------|
| `-s, --sample` | Sample name (for single sample) |
| `-S, --sample-tsv` | Sample sheet file (CSV/TSV) |
| `-1, --r1` | First read file |
| `-2, --r2` | Second read file |
| `-g, --genome-fasta` | Genome FASTA file |
| `-gtf, --gtf-file` | GTF annotation file |
| `-gf, --transcript-fasta` | Transcript FASTA file |
| `-o, --output-dir` | Output directory |
| `--index-dir` | STAR index directory (default: STAR_index) |
| `-t, --threads` | Number of threads (default: 8) |
| `--trim` | Trim reads with fastp |
| `--force-index` | Force rebuild index |

### rskit deseq2

DESeq2 differential expression analysis.

| Option | Description |
|--------|-------------|
| `--salmon-dir` | Directory containing Salmon quant folders |
| `--gene-counts` | Gene counts matrix file |
| `--coldata` | Sample metadata file (required) |
| `--gtf` | GTF annotation file |
| `--tx2gene` | Transcript-to-gene mapping file |
| `--design` | Design formula (default: ~condition) |
| `--contrast` | Contrast specification |
| `--alpha` | Significance threshold (default: 0.05) |
| `-o, --output-dir` | Output directory |
| `-t, --threads` | Number of threads |

### rskit wgcna

WGCNA co-expression network analysis.

| Option | Description |
|--------|-------------|
| `-e, --expression` | Expression matrix file (required) |
| `-o, --output-dir` | Output directory (required) |
| `--sample-info` | Sample metadata file |
| `--gene-info` | Gene metadata file |
| `--name` | Analysis name (default: WGCNA) |
| `--network-type` | Network type: unsigned, signed, signed hybrid (default: signed hybrid) |
| `--tom-type` | TOM type: unsigned, signed (default: signed) |
| `--min-module-size` | Minimum module size (default: 50) |
| `--power` | Soft thresholding power (auto-detected if not specified) |
| `--rsquared-cut` | R-squared cutoff (default: 0.9) |
| `--mean-cut` | Mean connectivity cutoff (default: 100) |
| `--mediss-thresh` | Module merging threshold (default: 0.2) |
| `--tpm-cutoff` | TPM cutoff for filtering (default: 1) |
| `--species` | Species for enrichment analysis |
| `--level` | Data level: gene, transcript (default: gene) |

## Python API

### Quantification Pipeline

```python
from rskit import RNAseqPipeline, PipelineConfig

config = PipelineConfig()
pipeline = RNAseqPipeline(config)

samples = {
    "sample1": {
        "fq1": "data/sample1_R1.fq",
        "fq2": "data/sample1_R2.fq"
    }
}

results = pipeline.run(
    samples=samples,
    genome_fasta="genome.fa",
    gtf_file="annotation.gtf",
    transcript_fasta="transcripts.fa",
    index_dir="STAR_index",
    output_dir="results/02_bam",
    quant_output_dir="results/03_quant"
)
```

### DESeq2 Analysis

```python
from rskit.core.deseq2 import Deseq2Analyzer
from rskit.config import DESeq2Config

config = DESeq2Config()
analyzer = Deseq2Analyzer(config)

# Load data
counts_df = analyzer.load_counts_from_file("counts.csv")
metadata_df = analyzer.load_metadata("coldata.csv")

# Run analysis
results_df = analyzer.analyze(
    counts_df=counts_df,
    metadata_df=metadata_df,
    contrast=["condition", "treatment", "control"]
)
```

### WGCNA Analysis

```python
from rskit.core.wgcna import WGCNAAnalyzer

analyzer = WGCNAAnalyzer(
    output_dir="./wgcna_results",
    name="MyWGCNA",
    network_type="signed hybrid",
    min_module_size=50
)

# Load data
analyzer.load_data(
    expression_file="expression.csv",
    sample_info_file="sample_info.csv",
    gene_info_file="gene_info.csv"
)

# Run analysis
wgcna_obj = analyzer.run_analysis()

# Save results
analyzer.save_results()
```

## Output Structure

### Quantification Output
```
results/
├── 00_index/          # STAR index
├── 01_clean_data/     # Trimmed reads (if --trim)
├── 02_bam/           # Aligned BAM files
└── 03_quant/         # Salmon quantification
```

### DESeq2 Output
```
deseq2_results/
├── deseq2_results.csv      # Differential expression results
├── normalized_counts.csv   # Normalized counts
├── pca_plot.pdf           # PCA visualization
├── volcano_plot.pdf       # Volcano plot
└── ma_plot.pdf            # MA plot
```

### WGCNA Output
```
wgcna_results/
├── figures/              # All generated figures
├── WGCNA.p              # Pickled WGCNA object
└── module_info.csv       # Module membership
```

## File Format Notes

- **Automatic separator detection**: All subcommands automatically detect CSV (comma-separated) and TSV (tab-separated) files based on file extension
- **Gene expression matrix**: Rows = samples, Columns = genes
- **Sample metadata**: First column should contain sample identifiers matching the expression matrix
- **Gene metadata**: First column should contain gene identifiers matching the expression matrix columns
- **Unified sample sheet**: Use a single CSV with columns `sample,id,condition,r1,r2` for both quant and deseq2 workflows
- **Cross-tool compatibility**: Sample metadata files (coldata) are compatible between deseq2 and wgcna - use `sample,id,condition` format for both

## Citation

If you use rskit in your research, please cite:

```
[Your citation here]
```

## License

[Your license here]

## Contact

[Your contact information here]