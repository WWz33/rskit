# rskit - RNA-seq Analysis Toolkit

A Python toolkit for RNA-seq analysis, providing command-line interface and Python API for common bioinformatics tasks including read alignment, quantification, differential expression analysis, and co-expression network analysis.

## Features

- **Quantification Pipeline**: STAR alignment + Salmon quantification
- **Differential Expression Analysis**: DESeq2-based analysis with automatic data processing
- **Co-expression Network Analysis**: WGCNA for gene module detection
- **Flexible Input Formats**: Automatic detection of CSV/TSV files
- **Unified `--coldata`**: Single sample metadata file works across all subcommands
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

**Coldata format** (CSV or TSV):

Basic format (quant only):
```csv
sample,r1,r2
sample1,sample1_R1.fq.gz,sample1_R2.fq.gz
sample2,sample2_R1.fq.gz,sample2_R2.fq.gz
```

Full format (compatible with quant, deseq2, and wgcna):
```csv
sample,id,condition,r1,r2
sample1,ctrl,control,sample1_R1.fq.gz,sample1_R2.fq.gz
sample2,ctrl,control,sample2_R1.fq.gz,sample2_R2.fq.gz
sample3,treat,treatment,sample3_R1.fq.gz,sample3_R2.fq.gz
sample4,treat,treatment,sample4_R1.fq.gz,sample4_R2.fq.gz
```

> **Tip**: Use the full format with a single `coldata.csv` for your entire workflow. Each subcommand reads only the columns it needs.

### 2. Differential Expression Analysis

Perform DESeq2 analysis on quantified data:

```bash
# From Salmon quantification directory (long options)
rskit deseq2 --salmon-dir ./03_quant --coldata coldata.csv --gtf annotation.gtf

# From Salmon quantification directory (short options)
rskit deseq2 -sd ./03_quant -S coldata.csv -gtf annotation.gtf

# From gene counts matrix
rskit deseq2 --gene-counts counts.csv --coldata coldata.csv
```

**Coldata format** (CSV or TSV):

Basic format (deseq2 only):
```csv
sample,id,condition
sample1,ctrl,control
sample2,ctrl,control
sample3,treat,treatment
sample4,treat,treatment
```

Full format (compatible with quant, deseq2, and wgcna):
```csv
sample,id,condition,r1,r2
sample1,ctrl,control,sample1_R1.fq.gz,sample1_R2.fq.gz
sample2,ctrl,control,sample2_R1.fq.gz,sample2_R2.fq.gz
sample3,treat,treatment,sample3_R1.fq.gz,sample3_R2.fq.gz
sample4,treat,treatment,sample4_R1.fq.gz,sample4_R2.fq.gz
```

> **Tip**: Use the full format to share one `coldata.csv` across quant, deseq2, and wgcna. Each subcommand reads only the columns it needs.

**Advanced options:**
```bash
# Multi-factor design (column names must exist in coldata)
rskit deseq2 -sd ./03_quant -S coldata.csv -gtf annotation.gtf --design "~batch + condition"

# Specify contrast
rskit deseq2 -sd ./03_quant -S coldata.csv -gtf annotation.gtf --contrast "condition,treatment,control"
```

> **Note**: The `--design` formula uses R-style syntax. Column names in the formula (e.g., `batch`, `condition`) must exist in your coldata file. PyDESeq2 will validate these column names when creating the design matrix.

### 3. WGCNA Co-expression Network Analysis

Perform weighted gene co-expression network analysis:

```bash
# Basic analysis (long options)
rskit wgcna --expression expression.csv --output-dir ./wgcna_results

# Basic analysis (short options)
rskit wgcna -e expression.csv -o ./wgcna_results

# With sample and gene metadata
rskit wgcna \
  -e expression.csv \
  -S sample_info.csv \
  -G gene_info.csv \
  -o ./wgcna_results

# Custom network parameters
rskit wgcna \
  -e expression.csv \
  -o ./wgcna_results \
  -nt signed \
  -min 30 \
  -p 6
```

**Input file formats:**

Expression matrix (CSV, samples x genes):
```csv
sample,GENE1,GENE2,GENE3
sample1,12.5,8.3,15.2
sample2,11.8,9.1,14.8
```

**Coldata format** (CSV or TSV, compatible with quant and deseq2):

Basic format (wgcna only):
```csv
sample,condition
sample1,control
sample2,treatment
```

Full format (compatible with quant, deseq2, and wgcna):
```csv
sample,id,condition,r1,r2
sample1,ctrl,control,sample1_R1.fq.gz,sample1_R2.fq.gz
sample2,treat,treatment,sample2_R1.fq.gz,sample2_R2.fq.gz
```

> **Tip**: Use the full format to share one `coldata.csv` across quant, deseq2, and wgcna.

Gene metadata (CSV):
```csv
gene_id,gene_name,gene_type
GENE1,BRAF,protein_coding
GENE2,TP53,protein_coding
```

## Command Reference

### rskit quant

Complete quantification pipeline (index → align → quant).

| Option | Short | Description |
|--------|-------|-------------|
| `--sample` | `-s` | Sample name (for single sample) |
| `--coldata` | `-S` | Sample sheet file (CSV/TSV) |
| `--r1` | `-1` | First read file |
| `--r2` | `-2` | Second read file |
| `--genome-fasta` | `-g` | Genome FASTA file |
| `--gtf-file` | `-gtf` | GTF annotation file |
| `--transcript-fasta` | `-gf` | Transcript FASTA file |
| `--output-dir` | `-o` | Output directory |
| `--index-dir` | | STAR index directory (default: STAR_index) |
| `--threads` | `-t` | Number of threads (default: 8) |
| `--trim` | | Trim reads with fastp |
| `--force-index` | | Force rebuild index |

### rskit deseq2

DESeq2 differential expression analysis.

| Option | Short | Description |
|--------|-------|-------------|
| `--salmon-dir` | `-sd` | Directory containing Salmon quant folders |
| `--gene-counts` | `-gc` | Gene counts matrix file |
| `--coldata` | `-S` | Sample metadata file (required) |
| `--gtf` | `-gtf` | GTF annotation file |
| `--tx2gene` | `-t2g` | Transcript-to-gene mapping file |
| `--design` | | Design formula (default: ~condition) |
| `--contrast` | | Contrast specification |
| `--alpha` | | Significance threshold (default: 0.05) |
| `-o, --output-dir` | `-o` | Output directory |
| `-t, --threads` | `-t` | Number of threads |

### rskit wgcna

WGCNA co-expression network analysis.

| Option | Short | Description |
|--------|-------|-------------|
| `-e, --expression` | `-e` | Expression matrix file (required) |
| `-o, --output-dir` | `-o` | Output directory (required) |
| `--coldata` | `-S` | Sample metadata file |
| `--gene-info` | `-G` | Gene metadata file |
| `--sep` | `-sep` | Separator for input files (default: ,) |
| `--name` | `-n` | Analysis name (default: WGCNA) |
| `--species` | `-s` | Species for enrichment analysis |
| `--level` | `-l` | Data level: gene, transcript (default: gene) |
| `--network-type` | `-nt` | Network type: unsigned, signed, signed hybrid (default: signed hybrid) |
| `--tom-type` | `-tom` | TOM type: unsigned, signed (default: signed) |
| `--min-module-size` | `-min` | Minimum module size (default: 50) |
| `--power` | `-p` | Soft thresholding power (auto-detected if not specified) |
| `--rsquared-cut` | `-rsquared` | R-squared cutoff (default: 0.9) |
| `--mean-cut` | `-mean` | Mean connectivity cutoff (default: 100) |
| `--mediss-thresh` | `-mediss` | Module merging threshold (default: 0.2) |
| `--tpm-cutoff` | `-tpm` | TPM cutoff for filtering (default: 1) |

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
    coldata="coldata.csv",
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
- **Gene metadata**: First column should contain gene identifiers matching the expression matrix columns

### Unified `--coldata` Parameter

All three subcommands use the same `--coldata` / `-S` parameter for sample metadata:

| 子命令 | coldata 用途 | 必需列 |
|--------|-------------|--------|
| quant | 样本文件（含 reads 路径） | `sample, r1, r2` |
| deseq2 | 样本元数据 | `sample, id, condition` |
| wgcna | 样本元数据 | `sample` + 任意元数据列 |

### Coldata Format Compatibility

**Full format** (compatible with all three subcommands):
```csv
sample,id,condition,r1,r2
sample1,ctrl,control,sample1_R1.fq.gz,sample1_R2.fq.gz
sample2,ctrl,control,sample2_R1.fq.gz,sample2_R2.fq.gz
sample3,treat,treatment,sample3_R1.fq.gz,sample3_R2.fq.gz
sample4,treat,treatment,sample4_R1.fq.gz,sample4_R2.fq.gz
```

**Usage scenarios**:

```bash
# 1. quant only (reads sample, r1, r2 columns)
rskit quant -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/

# 2. deseq2 only (reads sample, id, condition columns)
rskit deseq2 -sd ./03_quant -S coldata.csv -gtf annotation.gtf

# 3. wgcna only (reads sample + any metadata columns)
rskit wgcna -e expression.csv -S coldata.csv -o ./wgcna_results

# 4. Full pipeline: quant → deseq2 → wgcna (same coldata file for all)
rskit quant -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/
rskit deseq2 -sd results/03_quant -S coldata.csv -gtf annotation.gtf
rskit wgcna -e results/03_quant/gene_counts.csv -S coldata.csv -o ./wgcna_results
```

> **Key Point**: Each subcommand reads only the columns it needs from `--coldata`, so you can maintain a single coldata file for your entire analysis workflow.

## Citation

If you use rskit in your research, please cite:

```
[Your citation here]
```

## License

[Your license here]

## Contact

[Your contact information here]