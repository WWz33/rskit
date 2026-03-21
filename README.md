# rskit - RNA-seq Analysis Toolkit

A Python toolkit for RNA-seq analysis, providing command-line interface and Python API for common bioinformatics tasks including read alignment, quantification, differential expression analysis, and co-expression network analysis.

## Features

- **Quantification Pipeline**: STAR alignment + Salmon quantification
- **Differential Expression Analysis**: DESeq2-based analysis with automatic data processing
- **Co-expression Network Analysis**: WGCNA for gene module detection
- **Complete Pipeline**: One command from reads to differential expression (`rskit all`)
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

### 1. Complete Pipeline (quant → deseq2)

Run the complete analysis in one command:

```bash
rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/

# With custom design formula
rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ --design "~batch + condition"
```

### 2. Quantification Pipeline

Run complete quantification from raw reads to gene counts:

```bash
# Single sample
rskit quant -s sample1 -1 sample1_R1.fq.gz -2 sample1_R2.fq.gz -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/

# Multiple samples using coldata
rskit quant -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/
```

### 3. Differential Expression Analysis

```bash
# From Salmon quantification directory
rskit deseq2 -sd ./03_quant -S coldata.csv -gtf annotation.gtf

# From gene counts matrix
rskit deseq2 -gc counts.csv -S coldata.csv

# Multi-factor design
rskit deseq2 -sd ./03_quant -S coldata.csv -gtf annotation.gtf --design "~batch + condition"
```

### 4. WGCNA Co-expression Network Analysis

```bash
# Basic analysis
rskit wgcna -e expression.csv -o ./wgcna_results

# With metadata
rskit wgcna -e expression.csv -S coldata.csv -G gene_info.csv -o ./wgcna_results
```

## Coldata Format

All subcommands use the same `--coldata` / `-S` parameter. Use one file for your entire workflow:

**Full format** (compatible with all subcommands):
```csv
sample,id,condition,r1,r2
sample1,ctrl,control,sample1_R1.fq.gz,sample1_R2.fq.gz
sample2,ctrl,control,sample2_R1.fq.gz,sample2_R2.fq.gz
sample3,treat,treatment,sample3_R1.fq.gz,sample3_R2.fq.gz
sample4,treat,treatment,sample4_R1.fq.gz,sample4_R2.fq.gz
```

Each subcommand reads only the columns it needs:
- **quant**: reads `sample`, `r1`, `r2`
- **deseq2**: reads `sample`, `id`, `condition`
- **wgcna**: reads `sample` + any metadata columns

## Command Reference

### rskit all

Complete pipeline: quantification + DESeq2 analysis.

| Option | Description |
|--------|-------------|
| `-S, --coldata` | Coldata file (CSV/TSV) with columns: sample,id,condition,r1,r2 (required) |
| `-g, --genome-fasta` | Genome FASTA file (required) |
| `-gtf, --gtf-file` | GTF annotation file (required) |
| `-gf, --transcript-fasta` | Transcript FASTA file (required) |
| `-o, --output-dir` | Output directory (required) |
| `--index-dir` | STAR index directory (default: STAR_index) |
| `-t2g, --tx2gene` | Transcript-to-gene mapping file |
| `-t, --threads` | Number of threads (default: 8) |
| `--trim` | Trim reads with fastp |
| `--force-index` | Force rebuild index |
| `--design` | Design formula (default: ~condition) |
| `--contrast` | Contrast specification (e.g., 'condition,treatment,control') |
| `--alpha` | Significance threshold (default: 0.05) |

### rskit quant

Complete quantification pipeline (index → align → quant).

| Option | Description |
|--------|-------------|
| `-s, --sample` | Sample name (for single sample) |
| `-S, --coldata` | Sample file (CSV/TSV) with columns: sample,r1,r2 |
| `-1, --r1` | First read file |
| `-2, --r2` | Second read file |
| `-g, --genome-fasta` | Genome FASTA file (required) |
| `-gtf, --gtf-file` | GTF annotation file (required) |
| `-gf, --transcript-fasta` | Transcript FASTA file (required) |
| `-o, --output-dir` | Output directory (required) |
| `--index-dir` | STAR index directory (default: STAR_index) |
| `-t, --threads` | Number of threads (default: 8) |
| `--trim` | Trim reads with fastp |
| `--force-index` | Force rebuild index |

### rskit deseq2

DESeq2 differential expression analysis.

| Option | Description |
|--------|-------------|
| `-sd, --salmon-dir` | Directory containing Salmon quant folders |
| `-gc, --gene-counts` | Gene counts matrix file |
| `-S, --coldata` | Sample metadata file (required) |
| `-gtf, --gtf` | GTF annotation file |
| `-t2g, --tx2gene` | Transcript-to-gene mapping file |
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
| `-S, --coldata` | Sample metadata file |
| `-G, --gene-info` | Gene metadata file |
| `-sep, --sep` | Separator for input files (default: ,) |
| `-n, --name` | Analysis name (default: WGCNA) |
| `-s, --species` | Species for enrichment analysis |
| `-l, --level` | Data level: gene, transcript (default: gene) |
| `-nt, --network-type` | Network type: unsigned, signed, signed hybrid (default: signed hybrid) |
| `-tom, --tom-type` | TOM type: unsigned, signed (default: signed) |
| `-min, --min-module-size` | Minimum module size (default: 50) |
| `-p, --power` | Soft thresholding power (auto-detected if not specified) |
| `-rsquared, --rsquared-cut` | R-squared cutoff (default: 0.9) |
| `-mean, --mean-cut` | Mean connectivity cutoff (default: 100) |
| `-mediss, --mediss-thresh` | Module merging threshold (default: 0.2) |
| `-tpm, --tpm-cutoff` | TPM cutoff for filtering (default: 1) |

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

### Complete Pipeline Output (rskit all)
```
results/
├── 00_index/          # STAR index
├── 01_clean_data/     # Trimmed reads (if --trim)
├── 02_bam/           # Aligned BAM files
├── 03_quant/         # Salmon quantification
└── 04_deseq2/        # DESeq2 results
```

### DESeq2 Output
```
04_deseq2/
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

## Citation

If you use rskit in your research, please cite:

```
[Your citation here]
```

## License

[Your license here]

## Contact

[Your contact information here]