# rskit - RNA-seq Analysis Toolkit

<!-- README-I18N:START -->

**English** | [中文](./README.zh.md)

<!-- README-I18N:END -->

A Python toolkit for RNA-seq analysis with a CLI and Python API for common workflows including read alignment, Salmon quantification, DESeq2 differential expression, and WGCNA.

## Features

- Quantification pipeline: STAR alignment + Salmon quantification
- Gene-level expression export during `quant`
- Differential expression analysis with DESeq2
- Co-expression network analysis with WGCNA
- End-to-end workflow with `rskit all`
- Automatic CSV/TSV detection
- Shared `--coldata` metadata format across subcommands
- Early `tx2gene.tsv` and gene-level matrix export before DESeq2
- Strict sample ID validation between metadata and count/expression matrices
- Input preflight checks with `rskit validate` / `rskit doctor`

## Installation

### Install from source
```bash
git clone https://github.com/WWz33/rskit.git
cd rskit
pip install -e .
```

### Dependencies

- Python >= 3.8
- pandas
- numpy
- pydeseq2
- pytximport
- PyWGCNA
- STAR
- Salmon
- fastp (optional)

## Quick Start

### 1. Complete pipeline

```bash
rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/

# With custom design formula
rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ --design "~batch + condition"
```

### 2. Quantification pipeline

`quant` writes per-sample Salmon output plus `tx2gene.tsv`, `gene_counts.csv`, `gene_tpm.csv`, and `gene_log2_tpm_plus1.csv` into `03_quant/`. The transcript-to-gene map is created during quantification, before DESeq2 needs gene-level counts.

```bash
# Single sample
rskit quant -s sample1 -1 sample1_R1.fq.gz -2 sample1_R2.fq.gz -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/

# Multiple samples using coldata
rskit quant -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/
```

### 3. Differential expression analysis

```bash
# Prefer precomputed gene counts from a quant directory
rskit deseq2 -sd ./03_quant -S coldata.csv -gtf annotation.gtf

# Or provide a counts matrix directly
rskit deseq2 -gc counts.csv -S coldata.csv

# Multi-factor design
rskit deseq2 -sd ./03_quant -S coldata.csv -gtf annotation.gtf --design "~batch + condition"
```

### 4. WGCNA

```bash
rskit wgcna -e expression.csv -o ./wgcna_results
rskit wgcna -e expression.csv -S coldata.csv -G gene_info.csv -o ./wgcna_results
```

### 5. Input validation

```bash
# Check metadata columns, read paths, and counts/sample alignment before running analysis
rskit validate -S coldata.csv --check-reads -gc counts.csv --design "~batch + condition"

# doctor is an alias for validate
rskit doctor -S coldata.csv -e expression.csv
```

## Coldata Format

All subcommands use the same `--coldata` / `-S` parameter.

```csv
sample,id,condition,r1,r2
sample1,ctrl,control,sample1_R1.fq.gz,sample1_R2.fq.gz
sample2,ctrl,control,sample2_R1.fq.gz,sample2_R2.fq.gz
sample3,treat,treatment,sample3_R1.fq.gz,sample3_R2.fq.gz
sample4,treat,treatment,sample4_R1.fq.gz,sample4_R2.fq.gz
```

Each subcommand reads only the columns it needs:

- `quant`: `sample`, `r1`, `r2`
- `deseq2`: `sample` plus every metadata column referenced by `--design`; the default `~condition` requires `condition`
- `wgcna`: `sample` plus any metadata columns

Relative `r1` and `r2` paths are resolved relative to the coldata file. Count and expression matrices are expected as rows = samples and columns = genes. If a DESeq2 counts file is genes x samples, it is transposed only when the sample IDs match the coldata columns. After orientation, sample IDs must match coldata exactly; rskit does not silently drop or intersect samples.

## Command Reference

### `rskit all`

Complete pipeline: quantification + DESeq2 analysis.

| Option | Description |
|--------|-------------|
| `-S, --coldata` | Coldata file with `sample,id,condition,r1,r2` |
| `-g, --genome-fasta` | Genome FASTA file |
| `-gtf, --gtf-file` | GTF annotation file |
| `-gf, --transcript-fasta` | Transcript FASTA file |
| `-o, --output-dir` | Output directory |
| `-idx, --index-dir` | STAR index directory |
| `-t2g, --tx2gene` | Existing transcript-to-gene mapping file; otherwise `tx2gene.tsv` is created from `--gtf-file` during quantification |
| `-t, --threads` | Threads per sample |
| `-p, --parallel` | Total cores for parallel execution |
| `--trim` | Trim reads with fastp |
| `--force-index` | Force STAR index rebuild |
| `--skip-existing` | Skip sample work when output already exists |
| `--design` | DESeq2 design formula |
| `--contrast` | Contrast specification |
| `--alpha` | Adjusted p-value threshold |
| `--lfc` | Log2 fold-change threshold |

### `rskit quant`

Complete quantification pipeline: index -> align -> quant -> gene-level table export.

| Option | Description |
|--------|-------------|
| `-s, --sample` | Sample name for single-sample mode |
| `-S, --coldata` | Sample file with `sample,r1,r2` |
| `-1, --r1` | First read file |
| `-2, --r2` | Second read file |
| `-g, --genome-fasta` | Genome FASTA file |
| `-gtf, --gtf-file` | GTF annotation file |
| `-gf, --transcript-fasta` | Transcript FASTA file |
| `-o, --output-dir` | Output directory |
| `-idx, --index-dir` | STAR index directory |
| `-t2g, --tx2gene` | Existing transcript-to-gene mapping file for gene-level export; otherwise `tx2gene.tsv` is created from `--gtf-file` |
| `-t, --threads` | Threads per sample |
| `-p, --parallel` | Total cores for parallel execution |
| `--trim` | Trim reads with fastp |
| `--force-index` | Force STAR index rebuild |
| `--skip-existing` | Skip sample work when output already exists |

### `rskit deseq2`

DESeq2 differential expression analysis.

| Option | Description |
|--------|-------------|
| `-sd, --salmon-dir` | Directory containing Salmon quant folders |
| `-gc, --gene-counts` | Gene counts matrix file |
| `-S, --coldata` | Sample metadata file with `sample` and all columns referenced by `--design` |
| `-gtf, --gtf` | GTF annotation file |
| `-t2g, --tx2gene` | Transcript-to-gene mapping file |
| `--design` | Design formula |
| `--contrast` | Contrast specification |
| `--alpha` | Adjusted p-value threshold |
| `--lfc` | Log2 fold-change threshold |
| `-o, --output-dir` | Output directory |
| `-t, --threads` | Number of threads |

When `--salmon-dir` points at a `quant` output directory, `deseq2` reuses `gene_counts.csv` or `gene_counts.tsv` if present. It falls back to importing from `quant.sf` only when no precomputed gene counts are available. Metadata sample IDs and count matrix sample IDs must match exactly.

### `rskit validate` / `rskit doctor`

Validate input files without running analysis tools.

| Option | Description |
|--------|-------------|
| `-S, --coldata` | Coldata file with a `sample` column |
| `--design` | DESeq2 design formula used to validate required metadata columns |
| `--check-reads` | Require `r1`/`r2` columns and verify read files exist |
| `-gc, --gene-counts` | Optional gene counts matrix to validate against coldata |
| `-e, --expression` | Optional expression matrix to validate against coldata |

### `rskit wgcna`

WGCNA co-expression network analysis.

| Option | Description |
|--------|-------------|
| `-e, --expression` | Expression matrix file |
| `-o, --output-dir` | Output directory |
| `-S, --coldata` | Sample metadata file |
| `-G, --gene-info` | Gene metadata file |
| `-sep, --sep` | Override separator for input files; CSV/TSV is auto-detected by default |
| `-n, --name` | Analysis name |
| `-s, --species` | Species for enrichment analysis |
| `-l, --level` | Data level: `gene` or `transcript` |
| `-nt, --network-type` | Network type |
| `-tom, --tom-type` | TOM type |
| `-min, --min-module-size` | Minimum module size |
| `-p, --power` | Soft thresholding power |
| `-rsquared, --rsquared-cut` | R-squared cutoff |
| `-mean, --mean-cut` | Mean connectivity cutoff |
| `-mediss, --mediss-thresh` | Module merging threshold |
| `-tpm, --tpm-cutoff` | TPM cutoff for filtering |

## Python API

### Quantification

```python
from rskit import RNAseqPipeline, PipelineConfig

config = PipelineConfig()
pipeline = RNAseqPipeline(config)

samples = {
    "sample1": {
        "fq1": "data/sample1_R1.fq",
        "fq2": "data/sample1_R2.fq",
    }
}

results = pipeline.run(
    samples=samples,
    genome_fasta="genome.fa",
    gtf_file="annotation.gtf",
    transcript_fasta="transcripts.fa",
    index_dir="STAR_index",
    output_dir="results/02_bam",
    quant_output_dir="results/03_quant",
)
```

### DESeq2

```python
from rskit.core.deseq2 import Deseq2Analyzer
from rskit.config import DESeq2Config

config = DESeq2Config(alpha=0.05, lfc_threshold=2.0)
analyzer = Deseq2Analyzer(config)

metadata_df = analyzer.load_metadata("coldata.csv", required_columns=["condition"])
counts_df = analyzer.load_counts_from_file("counts.csv", metadata_df=metadata_df)

results_df = analyzer.analyze(
    counts_df=counts_df,
    metadata_df=metadata_df,
    contrast=["condition", "treatment", "control"],
)

summary = analyzer.get_summary()
print(f"Significant genes: {summary['significant_genes']}")
```

### WGCNA

```python
from rskit.core.wgcna import WGCNAAnalyzer

analyzer = WGCNAAnalyzer(
    output_dir="./wgcna_results",
    name="MyWGCNA",
    network_type="signed hybrid",
    min_module_size=50,
)

analyzer.load_data(
    expression_file="expression.csv",
    coldata="coldata.csv",
    gene_info_file="gene_info.csv",
)

wgcna_obj = analyzer.run_analysis()
analyzer.save_results()
```

## Output Structure

### Quant output

```text
03_quant/
├── <sample>/quant.sf
├── gene_counts.csv
├── gene_tpm.csv
├── gene_log2_tpm_plus1.csv
└── tx2gene.tsv
```

### Complete pipeline output

```text
results/
├── 00_index/
├── 01_clean_data/
├── 02_bam/
├── 03_quant/
└── 04_deseq2/
```

### DESeq2 output

```text
04_deseq2/
├── deseq2_results.csv
├── gene_counts.csv
├── pca_plot.pdf
├── volcano_plot.pdf
└── ma_plot.pdf
```

### WGCNA output

```text
wgcna_results/
├── figures/
├── WGCNA.p
└── module_info.csv
```

## File Format Notes

- CSV/TSV separators are detected from file extension
- Gene expression matrices are stored as rows = samples, columns = genes
- DESeq2 count matrices may be provided as genes x samples only when coldata sample IDs identify the count columns
- DESeq2 metadata must contain `sample` plus every column referenced by `--design`
- Matrix sample IDs must match coldata sample IDs exactly
- Gene metadata should use the first column as the gene identifier
