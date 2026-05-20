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
- Example input templates with `rskit template`

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

## User Scenarios

### I need a valid input file to start from

Generate a small `coldata` template instead of guessing column names.

```bash
rskit template coldata -o coldata.csv
rskit template contrast -o contrast.tsv
```

### I want to catch input problems before a long run

Validate metadata columns, read paths, and count or expression matrix sample IDs without running STAR, Salmon, DESeq2, or WGCNA.

```bash
rskit validate -S coldata.csv -r -gc counts.csv -d "~batch + condition"
rskit doctor -S coldata.csv -e expression.csv
```

### I have FASTQ files and want the full RNA-seq workflow

Run alignment, Salmon quantification, gene-level export, and DESeq2. `quant` creates `tx2gene.tsv`, `gene_counts.csv`, `gene_tpm.csv`, and `gene_log2_tpm_plus1.csv` before DESeq2 needs gene-level counts.

```bash
rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/
rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ -d "~batch + condition"
```

### I only need quantification

Use single-sample mode for one library or `--coldata` for a batch.

```bash
rskit quant -s sample1 -1 sample1_R1.fq.gz -2 sample1_R2.fq.gz -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/
rskit quant -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/
```

### I already have gene counts and want DESeq2

Use a matrix directly, or point at a `quant` output directory and let rskit reuse `gene_counts.csv` or `gene_counts.tsv` when present.

```bash
rskit deseq2 -gc counts.csv -S coldata.csv -c condition,treatment,control
rskit deseq2 -sd ./03_quant -S coldata.csv -gtf annotation.gtf -d "~batch + condition"
```

### I want co-expression modules

Provide an expression matrix as rows = samples and columns = genes. Add coldata and gene metadata when available.

```bash
rskit wgcna -e expression.csv -o ./wgcna_results
rskit wgcna -e expression.csv -S coldata.csv -G gene_info.csv -o ./wgcna_results
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

| Short | Long | Description |
|-------|------|-------------|
| `-S` | `--coldata` | Required coldata file with `sample,id,condition,r1,r2`. Relative `r1`/`r2` paths are resolved from this file. |
| `-g` | `--genome-fasta` | Required genome FASTA used to build or check the STAR index. |
| `-gtf` | `--gtf-file` | Required GTF/GFF annotation used by STAR and to create `tx2gene.tsv` when `--tx2gene` is not provided. |
| `-gf` | `--transcript-fasta` | Required transcript FASTA used by Salmon quantification. |
| `-o` | `--output-dir` | Required workflow output directory; rskit creates `00_index/`, `02_bam/`, `03_quant/`, and `04_deseq2/` under it. |
| `-idx` | `--index-dir` | Optional existing STAR index directory; defaults to `<output-dir>/00_index`. |
| `-t2g` | `--tx2gene` | Optional transcript-to-gene mapping file; if omitted, rskit writes `03_quant/tx2gene.tsv` from `--gtf-file`. |
| `-t` | `--threads` | Threads per sample when running sequentially; also used for index and DESeq2 defaults. |
| `-p` | `--parallel` | Total cores for multi-sample parallel execution; rskit divides this across samples. |
| `-tr` | `--trim` | Run fastp before alignment and use trimmed FASTQ files. |
| `-fi` | `--force-index` | Rebuild the STAR index even when an index directory already exists. |
| `-se` | `--skip-existing` | Skip sample-level work when expected output files already exist. |
| `-d` | `--design` | DESeq2 design formula; every referenced column must exist in coldata. Default: `~condition`. |
| `-c` | `--contrast` | DESeq2 contrast as `factor,level1,level2`; factor and levels are validated against coldata. |
| `-a` | `--alpha` | Adjusted p-value threshold used for significance summaries. Default: `0.05`. |
| `-l` | `--lfc` | Absolute log2 fold-change threshold used for significance summaries. Default: `2.0`. |

### `rskit quant`

Complete quantification pipeline: index -> align -> quant -> gene-level table export.

| Short | Long | Description |
|-------|------|-------------|
| `-s` | `--sample` | Sample name for single-sample mode; use with `-1` and `-2`. |
| `-S` | `--coldata` | Batch sample table with `sample,r1,r2`; replaces `--sample`, `--r1`, and `--r2`. |
| `-1` | `--r1` | First read file for single-sample mode. |
| `-2` | `--r2` | Second read file for single-sample mode. |
| `-g` | `--genome-fasta` | Required genome FASTA used to build or check the STAR index. |
| `-gtf` | `--gtf-file` | Required annotation used by STAR and `tx2gene.tsv` generation. |
| `-gf` | `--transcript-fasta` | Required transcript FASTA used by Salmon. |
| `-o` | `--output-dir` | Required output/work directory. |
| `-idx` | `--index-dir` | Optional existing STAR index directory; defaults to `<output-dir>/00_index`. |
| `-t2g` | `--tx2gene` | Optional transcript-to-gene mapping for gene-level export; otherwise generated from `--gtf-file`. |
| `-t` | `--threads` | Threads per sample. Default: `8`. |
| `-p` | `--parallel` | Total cores for multi-sample parallel execution. |
| `-tr` | `--trim` | Run fastp before alignment. |
| `-fi` | `--force-index` | Rebuild the STAR index even if it exists. |
| `-se` | `--skip-existing` | Skip sample work when expected output already exists. |

### `rskit deseq2`

DESeq2 differential expression analysis.

| Short | Long | Description |
|-------|------|-------------|
| `-sd` | `--salmon-dir` | Directory containing Salmon `quant.sf` folders or precomputed `gene_counts.csv`/`gene_counts.tsv`. Mutually exclusive with `--gene-counts`. |
| `-gc` | `--gene-counts` | Gene counts matrix file. Mutually exclusive with `--salmon-dir`. |
| `-S` | `--coldata` | Required metadata file with `sample` and every column referenced by `--design`. |
| `-gtf` | `--gtf` | GTF/GFF annotation; required when importing `quant.sf` without `--tx2gene`. |
| `-t2g` | `--tx2gene` | Optional transcript-to-gene mapping used when importing from Salmon outputs. |
| `-w` | `--work-dir` | Work directory used to place the default `04_deseq2/` output directory. Default: current directory. |
| `-o` | `--output-dir` | Custom DESeq2 output directory; overrides `<work-dir>/04_deseq2`. |
| `-d` | `--design` | DESeq2 design formula. Default: `~condition`. |
| `-c` | `--contrast` | Contrast as `factor,level1,level2`; validated against coldata before counts are loaded. |
| `-a` | `--alpha` | Adjusted p-value threshold used for result summaries. Default: `0.05`. |
| `-l` | `--lfc` | Absolute log2 fold-change threshold used for result summaries. Default: `2.0`. |
| `-t` | `--threads` | Number of CPUs for PyDESeq2 inference. |

When `--salmon-dir` points at a `quant` output directory, `deseq2` reuses `gene_counts.csv` or `gene_counts.tsv` if present. It falls back to importing from `quant.sf` only when no precomputed gene counts are available. Metadata sample IDs and count matrix sample IDs must match exactly.

`--contrast` must use `factor,level1,level2`. The factor must be a coldata column, and both levels must exist in that column; rskit validates this before loading counts or running DESeq2.

### `rskit validate` / `rskit doctor`

Validate input files without running analysis tools.

| Short | Long | Description |
|-------|------|-------------|
| `-S` | `--coldata` | Required coldata file with a `sample` column. |
| `-d` | `--design` | DESeq2 design formula used to check required metadata columns. Default: `~condition`. |
| `-r` | `--check-reads` | Require `r1`/`r2` columns and verify that read files exist. |
| `-gc` | `--gene-counts` | Optional gene counts matrix to validate against coldata sample IDs. |
| `-e` | `--expression` | Optional expression matrix to validate against coldata sample IDs. |

### `rskit template`

Write example input template files.

| Short | Long | Description |
|-------|------|-------------|
| n/a | `template_name` | Positional template type: `coldata` or `contrast`. |
| `-o` | `--output` | Required output path; `.csv` writes CSV and `.tsv`/`.txt` writes TSV. |
| `-f` | `--force` | Overwrite the output file if it already exists. |

### `rskit wgcna`

WGCNA co-expression network analysis.

| Short | Long | Description |
|-------|------|-------------|
| `-e` | `--expression` | Required expression matrix; rows must be samples and columns genes. |
| `-o` | `--output-dir` | Required WGCNA output directory. |
| `-S` | `--coldata` | Optional sample metadata file; sample IDs must match expression rows. |
| `-G` | `--gene-info` | Optional gene metadata file indexed by gene ID. |
| `-sp` (`-sep`) | `--sep` | Optional separator override; by default `.csv` uses comma and `.tsv`/`.txt` use tab. |
| `-n` | `--name` | Analysis name stored in the PyWGCNA object. Default: `WGCNA`. |
| `-s` | `--species` | Species label used by PyWGCNA enrichment analysis. |
| `-l` | `--level` | Data level passed to PyWGCNA: `gene` or `transcript`. Default: `gene`. |
| `-nw` (`-nt`) | `--network-type` | WGCNA network type: `unsigned`, `signed`, or `signed hybrid`. |
| `-tt` (`-tom`) | `--tom-type` | TOM type passed to PyWGCNA: `unsigned` or `signed`. |
| `-ms` (`-min`) | `--min-module-size` | Minimum module size for module detection. Default: `50`. |
| `-p` | `--power` | Soft-thresholding power; omit to let PyWGCNA auto-detect. |
| `-rs` (`-rsquared`) | `--rsquared-cut` | R-squared cutoff for power selection. Default: `0.9`. |
| `-mc` (`-mean`) | `--mean-cut` | Mean connectivity cutoff. Default: `100`. |
| `-md` (`-mediss`) | `--mediss-thresh` | Module eigengene dissimilarity threshold for merging modules. Default: `0.2`. |
| `-tc` (`-tpm`) | `--tpm-cutoff` | TPM cutoff used by PyWGCNA filtering. Default: `1`. |

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
├── manifest.json
├── pca_plot.pdf
├── volcano_plot.pdf
└── ma_plot.pdf
```

`manifest.json` records the DESeq2 inputs, resolved counts file, sample IDs, design, contrast, summary, and key output files.

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
