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

| Dependency | Version |
|------------|---------|
| Python | >= 3.8 |
| pandas | 3.0.3 |
| numpy | 2.5.0 |
| pydeseq2 | 0.5.4 |
| pytximport | 0.13.0 |
| PyWGCNA | 2.2.1 |
| STAR | 2.7.11b |
| Salmon | 2.2.1 |
| fastp | 1.3.6 (optional) |

## User Scenarios

### I need a valid input file to start from

Generate a small `coldata` template instead of guessing column names.

```bash
rskit template coldata -o coldata.csv
```

### I want to catch input problems before a long run

Validate metadata columns, read paths, and count or expression matrix sample IDs without running STAR, Salmon, DESeq2, or WGCNA.

```bash
rskit validate -S coldata.csv -r -gc counts.csv -d "~batch + condition"
rskit doctor -S coldata.csv -e expression.csv
```

### I have FASTQ files and want the full RNA-seq workflow

Run alignment, Salmon quantification, gene-level export, and DESeq2. `quant` creates `tx2gene.tsv`, `gene_counts.csv`, `gene_tpm.csv`, and `gene_log2_tpm.csv` before DESeq2 needs gene-level counts.

```bash
rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/
rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ -t 100 -j 20
rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ -d "~batch + condition"
```

### I only need quantification

Use single-sample mode for one library or `--coldata` for a batch.

```bash
rskit quant -s sample1 -1 sample1_R1.fq.gz -2 sample1_R2.fq.gz -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/
rskit quant -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/
rskit quant -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ -t 100 -j 20
rskit quant -s sample1 -1 sample1_R1.fq.gz -2 sample1_R2.fq.gz -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ -ms
```

Default batch quantification uses `-j 1`. `-t/--threads` is the total thread budget, and `-j/--jobs` is the sample concurrency. For example, `-t 100 -j 20` gives each sample 5 threads. Use `-ms/--merge-sf` to regenerate gene-level CSV files from all `03_quant/*/quant.sf` folders.

### I need to tune STAR, Salmon, or fastp

Pass advanced tool arguments as a quoted string. If an allowed argument conflicts with an rskit default, the user-provided value replaces the default. rskit rejects arguments that would change managed inputs, outputs, reports, indexes, library type, or thread counts.

```bash
rskit quant -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ --star-args "--outFilterMultimapNmax 8"
rskit quant -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ --salmon-args "--validateMappings --minScoreFraction 0.95"
rskit quant -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ -tr --fastp-args "--length_required 30 --cut_front"
```

### I already have gene counts and want DESeq2

Use a matrix directly, or point at a `quant` output directory and let rskit reuse `gene_counts.csv` or `gene_counts.tsv` when present.

```bash
rskit deseq2 -gc counts.csv -S coldata.csv -c condition,treatment,control
rskit deseq2 -sd ./03_quant -S coldata.csv -gtf annotation.gtf -d "~batch + condition"
```

### I want co-expression modules

Provide an expression matrix as rows = genes and columns = samples. Add coldata and gene metadata when available.

```bash
rskit wgcna -e expression.csv -o ./wgcna_results
rskit wgcna -e expression.csv -S coldata.csv -G gene_info.csv -o ./wgcna_results
```

## Format

### `coldata.csv`

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

`r1` and `r2` are FASTQ paths and can be relative to `coldata.csv`.

### `counts.csv`

Used by `rskit deseq2 -gc`.

```csv
gene_id,sample1,sample2,sample3,sample4
geneA,10,12,80,77
geneB,0,1,4,5
```

- rows: genes
- columns: samples
- first column: gene ID

### `expression.csv`

Used by `rskit wgcna -e`. Values are **TPM** (or any normalized expression measure), not raw counts — `wgcna` does not require integers. Use `gene_tpm.csv` from `rskit quant`/`rskit all`, or a `gene_log2_tpm.csv`/log-transformed matrix; raw `gene_counts.csv` is not recommended for WGCNA.

```csv
gene_id,sample1,sample2,sample3,sample4
geneA,3.46,3.70,6.34,6.29
geneB,0.00,1.00,2.32,2.58
```

- rows: genes
- columns: samples
- first column: gene ID
- values: normalized expression (TPM recommended)

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
| `-t` | `--threads` | Total thread budget for sample processing. Default: `8`. |
| `-j` | `--jobs` | Maximum number of samples to process concurrently. Default: `1`. |
| `-ms` | `--merge-sf` | Scan all `03_quant/*/quant.sf` files and regenerate gene-level CSV files. |
| `-tr` | `--trim` | Run fastp before alignment and use trimmed FASTQ files. |
| `-fi` | `--force-index` | Rebuild the STAR index even when an index directory already exists. |
| `-se` | `--skip-existing` | Skip sample-level work when expected output files already exist. |
| n/a | `--star-args` | Advanced STAR arguments. Allowed conflicts replace rskit defaults; protected options include `--runThreadN`, `--genomeDir`, `--readFilesIn`, `--readFilesCommand`, `--outFileNamePrefix`, `--outSAMtype`, `--quantMode`, `--genomeFastaFiles`, and `--sjdbGTFfile`. |
| n/a | `--salmon-args` | Advanced `salmon quant` arguments. Allowed conflicts replace rskit defaults; protected options include `-t`/`--targets`, `-a`/`--alignments`, `-o`/`--output`, `-p`/`--threads`, and `-l`/`--libType`. |
| n/a | `--fastp-args` | Advanced fastp arguments used only with `--trim`. Allowed conflicts replace rskit defaults; protected options include `-i`/`--in1`, `-I`/`--in2`, `-o`/`--out1`, `-O`/`--out2`, `-w`/`--thread`, report paths, STDIN/STDOUT, and extra output-file options. |
| `-d` | `--design` | DESeq2 design formula; every referenced column must exist in coldata. Default: `~condition`. |
| `-c` | `--contrast` | DESeq2 contrast as `factor,level1,level2`; factor and levels are validated against coldata. |
| `-a` | `--alpha` | Adjusted p-value threshold used for significance summaries. Default: `0.05`. |
| `-l` | `--lfc` | Absolute log2 fold-change threshold used for significance summaries. Default: `2.0`. |
| `-F` | `--min-count` | Minimum total count for DESeq2 gene prefiltering. Default: `10`; use `0` to disable. |

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
| `-t` | `--threads` | Total thread budget for sample processing. Default: `8`. |
| `-j` | `--jobs` | Maximum number of samples to process concurrently. Default: `1`. |
| `-ms` | `--merge-sf` | Scan all `03_quant/*/quant.sf` files and regenerate gene-level CSV files. |
| `-tr` | `--trim` | Run fastp before alignment. |
| `-fi` | `--force-index` | Rebuild the STAR index even if it exists. |
| `-se` | `--skip-existing` | Skip sample work when expected output already exists. |
| n/a | `--star-args` | Advanced STAR arguments. Allowed conflicts replace rskit defaults; protected options include `--runThreadN`, `--genomeDir`, `--readFilesIn`, `--readFilesCommand`, `--outFileNamePrefix`, `--outSAMtype`, `--quantMode`, `--genomeFastaFiles`, and `--sjdbGTFfile`. |
| n/a | `--salmon-args` | Advanced `salmon quant` arguments. Allowed conflicts replace rskit defaults; protected options include `-t`/`--targets`, `-a`/`--alignments`, `-o`/`--output`, `-p`/`--threads`, and `-l`/`--libType`. |
| n/a | `--fastp-args` | Advanced fastp arguments used only with `--trim`. Allowed conflicts replace rskit defaults; protected options include `-i`/`--in1`, `-I`/`--in2`, `-o`/`--out1`, `-O`/`--out2`, `-w`/`--thread`, report paths, STDIN/STDOUT, and extra output-file options. |

### `rskit deseq2`

DESeq2 differential expression analysis.

| Short | Long | Description |
|-------|------|-------------|
| `-sd` | `--salmon-dir` | Directory containing Salmon `quant.sf` folders or precomputed `gene_counts.csv`/`gene_counts.tsv`. Mutually exclusive with `--gene-counts`. |
| `-gc` | `--gene-counts` | Gene counts matrix file with rows = genes and columns = samples. Mutually exclusive with `--salmon-dir`. |
| `-S` | `--coldata` | Required metadata file with `sample` and every column referenced by `--design`. |
| `-gtf` | `--gtf` | GTF/GFF annotation; required when importing `quant.sf` without `--tx2gene`. |
| `-t2g` | `--tx2gene` | Optional transcript-to-gene mapping used when importing from Salmon outputs. |
| `-w` | `--work-dir` | Work directory used to place the default `04_deseq2/` output directory. Default: current directory. |
| `-o` | `--output-dir` | Custom DESeq2 output directory; overrides `<work-dir>/04_deseq2`. |
| `-d` | `--design` | DESeq2 design formula. Default: `~condition`. |
| `-c` | `--contrast` | Contrast as `factor,level1,level2`; validated against coldata before counts are loaded. |
| `-a` | `--alpha` | Adjusted p-value threshold used for result summaries. Default: `0.05`. |
| `-l` | `--lfc` | Absolute log2 fold-change threshold used for result summaries. Default: `2.0`. |
| `-F` | `--min-count` | Minimum total count for DESeq2 gene prefiltering. Default: `10`; use `0` to disable. |
| `-t` | `--threads` | Number of CPUs for PyDESeq2 inference. |

### `rskit validate` / `rskit doctor`

Validate input files without running analysis tools.

| Short | Long | Description |
|-------|------|-------------|
| `-S` | `--coldata` | Required coldata file with a `sample` column. |
| `-d` | `--design` | DESeq2 design formula used to check required metadata columns. Default: `~condition`. |
| `-r` | `--check-reads` | Require `r1`/`r2` columns and verify that read files exist. |
| `-gc` | `--gene-counts` | Optional genes x samples count matrix to validate against coldata sample IDs. |
| `-e` | `--expression` | Optional genes x samples expression matrix to validate against coldata sample IDs. |

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
| `-e` | `--expression` | Required expression matrix; rows must be genes and columns samples. |
| `-o` | `--output-dir` | Required WGCNA output directory. |
| `-S` | `--coldata` | Optional sample metadata file; sample IDs must match expression columns. |
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

config = DESeq2Config(alpha=0.05, lfc_threshold=2.0, prefilter_min_count=10)
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
├── gene_log2_tpm.csv
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
