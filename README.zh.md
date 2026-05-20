# rskit - RNA-seq 分析工具包

<!-- README-I18N:START -->

[English](./README.md) | **中文**

<!-- README-I18N:END -->

rskit 是一个用于 RNA-seq 分析的 Python 工具包，提供 CLI 和 Python API，覆盖 read alignment、Salmon quantification、DESeq2 differential expression 和 WGCNA 等常见流程。

## 功能

- 定量流程：STAR alignment + Salmon quantification
- 在 `quant` 阶段导出基因级表达矩阵
- 使用 DESeq2 进行差异表达分析
- 使用 WGCNA 进行共表达网络分析
- 使用 `rskit all` 运行端到端工作流
- 自动识别 CSV/TSV
- 子命令共享统一的 `--coldata` 元数据格式
- 在 DESeq2 之前提前导出 `tx2gene.tsv` 和基因级矩阵
- 严格校验 metadata 与 count/expression matrix 的 sample ID
- 使用 `rskit validate` / `rskit doctor` 进行输入预检

## 安装

### 从源码安装
```bash
git clone https://github.com/WWz33/rskit.git
cd rskit
pip install -e .
```

### 依赖

- Python >= 3.8
- pandas
- numpy
- pydeseq2
- pytximport
- PyWGCNA
- STAR
- Salmon
- fastp (optional)

## 快速开始

### 1. 完整流程

```bash
rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/

# With custom design formula
rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ --design "~batch + condition"
```

### 2. 定量流程

`quant` 会在 `03_quant/` 中写入每个样本的 Salmon 输出，以及 `tx2gene.tsv`、`gene_counts.csv`、`gene_tpm.csv` 和 `gene_log2_tpm_plus1.csv`。transcript-to-gene 映射在定量阶段生成，早于 DESeq2 对基因级 counts 的使用。

```bash
# Single sample
rskit quant -s sample1 -1 sample1_R1.fq.gz -2 sample1_R2.fq.gz -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/

# Multiple samples using coldata
rskit quant -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/
```

### 3. 差异表达分析

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

### 5. 输入预检

```bash
# 在运行分析前检查 metadata 列、read 路径和 counts/sample 对齐
rskit validate -S coldata.csv --check-reads -gc counts.csv --design "~batch + condition"

# doctor 是 validate 的别名
rskit doctor -S coldata.csv -e expression.csv
```

## Coldata 格式

所有子命令都使用同一个 `--coldata` / `-S` 参数。

```csv
sample,id,condition,r1,r2
sample1,ctrl,control,sample1_R1.fq.gz,sample1_R2.fq.gz
sample2,ctrl,control,sample2_R1.fq.gz,sample2_R2.fq.gz
sample3,treat,treatment,sample3_R1.fq.gz,sample3_R2.fq.gz
sample4,treat,treatment,sample4_R1.fq.gz,sample4_R2.fq.gz
```

每个子命令只读取自己需要的列：

- `quant`: `sample`, `r1`, `r2`
- `deseq2`: `sample` 加上 `--design` 引用的所有 metadata 列；默认 `~condition` 需要 `condition`
- `wgcna`: `sample` 加任意 metadata 列

相对路径形式的 `r1` 和 `r2` 会按 coldata 文件所在目录解析。Count 和 expression matrix 期望使用 rows = samples、columns = genes。如果 DESeq2 counts 文件是 genes x samples，只有当 coldata 的 sample ID 与 counts 的列名匹配时才会转置。定向之后，sample ID 必须与 coldata 完全一致；rskit 不会静默丢弃或取交集。

## 命令参考

### `rskit all`

完整流程：定量 + DESeq2 分析。

| Option | Description |
|--------|-------------|
| `-S, --coldata` | 包含 `sample,id,condition,r1,r2` 的 coldata 文件 |
| `-g, --genome-fasta` | Genome FASTA 文件 |
| `-gtf, --gtf-file` | GTF annotation 文件 |
| `-gf, --transcript-fasta` | Transcript FASTA 文件 |
| `-o, --output-dir` | 输出目录 |
| `-idx, --index-dir` | STAR index 目录 |
| `-t2g, --tx2gene` | 已有 transcript-to-gene mapping 文件；否则在 quant 阶段从 `--gtf-file` 创建 `tx2gene.tsv` |
| `-t, --threads` | 每个样本使用的线程数 |
| `-p, --parallel` | 并行执行使用的总核心数 |
| `--trim` | 使用 fastp 修剪 reads |
| `--force-index` | 强制重建 STAR index |
| `--skip-existing` | 已有输出时跳过样本任务 |
| `--design` | DESeq2 design formula |
| `--contrast` | Contrast specification |
| `--alpha` | Adjusted p-value threshold |
| `--lfc` | Log2 fold-change threshold |

### `rskit quant`

完整定量流程：index -> align -> quant -> 基因级表格导出。

| Option | Description |
|--------|-------------|
| `-s, --sample` | 单样本模式的样本名 |
| `-S, --coldata` | 包含 `sample,r1,r2` 的样本文件 |
| `-1, --r1` | First read 文件 |
| `-2, --r2` | Second read 文件 |
| `-g, --genome-fasta` | Genome FASTA 文件 |
| `-gtf, --gtf-file` | GTF annotation 文件 |
| `-gf, --transcript-fasta` | Transcript FASTA 文件 |
| `-o, --output-dir` | 输出目录 |
| `-idx, --index-dir` | STAR index 目录 |
| `-t2g, --tx2gene` | 用于基因级导出的已有 transcript-to-gene mapping 文件；否则从 `--gtf-file` 创建 `tx2gene.tsv` |
| `-t, --threads` | 每个样本使用的线程数 |
| `-p, --parallel` | 并行执行使用的总核心数 |
| `--trim` | 使用 fastp 修剪 reads |
| `--force-index` | 强制重建 STAR index |
| `--skip-existing` | 已有输出时跳过样本任务 |

### `rskit deseq2`

DESeq2 差异表达分析。

| Option | Description |
|--------|-------------|
| `-sd, --salmon-dir` | 包含 Salmon quant 文件夹的目录 |
| `-gc, --gene-counts` | Gene counts matrix 文件 |
| `-S, --coldata` | 包含 `sample` 和 `--design` 引用列的 sample metadata 文件 |
| `-gtf, --gtf` | GTF annotation 文件 |
| `-t2g, --tx2gene` | Transcript-to-gene mapping 文件 |
| `--design` | Design formula |
| `--contrast` | Contrast specification |
| `--alpha` | Adjusted p-value threshold |
| `--lfc` | Log2 fold-change threshold |
| `-o, --output-dir` | 输出目录 |
| `-t, --threads` | 线程数 |

当 `--salmon-dir` 指向 `quant` 输出目录时，`deseq2` 会优先复用已有的 `gene_counts.csv` 或 `gene_counts.tsv`。只有不存在预计算 gene counts 时，才会从 `quant.sf` 重新导入。Metadata sample IDs 和 count matrix sample IDs 必须完全匹配。

### `rskit validate` / `rskit doctor`

不运行分析工具，只校验输入文件。

| Option | Description |
|--------|-------------|
| `-S, --coldata` | 包含 `sample` 列的 coldata 文件 |
| `--design` | 用于校验 metadata 必需列的 DESeq2 design formula |
| `--check-reads` | 要求存在 `r1`/`r2` 列，并检查 read 文件是否存在 |
| `-gc, --gene-counts` | 可选 gene counts matrix，用于和 coldata 校验 |
| `-e, --expression` | 可选 expression matrix，用于和 coldata 校验 |

### `rskit wgcna`

WGCNA 共表达网络分析。

| Option | Description |
|--------|-------------|
| `-e, --expression` | Expression matrix 文件 |
| `-o, --output-dir` | 输出目录 |
| `-S, --coldata` | Sample metadata 文件 |
| `-G, --gene-info` | Gene metadata 文件 |
| `-sep, --sep` | 输入文件分隔符 |
| `-n, --name` | 分析名称 |
| `-s, --species` | 用于 enrichment analysis 的物种 |
| `-l, --level` | 数据层级：`gene` 或 `transcript` |
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

## 输出结构

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

## 文件格式说明

- CSV/TSV 分隔符由文件扩展名自动识别
- Gene expression matrices 使用 rows = samples、columns = genes
- DESeq2 count matrices 只有在 coldata sample IDs 能识别 count columns 时，才可以用 genes x samples 形式提供
- DESeq2 metadata 必须包含 `sample` 和 `--design` 引用的所有列
- Matrix sample IDs 必须与 coldata sample IDs 完全一致
- Gene metadata 应使用第一列作为 gene identifier
