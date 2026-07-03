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
- 使用 `rskit template` 生成示例输入模板

## 安装

### 从源码安装
```bash
git clone https://github.com/WWz33/rskit.git
cd rskit
pip install -e .
```

### 依赖

| 依赖 | 版本 |
|------|------|
| Python | >= 3.8 |
| pandas | 3.0.3 |
| numpy | 2.5.0 |
| pydeseq2 | 0.5.4 |
| pytximport | 0.13.0 |
| PyWGCNA | 2.2.1 |
| STAR | 2.7.11b |
| Salmon | 2.2.1 |
| fastp | 1.3.6 (optional) |

## 用户场景

### 我需要一个可用的输入文件起点

先生成一个小的 `coldata` 模板，不需要猜字段名。

```bash
rskit template coldata -o coldata.csv
rskit template contrast -o contrast.tsv
```

### 我想在长时间运行前发现输入问题

只校验 metadata 列、read 路径、count 或 expression matrix 的 sample ID，不运行 STAR、Salmon、DESeq2 或 WGCNA。

```bash
rskit validate -S coldata.csv -r -gc counts.csv -d "~batch + condition"
rskit doctor -S coldata.csv -e expression.csv
```

### 我有 FASTQ 文件，想跑完整 RNA-seq 流程

运行 alignment、Salmon quantification、基因级导出和 DESeq2。`quant` 会在 DESeq2 使用基因级 counts 前创建 `tx2gene.tsv`、`gene_counts.csv`、`gene_tpm.csv` 和 `gene_log2_tpm.csv`。

```bash
rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/
rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ -t 100 -j 20
rskit all -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ -d "~batch + condition"
```

### 我只需要定量

单个文库用 single-sample 模式，批量样本用 `--coldata`。

```bash
rskit quant -s sample1 -1 sample1_R1.fq.gz -2 sample1_R2.fq.gz -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/
rskit quant -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/
rskit quant -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ -t 100 -j 20
rskit quant -s sample1 -1 sample1_R1.fq.gz -2 sample1_R2.fq.gz -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ -ms
```

默认批量定量为 `-j 1`。`-t/--threads` 是总线程数，`-j/--jobs` 是并发样本数；例如 `-t 100 -j 20` 为每个样本分配 5 线程。使用 `-ms/--merge-sf` 时，会从所有 `03_quant/*/quant.sf` 重新生成基因级 CSV。

### 我需要调整 STAR、Salmon 或 fastp 参数

把底层软件的高级参数作为一个带引号的字符串传入。如果允许覆盖的参数与 rskit 默认值冲突，用户传入的值会替换默认值。rskit 会拒绝会改变输入、输出、报告、index、library type 或线程数的参数。

```bash
rskit quant -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ --star-args "--outFilterMultimapNmax 8"
rskit quant -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ --salmon-args "--validateMappings --minScoreFraction 0.95"
rskit quant -S coldata.csv -g genome.fa -gtf annotation.gtf -gf transcripts.fa -o results/ -tr --fastp-args "--length_required 30 --cut_front"
```

### 我已经有 gene counts，想直接做 DESeq2

可以直接提供矩阵，也可以指向 `quant` 输出目录，让 rskit 优先复用已有的 `gene_counts.csv` 或 `gene_counts.tsv`。

```bash
rskit deseq2 -gc counts.csv -S coldata.csv -c condition,treatment,control
rskit deseq2 -sd ./03_quant -S coldata.csv -gtf annotation.gtf -d "~batch + condition"
```

### 我想分析共表达模块

expression matrix 使用 rows = genes、columns = samples。有 coldata 和 gene metadata 时一起提供。

```bash
rskit wgcna -e expression.csv -o ./wgcna_results
rskit wgcna -e expression.csv -S coldata.csv -G gene_info.csv -o ./wgcna_results
```

## 格式说明

### `coldata.csv`

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

`r1` 和 `r2` 是 FASTQ 路径，可以相对 `coldata.csv`。

### `counts.csv`

用于 `rskit deseq2 -gc`。

```csv
gene_id,sample1,sample2,sample3,sample4
geneA,10,12,80,77
geneB,0,1,4,5
```

- rows: genes
- columns: samples
- 第一列：gene ID

### `expression.csv`

用于 `rskit wgcna -e`。数值为 **TPM**（或任意归一化表达量），不是原始 counts —— `wgcna` 不要求整数。可直接用 `rskit quant`/`rskit all` 产出的 `gene_tpm.csv`，或 `gene_log2_tpm.csv`/对数转换矩阵；不建议用原始 `gene_counts.csv` 跑 WGCNA。

```csv
gene_id,sample1,sample2,sample3,sample4
geneA,3.46,3.70,6.34,6.29
geneB,0.00,1.00,2.32,2.58
```

- rows: genes
- columns: samples
- 第一列：gene ID
- 数值：归一化表达量（推荐 TPM）

## 命令参考

### `rskit all`

完整流程：定量 + DESeq2 分析。

| 简写 | 长参数 | 说明 |
|------|--------|------|
| `-S` | `--coldata` | 必需的 coldata 文件，包含 `sample,id,condition,r1,r2`；相对 `r1`/`r2` 路径按该文件所在目录解析。 |
| `-g` | `--genome-fasta` | 必需的 genome FASTA，用于构建或检查 STAR index。 |
| `-gtf` | `--gtf-file` | 必需的 GTF/GFF annotation，用于 STAR，并在未提供 `--tx2gene` 时生成 `tx2gene.tsv`。 |
| `-gf` | `--transcript-fasta` | 必需的 transcript FASTA，用于 Salmon 定量。 |
| `-o` | `--output-dir` | 必需的流程输出目录；rskit 会在其中创建 `00_index/`、`02_bam/`、`03_quant/` 和 `04_deseq2/`。 |
| `-idx` | `--index-dir` | 可选的已有 STAR index 目录；默认使用 `<output-dir>/00_index`。 |
| `-t2g` | `--tx2gene` | 可选 transcript-to-gene mapping；未提供时从 `--gtf-file` 写出 `03_quant/tx2gene.tsv`。 |
| `-t` | `--threads` | 样本处理使用的总线程预算。默认：`8`。 |
| `-j` | `--jobs` | 同时处理的最大样本数。默认：`1`。 |
| `-ms` | `--merge-sf` | 扫描所有 `03_quant/*/quant.sf` 并重新生成基因级 CSV。 |
| `-tr` | `--trim` | alignment 前运行 fastp，并使用修剪后的 FASTQ。 |
| `-fi` | `--force-index` | 即使 index 目录已存在，也强制重建 STAR index。 |
| `-se` | `--skip-existing` | 目标输出已存在时跳过对应样本任务。 |
| n/a | `--star-args` | STAR 高级参数。允许覆盖的冲突参数会替换 rskit 默认值；受保护参数包括 `--runThreadN`、`--genomeDir`、`--readFilesIn`、`--readFilesCommand`、`--outFileNamePrefix`、`--outSAMtype`、`--quantMode`、`--genomeFastaFiles` 和 `--sjdbGTFfile`。 |
| n/a | `--salmon-args` | `salmon quant` 高级参数。允许覆盖的冲突参数会替换 rskit 默认值；受保护参数包括 `-t`/`--targets`、`-a`/`--alignments`、`-o`/`--output`、`-p`/`--threads` 和 `-l`/`--libType`。 |
| n/a | `--fastp-args` | fastp 高级参数，只在使用 `--trim` 时生效。允许覆盖的冲突参数会替换 rskit 默认值；受保护参数包括 `-i`/`--in1`、`-I`/`--in2`、`-o`/`--out1`、`-O`/`--out2`、`-w`/`--thread`、报告路径、STDIN/STDOUT 和额外输出文件参数。 |
| `-d` | `--design` | DESeq2 design formula；引用的每一列都必须存在于 coldata。默认：`~condition`。 |
| `-c` | `--contrast` | DESeq2 contrast，格式为 `factor,level1,level2`；factor 和 level 会按 coldata 校验。 |
| `-a` | `--alpha` | summary 使用的 adjusted p-value 阈值。默认：`0.05`。 |
| `-l` | `--lfc` | summary 使用的绝对 log2 fold-change 阈值。默认：`2.0`。 |
| `-F` | `--min-count` | DESeq2 预过滤的 gene total count 下限。默认：`10`；设为 `0` 可关闭。 |

### `rskit quant`

完整定量流程：index -> align -> quant -> 基因级表格导出。

| 简写 | 长参数 | 说明 |
|------|--------|------|
| `-s` | `--sample` | 单样本模式的样本名；与 `-1` 和 `-2` 一起使用。 |
| `-S` | `--coldata` | 批量样本表，包含 `sample,r1,r2`；用于替代 `--sample`、`--r1` 和 `--r2`。 |
| `-1` | `--r1` | 单样本模式的 first read 文件。 |
| `-2` | `--r2` | 单样本模式的 second read 文件。 |
| `-g` | `--genome-fasta` | 必需的 genome FASTA，用于构建或检查 STAR index。 |
| `-gtf` | `--gtf-file` | 必需 annotation，用于 STAR 和 `tx2gene.tsv` 生成。 |
| `-gf` | `--transcript-fasta` | 必需 transcript FASTA，用于 Salmon。 |
| `-o` | `--output-dir` | 必需输出/工作目录。 |
| `-idx` | `--index-dir` | 可选已有 STAR index 目录；默认使用 `<output-dir>/00_index`。 |
| `-t2g` | `--tx2gene` | 可选 transcript-to-gene mapping，用于基因级导出；未提供时从 `--gtf-file` 生成。 |
| `-t` | `--threads` | 样本处理使用的总线程预算。默认：`8`。 |
| `-j` | `--jobs` | 同时处理的最大样本数。默认：`1`。 |
| `-ms` | `--merge-sf` | 扫描所有 `03_quant/*/quant.sf` 并重新生成基因级 CSV。 |
| `-tr` | `--trim` | alignment 前运行 fastp。 |
| `-fi` | `--force-index` | 即使 STAR index 已存在也强制重建。 |
| `-se` | `--skip-existing` | 目标输出已存在时跳过样本任务。 |
| n/a | `--star-args` | STAR 高级参数。允许覆盖的冲突参数会替换 rskit 默认值；受保护参数包括 `--runThreadN`、`--genomeDir`、`--readFilesIn`、`--readFilesCommand`、`--outFileNamePrefix`、`--outSAMtype`、`--quantMode`、`--genomeFastaFiles` 和 `--sjdbGTFfile`。 |
| n/a | `--salmon-args` | `salmon quant` 高级参数。允许覆盖的冲突参数会替换 rskit 默认值；受保护参数包括 `-t`/`--targets`、`-a`/`--alignments`、`-o`/`--output`、`-p`/`--threads` 和 `-l`/`--libType`。 |
| n/a | `--fastp-args` | fastp 高级参数，只在使用 `--trim` 时生效。允许覆盖的冲突参数会替换 rskit 默认值；受保护参数包括 `-i`/`--in1`、`-I`/`--in2`、`-o`/`--out1`、`-O`/`--out2`、`-w`/`--thread`、报告路径、STDIN/STDOUT 和额外输出文件参数。 |

### `rskit deseq2`

DESeq2 差异表达分析。

| 简写 | 长参数 | 说明 |
|------|--------|------|
| `-sd` | `--salmon-dir` | 包含 Salmon `quant.sf` 文件夹或预计算 `gene_counts.csv`/`gene_counts.tsv` 的目录；与 `--gene-counts` 互斥。 |
| `-gc` | `--gene-counts` | Gene counts matrix 文件，rows = genes、columns = samples；与 `--salmon-dir` 互斥。 |
| `-S` | `--coldata` | 必需 metadata 文件，包含 `sample` 和 `--design` 引用的所有列。 |
| `-gtf` | `--gtf` | GTF/GFF annotation；从 `quant.sf` 导入且没有 `--tx2gene` 时需要。 |
| `-t2g` | `--tx2gene` | 从 Salmon 输出导入时使用的可选 transcript-to-gene mapping。 |
| `-w` | `--work-dir` | 用于放置默认 `04_deseq2/` 输出目录的工作目录。默认：当前目录。 |
| `-o` | `--output-dir` | 自定义 DESeq2 输出目录；覆盖 `<work-dir>/04_deseq2`。 |
| `-d` | `--design` | DESeq2 design formula。默认：`~condition`。 |
| `-c` | `--contrast` | Contrast，格式为 `factor,level1,level2`；加载 counts 前按 coldata 校验。 |
| `-a` | `--alpha` | summary 使用的 adjusted p-value 阈值。默认：`0.05`。 |
| `-l` | `--lfc` | summary 使用的绝对 log2 fold-change 阈值。默认：`2.0`。 |
| `-F` | `--min-count` | DESeq2 预过滤的 gene total count 下限。默认：`10`；设为 `0` 可关闭。 |
| `-t` | `--threads` | PyDESeq2 inference 使用的 CPU 数。 |

### `rskit validate` / `rskit doctor`

不运行分析工具，只校验输入文件。

| 简写 | 长参数 | 说明 |
|------|--------|------|
| `-S` | `--coldata` | 必需 coldata 文件，包含 `sample` 列。 |
| `-d` | `--design` | 用于检查 metadata 必需列的 DESeq2 design formula。默认：`~condition`。 |
| `-r` | `--check-reads` | 要求存在 `r1`/`r2` 列，并检查 read 文件是否存在。 |
| `-gc` | `--gene-counts` | 可选 genes x samples count matrix，用于按 coldata sample ID 校验。 |
| `-e` | `--expression` | 可选 genes x samples expression matrix，用于按 coldata sample ID 校验。 |

### `rskit template`

写出示例输入模板文件。

| 简写 | 长参数 | 说明 |
|------|--------|------|
| n/a | `template_name` | 位置参数，模板类型：`coldata` 或 `contrast`。 |
| `-o` | `--output` | 必需输出路径；`.csv` 写 CSV，`.tsv`/`.txt` 写 TSV。 |
| `-f` | `--force` | 输出文件已存在时允许覆盖。 |

### `rskit wgcna`

WGCNA 共表达网络分析。

| 简写 | 长参数 | 说明 |
|------|--------|------|
| `-e` | `--expression` | 必需 expression matrix；rows 必须是 genes，columns 是 samples。 |
| `-o` | `--output-dir` | 必需 WGCNA 输出目录。 |
| `-S` | `--coldata` | 可选 sample metadata 文件；sample ID 必须匹配 expression columns。 |
| `-G` | `--gene-info` | 可选 gene metadata 文件，按 gene ID 作为 index。 |
| `-sp` (`-sep`) | `--sep` | 可选分隔符覆盖；默认 `.csv` 用逗号，`.tsv`/`.txt` 用 tab。 |
| `-n` | `--name` | 存入 PyWGCNA object 的分析名称。默认：`WGCNA`。 |
| `-s` | `--species` | PyWGCNA enrichment analysis 使用的物种标签。 |
| `-l` | `--level` | 传给 PyWGCNA 的数据层级：`gene` 或 `transcript`。默认：`gene`。 |
| `-nw` (`-nt`) | `--network-type` | WGCNA network type：`unsigned`、`signed` 或 `signed hybrid`。 |
| `-tt` (`-tom`) | `--tom-type` | 传给 PyWGCNA 的 TOM type：`unsigned` 或 `signed`。 |
| `-ms` (`-min`) | `--min-module-size` | 模块检测的最小模块大小。默认：`50`。 |
| `-p` | `--power` | Soft-thresholding power；省略时由 PyWGCNA 自动检测。 |
| `-rs` (`-rsquared`) | `--rsquared-cut` | power 选择使用的 R-squared cutoff。默认：`0.9`。 |
| `-mc` (`-mean`) | `--mean-cut` | Mean connectivity cutoff。默认：`100`。 |
| `-md` (`-mediss`) | `--mediss-thresh` | 用于合并模块的 module eigengene dissimilarity threshold。默认：`0.2`。 |
| `-tc` (`-tpm`) | `--tpm-cutoff` | PyWGCNA 过滤使用的 TPM cutoff。默认：`1`。 |

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

## 输出结构

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

`manifest.json` 会记录 DESeq2 输入、解析后的 counts 文件、sample IDs、design、contrast、summary 和关键输出文件。

### WGCNA output

```text
wgcna_results/
├── figures/
├── WGCNA.p
└── module_info.csv
```
