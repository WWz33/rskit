# Changelog

All notable changes to rskit are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and versions follow
semantic versioning.

## [0.3.0] - 2026-10-01

### Added

- `--salmon-direct`: quantify from the reads with Salmon's selective alignment,
  skipping STAR, the genome, and the index. `-g`/`-gtf` become optional; `all`
  still needs an annotation so DESeq2 has gene-level counts.
- Multiple `-c/--contrast` flags run several contrasts against a single DESeq2
  model fit. Each contrast gets its own `<factor>_<level1>_vs_<level2>/`
  subdirectory, `deseq2_significant_all.csv` merges them with a `contrast`
  column, and the manifest lists every contrast. A single `-c` keeps the
  previous flat layout.
- `03_quant/manifest.json` records the quantification inputs, samples,
  parameters, tool versions, and the STAR index fingerprint.
- STAR index fingerprint (`00_index/.rskit_index.json`): rskit warns when the
  genome FASTA or GTF changed after the index was built instead of silently
  aligning against a stale index.
- `00_summary/summary.csv` gains `salmon_expected_format`, surfacing Salmon's
  detected library type per sample.
- `--keep-going` (continue after a failed sample, then report all failures),
  `--dry-run` (print STAR/Salmon/fastp commands without running them), and
  `--verbose` (full traceback on failure).
- Optional-dependency extras: `pip install -e ".[deseq2]"`, `".[wgcna]"`, or
  `".[all]"`; a minimal install still runs the quantification phase.
- GitHub Actions CI (lint + Python 3.11/3.12/3.13 + minimal-install leg) and an
  integration test that runs real STAR and fastp on a tiny genome when they are
  on PATH.

### Changed

- `--skip-existing` now reuses every completed artifact: `quant.sf` skips the
  sample, a finished transcriptome BAM (with STAR's `Log.final.out` as the
  completion marker) skips alignment, and clean reads skip fastp. Fresh
  gene-level tables are reused instead of re-running tximport.
- fastp writes `.fq.gz` instead of uncompressed `.fq`.
- The thread budget for index building and sample scheduling is capped to the
  `sched_getaffinity` allocation, matching what DESeq2 inference already did.
- `Deseq2Analyzer.analyze()` is split into `fit()` and `contrast_results()` so
  API callers can share one fit across contrasts; `analyze()` remains as the
  single-contrast wrapper.
- Tool preflight checks STAR/salmon/fastp up front, and user-facing errors exit
  with a one-line message (exit code 1) instead of a traceback.

### Fixed

- STAR index completeness check required `chrNameLength`, but STAR writes
  `chrNameLength.txt`; valid indexes never passed, so sequential runs rebuilt
  the index on every invocation and parallel runs failed right after building
  it.
- GFF3 annotations are usable: tx2gene generation accepts `mRNA` records (GFF3)
  in addition to `transcript` (GTF).
- `--skip-existing` honoured during trimming (previously fastp ran for every
  sample even when the work was done).
- Empty `sample` cells and missing values in design columns are rejected with
  clear messages instead of becoming a sample literally named `nan` or failing
  inside pydeseq2.
- Headerless `tx2gene` files with non-Ensembl IDs no longer lose their first
  row to the header.
- `~` in coldata read paths expands to `$HOME`; `--alpha`/`--lfc` are
  range-checked; NaN count matrices are rejected early.
- `run_with_deseq2()` propagates DESeq2 failures instead of returning them in a
  success-looking result dict.
- Volcano plots clip `pvalue == 0` instead of losing the whole plot.
- The QC summary is also written when a run fails mid-way, for the samples that
  finished.

## [0.2.0]

Earlier release: CLI and Python API for STAR/Salmon quantification, DESeq2
differential expression, and WGCNA, with input preflight checks, coldata
templates, and passthrough arguments for STAR/Salmon/fastp. Recorded here
without a detailed change log.
