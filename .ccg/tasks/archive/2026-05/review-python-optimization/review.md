# Python Optimization Review

## Verified

- `python -m compileall -q rskit` passes.
- `python -m pytest -q` passes: 3 tests.
- Bare `pytest -q` fails in this environment with `ModuleNotFoundError: No module named 'rskit'`; use `python -m pytest` or install the package editable for local development.
- `python setup.py --name` cannot run because this Python environment lacks `setuptools`.
- `pyproject.toml` parses and declares project `rskit` with script `rskit = rskit.cli:main`.
- Dependency architecture was checked with `C:\Program Files\GitHub CLI\gh.exe` and DeepWiki, per user request.

## Dependency Findings

- `scverse/PyDESeq2` latest release from GitHub: `v0.5.4` on 2026-01-23. DeepWiki and upstream example `examples/plot_pandas_io_example.py` both confirm `DeseqDataSet` expects counts as samples x genes and metadata as samples x variables. This matches `rskit.core.deseq2` after `load_counts_from_file()` or `SalmonExpressionExporter` produces samples x genes.
- `complextissue/pytximport` latest release from GitHub: `0.13.0` on 2026-02-25. Upstream `tximport()` accepts `file_paths`, `data_type="salmon"`, a `transcript_gene_map` DataFrame/file with `transcript_id` and `gene_id`, and `counts_from_abundance="length_scaled_tpm"`. For `output_type="xarray"`, upstream arrays are gene/transcript x file; `rskit.core.salmon._to_dataframe(...).T` is therefore directionally correct for PyDESeq2.
- `mortazavilab/PyWGCNA` latest release from GitHub: `V2.2.1` on 2024-11-05. DeepWiki and upstream `PyWGCNA/wgcna.py` confirm `geneExp` should have samples as rows and genes as columns, with `sampleInfo` index matching samples and `geneInfo` index matching gene columns. `rskit.core.wgcna` passes the expression file directly with `index_col=0`, so documentation/tests must make this orientation explicit.
- `COMBINE-lab/salmon` latest release from GitHub: `v1.11.4` on 2026-03-11. Upstream docs confirm alignment-based quantification uses `salmon quant -t transcripts.fa -l <LIBTYPE> -a aln.bam -o salmon_quant`; Salmon also notes alignment-mode speed usually saturates around 8-12 threads. DeepWiki confirms `--validateMappings` does not apply to alignment/BAM mode.

## Highest Priority

1. Packaging metadata conflict:
   - `pyproject.toml` declares `name = "rskit"` and `rskit = "rskit.cli:main"`.
   - `setup.py` declares `name="rnaseq-tools"` and `rnaseq-tools=rnaseq_tools.cli:main`, but no `rnaseq_tools` package exists.
   - Recommendation: remove `setup.py` if `pyproject.toml` is canonical, or make it delegate consistently to `rskit`. This is a correctness issue, not style.

2. CLI module is doing too much:
   - `rskit/cli.py` contains parser construction, working directory creation, FASTQ trimming, sample parsing, index handling, quantification orchestration, and complete pipeline orchestration.
   - Recommendation: keep parser construction in `cli.py`, move reusable workflow functions into focused modules such as `rskit/core/workflow.py` or `rskit/core/samples.py`. This will make unit tests cheaper and reduce import-time coupling.

3. Library code exits the interpreter:
   - `rskit/core/wgcna.py` calls `sys.exit(1)` from library methods when PyWGCNA is missing or data is not loaded.
   - Recommendation: raise `ImportError` or `RuntimeError` in library code; only CLI entry points should decide to exit.

4. Salmon alignment-mode options should reflect upstream Salmon:
   - `rskit.core.salmon.SalmonQuantifier` uses alignment input (`-a`) and optionally appends `--validateMappings`.
   - Upstream Salmon v1.11.4 docs describe `--validateMappings` as mapping/selective-alignment behavior, not alignment/BAM quantification behavior.
   - Recommendation: remove `validate_mappings` from the alignment-mode command path, or split Salmon support into explicit mapping-mode and alignment-mode config objects. Also consider capping default alignment-mode threads to a practical default such as 8-12 instead of 56.

5. External command wrapper loses useful failure context:
   - `rskit/core/base.py` captures stdout/stderr but only logs stderr on failure.
   - Recommendation: include command, return code, stdout, and stderr in a structured exception. Tests should assert error message content for failed commands.

## Medium Priority

6. Count matrix orientation heuristic is fragile:
   - `rskit/core/deseq2.py` transposes when row count is more than twice column count.
   - Dependency check confirms PyDESeq2 requires samples x genes, and pytximport conversion already produces that orientation in this package.
   - Recommendation: use metadata sample names to validate orientation where available; if ambiguous, fail with a clear error instead of guessing.

7. GTF utility shadows built-in `open` and has a broken CLI helper:
   - `rskit/utils/gtf.py` defines `open(...)`.
   - `gtf_tx2gene()` then calls `open(args.output, "w")`, which resolves to the local parser function, and references `gtf.open` without importing `gtf`.
   - Recommendation: rename parser entry point to `iter_records` or explicitly use `builtins.open`; add tests for `gtf_tx2gene` or remove it if unused.

8. PyWGCNA input orientation should be validated:
   - Upstream PyWGCNA expects samples as rows and genes as columns.
   - `rskit.core.wgcna` reads the expression file and passes it straight through.
   - Recommendation: validate that expression rows match `sampleInfo` rows when coldata is provided, and add a clear error if the user supplies genes x samples.

9. Type hints are partial and old-style:
   - Several public functions return bare `dict`/`Dict`, accept untyped args, or use `Optional`/`List` while the project could move toward Python 3.10+ syntax.
   - Recommendation: first decide supported Python version. Current metadata says `>=3.8`, so do not use `X | None` unless raising to `>=3.10`.

10. Tests cover recent Salmon/DESeq2 flow but not command construction:
   - Add tests for STAR command construction, Salmon command construction, sample coldata parsing, missing column errors, and WGCNA missing dependency behavior.
   - Mock `subprocess.run`; do not require STAR, Salmon, or PyWGCNA in unit tests.

## Low Priority

11. Remove unused imports:
   - Examples include `os` in `rskit/cli.py` and `rskit/core/wgcna.py`, `SalmonQuantifier` in `rskit/cli.py`, and unused typing imports in some modules.

12. Add dev tooling to `pyproject.toml`:
    - Add optional dependencies for tests/dev, plus `ruff` and formatter configuration.
    - Keep this separate from behavior changes.

## Suggested Order

1. Fix packaging metadata and editable install/test invocation.
2. Add tests around current command construction and GTF helper behavior.
3. Fix Salmon alignment-mode config so it no longer exposes ineffective `--validateMappings`.
4. Add PyWGCNA orientation validation.
5. Fix GTF helper and library `sys.exit` behavior.
6. Improve command failure diagnostics.
7. Split `cli.py` only after tests protect behavior.
8. Tighten types and add ruff/mypy incrementally.
