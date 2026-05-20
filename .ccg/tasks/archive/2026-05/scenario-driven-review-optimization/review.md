# User-Scenario-Driven Review Framework

## Core Idea

Review from user intent outward:

1. Define the user scenario.
2. Define what "success" means for that user.
3. Identify where the workflow can fail.
4. Map each failure point to code, dependency contracts, CLI behavior, tests, and docs.
5. Prioritize by user harm, not by code aesthetics.

## rskit Primary Scenarios

1. First-time local user:
   - Installs package.
   - Runs `rskit --help`.
   - Tries one small demo or existing example.
   - Biggest risks: packaging mismatch, missing dependencies, unclear install docs.

2. RNA-seq pipeline user:
   - Runs `rskit quant` or `rskit all`.
   - Uses STAR, Salmon, GTF, coldata, and output directories.
   - Biggest risks: bad paths, wrong sample sheet columns, external command failures, ineffective Salmon options.

3. Differential expression user:
   - Starts from Salmon quant output or gene count matrix.
   - Runs PyDESeq2.
   - Biggest risks: counts orientation, metadata/sample mismatch, contrast inference, ambiguous matrix transposition.

4. WGCNA user:
   - Starts from expression matrix and optional metadata.
   - Runs PyWGCNA.
   - Biggest risks: rows/columns orientation, library code calling `sys.exit`, memory/runtime surprises.

5. Maintainer/developer:
   - Runs tests, edits modules, releases package.
   - Biggest risks: weak test harness, missing dev dependencies, conflicting `setup.py` and `pyproject.toml`.

## Review Lenses

- Installability: Can a new user install and import the package cleanly?
- CLI contract: Are help text, required flags, defaults, and runtime behavior aligned?
- Data contract: Are matrix orientations and required columns validated before expensive work starts?
- Dependency contract: Does code match upstream APIs and current behavior?
- External-tool contract: Are STAR/Salmon/fastp commands correct and diagnosable?
- Recovery: Can the user understand and resume after failure?
- Reproducibility: Are versions, parameters, outputs, and logs enough to rerun or debug?
- Testability: Can each scenario be verified without requiring heavy bioinformatics tools in unit tests?

## Priority Rule

Prioritize by this order:

1. Blocks installation or CLI startup.
2. Produces wrong scientific results silently.
3. Wastes hours of compute before failing.
4. Gives unclear errors for common user mistakes.
5. Makes maintenance or release fragile.
6. Cosmetic style issues.

## Recommended Review Output

For each scenario, write:

- User goal.
- Happy path command.
- Expected inputs and outputs.
- Failure classes.
- Code locations responsible.
- Upstream dependency assumptions.
- Tests that should exist.
- Optimization actions ranked by user impact.
