# Requirements

User asked to inspect the repository code and identify how to optimize it using Python-oriented review skills.

Success criteria:
- Inspect package structure, packaging metadata, core modules, utilities, and tests.
- Identify concrete optimization opportunities with file/line evidence.
- Prefer correctness, maintainability, and project structure improvements over speculative micro-optimizations.
- Verify current test/compile status with local commands.
- Do not change product code during this review pass.

Constraints:
- No `.ccg/spec` directory exists in this repository.
- CCG external dual-model wrapper `~/.claude/bin/codeagent-wrapper` is unavailable in this environment, so Gemini/Claude analysis and review could not be run.
