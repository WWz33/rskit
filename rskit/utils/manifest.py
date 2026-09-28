"""Run manifests: JSON records of inputs, parameters, and outputs.

Written next to DESeq2 results and in the quant output directory so a
finished run can be reproduced or cited without re-reading console logs.
"""

import json
from pathlib import Path
from typing import Dict

import numpy as np


def _json_safe(value):
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


def write_manifest(output_dir, manifest: Dict) -> Path:
    """Write a JSON run manifest and return its path."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    manifest_path = output_path / "manifest.json"
    manifest_path.write_text(
        json.dumps(_json_safe(manifest), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest_path
