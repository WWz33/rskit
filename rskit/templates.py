from pathlib import Path
from typing import Dict, List

import pandas as pd

from rskit.input_contracts import detect_separator


TEMPLATE_ROWS: Dict[str, List[Dict[str, str]]] = {
    "coldata": [
        {
            "sample": "sample1",
            "id": "ctrl",
            "condition": "control",
            "r1": "reads/sample1_R1.fq.gz",
            "r2": "reads/sample1_R2.fq.gz",
        },
        {
            "sample": "sample2",
            "id": "treat",
            "condition": "treatment",
            "r1": "reads/sample2_R1.fq.gz",
            "r2": "reads/sample2_R2.fq.gz",
        },
    ],
    "contrast": [
        {
            "factor": "condition",
            "level1": "treatment",
            "level2": "control",
        },
    ],
}


def write_template(template_name: str, output: str, force: bool = False) -> Path:
    """Write an input template file and return its path."""
    if template_name not in TEMPLATE_ROWS:
        raise ValueError("Unknown template: " + template_name)

    output_path = Path(output)
    if output_path.exists() and not force:
        raise FileExistsError(f"Output file already exists: {output_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(TEMPLATE_ROWS[template_name])
    table.to_csv(output_path, sep=detect_separator(str(output_path)), index=False)
    return output_path
