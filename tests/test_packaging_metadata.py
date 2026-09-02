import ast
import tomllib
import unittest
from pathlib import Path

import rskit


ROOT = Path(__file__).resolve().parents[1]


class PackagingMetadataTests(unittest.TestCase):
    def test_pyproject_is_single_source_of_metadata(self) -> None:
        pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        project = pyproject["project"]

        self.assertEqual(project["name"], "rskit")
        self.assertEqual(project["scripts"]["rskit"], "rskit.cli:main")
        self.assertEqual(project["version"], rskit.__version__)

        # setup.py must stay a bare shim: no duplicated metadata
        setup_tree = ast.parse((ROOT / "setup.py").read_text(encoding="utf-8"))
        setup_call = next(
            node
            for node in ast.walk(setup_tree)
            if isinstance(node, ast.Call)
            and getattr(node.func, "id", None) == "setup"
        )
        self.assertEqual(setup_call.keywords, [])


if __name__ == "__main__":
    unittest.main()
