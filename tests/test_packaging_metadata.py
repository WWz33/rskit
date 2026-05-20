import ast
import tomllib
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class PackagingMetadataTests(unittest.TestCase):
    def test_setup_py_matches_pyproject_name_and_cli(self) -> None:
        pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        setup_tree = ast.parse((ROOT / "setup.py").read_text(encoding="utf-8"))

        setup_call = next(
            node
            for node in ast.walk(setup_tree)
            if isinstance(node, ast.Call)
            and getattr(node.func, "id", None) == "setup"
        )
        setup_kwargs = {
            keyword.arg: ast.literal_eval(keyword.value)
            for keyword in setup_call.keywords
            if keyword.arg in {"name", "entry_points"}
        }

        self.assertEqual(setup_kwargs["name"], pyproject["project"]["name"])
        self.assertEqual(
            setup_kwargs["entry_points"]["console_scripts"],
            ["rskit=rskit.cli:main"],
        )
        self.assertEqual(pyproject["project"]["scripts"]["rskit"], "rskit.cli:main")


if __name__ == "__main__":
    unittest.main()
