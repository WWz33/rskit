import argparse
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from rskit import cli
from rskit.templates import write_template


class TemplateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.root = Path(self.tempdir.name)

    def test_write_template_writes_coldata_csv(self) -> None:
        output = self.root / "coldata.csv"

        written = write_template("coldata", str(output))

        table = pd.read_csv(written)
        self.assertEqual(
            list(table.columns),
            ["sample", "id", "condition", "r1", "r2"],
        )
        self.assertEqual(table.loc[0, "sample"], "sample1")

    def test_write_template_writes_contrast_tsv(self) -> None:
        output = self.root / "contrast.tsv"

        written = write_template("contrast", str(output))

        table = pd.read_csv(written, sep="\t")
        self.assertEqual(list(table.columns), ["factor", "level1", "level2"])
        self.assertEqual(table.loc[0, "factor"], "condition")

    def test_write_template_does_not_overwrite_without_force(self) -> None:
        output = self.root / "coldata.csv"
        output.write_text("existing\n", encoding="utf-8")

        with self.assertRaisesRegex(FileExistsError, "already exists"):
            write_template("coldata", str(output))

        self.assertEqual(output.read_text(encoding="utf-8"), "existing\n")

    def test_template_cli_uses_write_template(self) -> None:
        args = argparse.Namespace(
            template_name="coldata",
            output="coldata.csv",
            force=True,
        )

        with mock.patch("rskit.cli.write_template", return_value=Path("coldata.csv")) as write:
            cli.main_template(args)

        write.assert_called_once_with("coldata", "coldata.csv", force=True)


if __name__ == "__main__":
    unittest.main()
