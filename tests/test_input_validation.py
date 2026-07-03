import argparse
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from rskit import cli
from rskit.input_validation import validate_input_files


class InputValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.root = Path(self.tempdir.name)

    def _write_coldata(self) -> Path:
        reads_dir = self.root / "reads"
        reads_dir.mkdir()
        (reads_dir / "sample1_R1.fq.gz").write_text("", encoding="utf-8")
        (reads_dir / "sample1_R2.fq.gz").write_text("", encoding="utf-8")
        (reads_dir / "sample2_R1.fq.gz").write_text("", encoding="utf-8")
        (reads_dir / "sample2_R2.fq.gz").write_text("", encoding="utf-8")

        coldata = self.root / "coldata.csv"
        coldata.write_text(
            "sample,condition,r1,r2\n"
            "sample1,A,reads/sample1_R1.fq.gz,reads/sample1_R2.fq.gz\n"
            "sample2,B,reads/sample2_R1.fq.gz,reads/sample2_R2.fq.gz\n",
            encoding="utf-8",
        )
        return coldata

    def test_validate_input_files_accepts_reads_and_gene_by_sample_counts(self) -> None:
        coldata = self._write_coldata()
        counts = self.root / "counts.csv"
        pd.DataFrame(
            {"sample1": [10, 20], "sample2": [12, 24]},
            index=["geneA", "geneB"],
        ).to_csv(counts)

        messages = validate_input_files(
            coldata=str(coldata),
            check_reads=True,
            gene_counts=str(counts),
        )

        self.assertEqual(messages[0], "coldata: 2 samples")
        self.assertIn("reads: r1/r2 paths exist", messages)
        self.assertIn("gene counts: 2 samples x 2 genes", messages)

    def test_validate_input_files_requires_design_columns(self) -> None:
        coldata = self.root / "coldata.csv"
        coldata.write_text("sample,condition\nsample1,A\n", encoding="utf-8")

        with self.assertRaisesRegex(ValueError, "batch"):
            validate_input_files(str(coldata), design="~batch + condition")

    def test_validate_input_files_reports_missing_reads(self) -> None:
        coldata = self.root / "coldata.csv"
        coldata.write_text(
            "sample,condition,r1,r2\n"
            "sample1,A,missing_R1.fq.gz,missing_R2.fq.gz\n",
            encoding="utf-8",
        )

        with self.assertRaisesRegex(ValueError, "Read files do not exist"):
            validate_input_files(str(coldata), check_reads=True)

    def test_validate_input_files_accepts_gene_by_sample_expression(self) -> None:
        coldata = self.root / "coldata.csv"
        coldata.write_text(
            "sample,condition\nsample1,A\nsample2,B\n",
            encoding="utf-8",
        )
        expression = self.root / "expression.csv"
        pd.DataFrame(
            {"sample1": [1.0, 2.0], "sample2": [3.0, 4.0]},
            index=["geneA", "geneB"],
        ).to_csv(expression)

        messages = validate_input_files(str(coldata), expression=str(expression))

        self.assertIn("expression: 2 samples x 2 genes", messages)

    def test_validate_input_files_warns_on_samples_by_genes_expression(self) -> None:
        coldata = self.root / "coldata.csv"
        coldata.write_text(
            "sample,condition\nsample1,A\nsample2,B\n",
            encoding="utf-8",
        )
        expression = self.root / "expression.csv"
        pd.DataFrame(
            {"geneA": [1.0, 2.0], "geneB": [3.0, 4.0], "geneC": [5.0, 6.0]},
            index=["sample1", "sample2"],
        ).to_csv(expression)

        with self.assertLogs("rskit.input_validation", level="WARNING") as logs:
            messages = validate_input_files(str(coldata), expression=str(expression))

        self.assertIn("genes x samples", "\n".join(logs.output))
        self.assertIn("expression: 2 samples x 3 genes", messages)

    def test_doctor_cli_uses_validate_input_files(self) -> None:
        args = argparse.Namespace(
            coldata="coldata.csv",
            design="~condition",
            check_reads=True,
            gene_counts="counts.csv",
            expression=None,
        )

        with mock.patch("rskit.cli.validate_input_files", return_value=["coldata: 2 samples"]) as validate:
            cli.main_validate(args)

        validate.assert_called_once_with(
            coldata="coldata.csv",
            design="~condition",
            check_reads=True,
            gene_counts="counts.csv",
            expression=None,
        )


if __name__ == "__main__":
    unittest.main()
