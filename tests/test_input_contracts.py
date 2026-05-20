import tempfile
import unittest
from pathlib import Path

import pandas as pd

from rskit.input_contracts import (
    detect_separator,
    load_coldata,
    read_table,
    validate_sample_alignment,
)


class InputContractTests(unittest.TestCase):
    def test_detect_separator_uses_extension_unless_overridden(self) -> None:
        self.assertEqual(detect_separator("samples.csv"), ",")
        self.assertEqual(detect_separator("samples.tsv"), "\t")
        self.assertEqual(detect_separator("samples.txt"), "\t")
        self.assertEqual(detect_separator("samples.tsv", sep=","), ",")

    def test_read_table_uses_detected_separator(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            table_path = Path(tempdir) / "samples.tsv"
            table_path.write_text("sample\tcondition\ns1\tA\n", encoding="utf-8")

            table = read_table(str(table_path))

        self.assertEqual(table.to_dict("records"), [{"sample": "s1", "condition": "A"}])

    def test_load_coldata_requires_sample_column_and_requested_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            coldata_path = Path(tempdir) / "coldata.csv"
            coldata_path.write_text("sample,condition\ns1,A\ns2,B\n", encoding="utf-8")

            coldata = load_coldata(str(coldata_path), required_columns=["condition"])

        self.assertEqual(list(coldata.index), ["s1", "s2"])
        self.assertEqual(list(coldata.columns), ["condition"])

    def test_load_coldata_rejects_missing_sample_column(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            coldata_path = Path(tempdir) / "coldata.csv"
            coldata_path.write_text("id,condition\ns1,A\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "sample"):
                load_coldata(str(coldata_path))

    def test_validate_sample_alignment_requires_exact_match(self) -> None:
        table = pd.DataFrame({"geneA": [1, 2]}, index=["s1", "s3"])
        metadata = pd.DataFrame({"condition": ["A", "B"]}, index=["s1", "s2"])

        with self.assertRaisesRegex(ValueError, "missing from counts"):
            validate_sample_alignment(table, metadata, table_name="counts")


if __name__ == "__main__":
    unittest.main()
