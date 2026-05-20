import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from rskit.core.wgcna import WGCNAAnalyzer


class WGCNAAnalyzerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.root = Path(self.tempdir.name)

    def _write_coldata(self) -> Path:
        coldata = self.root / "coldata.csv"
        pd.DataFrame(
            {"condition": ["control", "treated"]},
            index=pd.Index(["sample1", "sample2"], name="sample"),
        ).to_csv(coldata)
        return coldata

    def test_load_data_accepts_samples_as_expression_rows(self) -> None:
        expression = self.root / "expression.csv"
        pd.DataFrame(
            {"geneA": [1.0, 2.0], "geneB": [3.0, 4.0]},
            index=pd.Index(["sample1", "sample2"], name="sample"),
        ).to_csv(expression)
        coldata = self._write_coldata()

        fake_wgcna = mock.Mock(return_value="wgcna-object")
        fake_module = types.SimpleNamespace(WGCNA=fake_wgcna)

        with mock.patch.dict(sys.modules, {"PyWGCNA": fake_module}):
            analyzer = WGCNAAnalyzer(output_dir=str(self.root / "out"))
            result = analyzer.load_data(str(expression), coldata=str(coldata))

        self.assertEqual(result, "wgcna-object")
        passed_gene_expr = fake_wgcna.call_args.kwargs["geneExp"]
        passed_sample_info = fake_wgcna.call_args.kwargs["sampleInfo"]
        self.assertEqual(list(passed_gene_expr.index), ["sample1", "sample2"])
        self.assertEqual(list(passed_sample_info.index), ["sample1", "sample2"])

    def test_load_data_rejects_transposed_expression_matrix(self) -> None:
        expression = self.root / "expression.csv"
        pd.DataFrame(
            {"sample1": [1.0, 2.0], "sample2": [3.0, 4.0]},
            index=pd.Index(["geneA", "geneB"], name="gene"),
        ).to_csv(expression)
        coldata = self._write_coldata()

        analyzer = WGCNAAnalyzer(output_dir=str(self.root / "out"))
        with self.assertRaisesRegex(ValueError, "appears transposed"):
            analyzer.load_data(str(expression), coldata=str(coldata))

    def test_load_data_rejects_missing_metadata_samples(self) -> None:
        expression = self.root / "expression.csv"
        pd.DataFrame(
            {"geneA": [1.0], "geneB": [3.0]},
            index=pd.Index(["sample1"], name="sample"),
        ).to_csv(expression)
        coldata = self._write_coldata()

        analyzer = WGCNAAnalyzer(output_dir=str(self.root / "out"))
        with self.assertRaisesRegex(ValueError, "sample2"):
            analyzer.load_data(str(expression), coldata=str(coldata))


if __name__ == "__main__":
    unittest.main()
