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

    def test_load_data_accepts_genes_as_expression_rows(self) -> None:
        expression = self.root / "expression.csv"
        pd.DataFrame(
            {"sample1": [1.0, 3.0], "sample2": [2.0, 4.0]},
            index=pd.Index(["geneA", "geneB"], name="gene_id"),
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
        self.assertEqual(list(passed_gene_expr.columns), ["geneA", "geneB"])
        self.assertEqual(list(passed_sample_info.index), ["sample1", "sample2"])

    def test_load_data_auto_detects_tsv_inputs(self) -> None:
        expression = self.root / "expression.tsv"
        pd.DataFrame(
            {"sample1": [1.0, 3.0], "sample2": [2.0, 4.0]},
            index=pd.Index(["geneA", "geneB"], name="gene_id"),
        ).to_csv(expression, sep="\t")
        coldata = self.root / "coldata.tsv"
        pd.DataFrame(
            {"condition": ["control", "treated"]},
            index=pd.Index(["sample1", "sample2"], name="sample"),
        ).to_csv(coldata, sep="\t")
        gene_info = self.root / "gene_info.tsv"
        pd.DataFrame(
            {"symbol": ["A", "B"]},
            index=pd.Index(["geneA", "geneB"], name="gene"),
        ).to_csv(gene_info, sep="\t")

        fake_wgcna = mock.Mock(return_value="wgcna-object")
        fake_module = types.SimpleNamespace(WGCNA=fake_wgcna)

        with mock.patch.dict(sys.modules, {"PyWGCNA": fake_module}):
            analyzer = WGCNAAnalyzer(output_dir=str(self.root / "out"))
            result = analyzer.load_data(
                str(expression),
                coldata=str(coldata),
                gene_info_file=str(gene_info),
            )

        self.assertEqual(result, "wgcna-object")
        passed_gene_expr = fake_wgcna.call_args.kwargs["geneExp"]
        passed_sample_info = fake_wgcna.call_args.kwargs["sampleInfo"]
        passed_gene_info = fake_wgcna.call_args.kwargs["geneInfo"]
        self.assertEqual(list(passed_gene_expr.columns), ["geneA", "geneB"])
        self.assertEqual(list(passed_sample_info.columns), ["condition"])
        self.assertEqual(list(passed_gene_info.index), ["geneA", "geneB"])

    def test_load_data_rejects_expression_in_samples_by_genes_orientation(self) -> None:
        expression = self.root / "expression.csv"
        pd.DataFrame(
            {"geneA": [1.0, 2.0], "geneB": [3.0, 4.0], "geneC": [5.0, 6.0]},
            index=pd.Index(["sample1", "sample2"], name="sample"),
        ).to_csv(expression)
        coldata = self._write_coldata()

        fake_wgcna = mock.Mock(return_value="wgcna-object")
        fake_module = types.SimpleNamespace(WGCNA=fake_wgcna)

        with mock.patch.dict(sys.modules, {"PyWGCNA": fake_module}):
            analyzer = WGCNAAnalyzer(output_dir=str(self.root / "out"))
            with self.assertRaisesRegex(ValueError, "transpose"):
                analyzer.load_data(str(expression), coldata=str(coldata))

    def test_load_data_rejects_missing_metadata_samples(self) -> None:
        expression = self.root / "expression.csv"
        pd.DataFrame(
            {"sample1": [1.0, 3.0]},
            index=pd.Index(["geneA", "geneB"], name="gene_id"),
        ).to_csv(expression)
        coldata = self._write_coldata()

        analyzer = WGCNAAnalyzer(output_dir=str(self.root / "out"))
        with self.assertRaisesRegex(ValueError, "sample2"):
            analyzer.load_data(str(expression), coldata=str(coldata))

    def test_load_data_passes_power_to_pywgcna(self) -> None:
        expression = self.root / "expression.csv"
        pd.DataFrame(
            {"sample1": [1.0, 3.0], "sample2": [2.0, 4.0]},
            index=pd.Index(["geneA", "geneB"], name="gene_id"),
        ).to_csv(expression)
        coldata = self._write_coldata()

        fake_wgcna = mock.Mock(return_value="wgcna-object")
        fake_module = types.SimpleNamespace(WGCNA=fake_wgcna)

        with mock.patch.dict(sys.modules, {"PyWGCNA": fake_module}):
            WGCNAAnalyzer(output_dir=str(self.root / "out"), power=3).load_data(
                str(expression), coldata=str(coldata)
            )
        self.assertEqual(fake_wgcna.call_args.kwargs["powers"], [3])

        with mock.patch.dict(sys.modules, {"PyWGCNA": fake_module}):
            WGCNAAnalyzer(output_dir=str(self.root / "out")).load_data(
                str(expression), coldata=str(coldata)
            )
        self.assertIsNone(fake_wgcna.call_args.kwargs["powers"])

    def test_load_data_rejects_square_matrix_with_samples_as_rows(self) -> None:
        expression = self.root / "expression.csv"
        # square 2x2: sample IDs appear as both rows and columns; samples-as-rows
        # wins the orientation check, so this must be rejected, not silently kept
        pd.DataFrame(
            {"sample1": [1.0, 2.0], "sample2": [3.0, 4.0]},
            index=pd.Index(["sample1", "sample2"], name="sample"),
        ).to_csv(expression)
        coldata = self._write_coldata()

        fake_wgcna = mock.Mock(return_value="wgcna-object")
        fake_module = types.SimpleNamespace(WGCNA=fake_wgcna)

        with mock.patch.dict(sys.modules, {"PyWGCNA": fake_module}):
            analyzer = WGCNAAnalyzer(output_dir=str(self.root / "out"))
            with self.assertRaisesRegex(ValueError, "transpose"):
                analyzer.load_data(str(expression), coldata=str(coldata))


if __name__ == "__main__":
    unittest.main()
