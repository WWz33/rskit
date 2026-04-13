import argparse
import csv
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from rskit import cli
from rskit.core.deseq2 import run_deseq2_cli
from rskit.core.salmon import SalmonExpressionExporter


class _FakeArray:
    def __init__(self, frame: pd.DataFrame):
        self._frame = frame

    def to_pandas(self) -> pd.DataFrame:
        return self._frame.copy()


class _FakeDataset(dict):
    pass


class QuantExpressionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.root = Path(self.tempdir.name)

    def _write_quant_stub(self, sample_name: str) -> None:
        sample_dir = self.root / "03_quant" / sample_name
        sample_dir.mkdir(parents=True, exist_ok=True)
        (sample_dir / "quant.sf").write_text(
            "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
            "tx1\t100\t80\t1.0\t10\n",
            encoding="utf-8",
        )

    def _write_gtf(self) -> Path:
        gtf_path = self.root / "annotation.gtf"
        gtf_path.write_text(
            "\n".join(
                [
                    'chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id "geneA"; transcript_id "tx1";',
                    'chr1\tsrc\ttranscript\t101\t200\t.\t+\t.\tgene_id "geneB"; transcript_id "tx2";',
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        return gtf_path

    def test_exporter_writes_gene_expression_tables(self) -> None:
        self._write_quant_stub("sample1")
        self._write_quant_stub("sample2")
        gtf_path = self._write_gtf()

        counts = pd.DataFrame(
            {"0": [10.0, 20.0], "1": [30.0, 40.0]},
            index=pd.Index(["geneA", "geneB"], name="gene_id"),
        )
        abundance = pd.DataFrame(
            {"0": [1.0, 4.0], "1": [3.0, 15.0]},
            index=pd.Index(["geneA", "geneB"], name="gene_id"),
        )
        fake_module = types.SimpleNamespace(
            tximport=lambda **kwargs: _FakeDataset(
                counts=_FakeArray(counts),
                abundance=_FakeArray(abundance),
            )
        )

        with mock.patch.dict(sys.modules, {"pytximport": fake_module}):
            exporter = SalmonExpressionExporter()
            outputs = exporter.export_gene_tables(
                salmon_dir=str(self.root / "03_quant"),
                output_dir=str(self.root / "03_quant"),
                gtf_file=str(gtf_path),
                sample_names=["sample1", "sample2"],
            )

        counts_df = pd.read_csv(outputs["gene_counts"], index_col=0)
        tpm_df = pd.read_csv(outputs["gene_tpm"], index_col=0)
        log_tpm_df = pd.read_csv(outputs["gene_log2_tpm_plus1"], index_col=0)
        tx2gene_df = pd.read_csv(outputs["tx2gene"], sep="\t")

        self.assertEqual(list(counts_df.index), ["sample1", "sample2"])
        self.assertEqual(list(counts_df.columns), ["geneA", "geneB"])
        self.assertEqual(counts_df.loc["sample1", "geneA"], 10.0)
        self.assertEqual(tpm_df.loc["sample2", "geneB"], 15.0)
        self.assertAlmostEqual(log_tpm_df.loc["sample1", "geneB"], 2.321928, places=5)
        self.assertEqual(tx2gene_df.to_dict("records"), [
            {"transcript_id": "tx1", "gene_id": "geneA"},
            {"transcript_id": "tx2", "gene_id": "geneB"},
        ])

    def test_run_deseq2_cli_prefers_quant_gene_counts_file(self) -> None:
        quant_dir = self.root / "03_quant"
        quant_dir.mkdir(parents=True, exist_ok=True)
        precomputed_counts = quant_dir / "gene_counts.csv"
        pd.DataFrame({"geneA": [10, 12]}, index=["sample1", "sample2"]).to_csv(precomputed_counts)

        coldata_path = self.root / "coldata.csv"
        with coldata_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["sample", "condition"])
            writer.writeheader()
            writer.writerow({"sample": "sample1", "condition": "A"})
            writer.writerow({"sample": "sample2", "condition": "B"})

        output_dir = self.root / "04_deseq2"
        args = argparse.Namespace(
            salmon_dir=str(quant_dir),
            gene_counts=None,
            coldata=str(coldata_path),
            gtf=None,
            tx2gene=None,
            output_dir=str(output_dir),
            design="~condition",
            alpha=0.05,
            lfc_threshold=2.0,
            threads=None,
            contrast=None,
        )

        fake_results = pd.DataFrame(
            {
                "baseMean": [1.0],
                "log2FoldChange": [0.5],
                "pvalue": [0.1],
                "padj": [0.2],
            },
            index=["geneA"],
        )

        with mock.patch("rskit.core.deseq2.Deseq2Analyzer.load_counts_from_file", return_value=pd.DataFrame({"geneA": [10, 12]}, index=["sample1", "sample2"])) as load_counts_from_file, \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.load_counts_from_salmon") as load_counts_from_salmon, \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.analyze", return_value=fake_results), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.save_results", return_value={}), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.plot_volcano"), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.plot_pca"), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.get_summary", return_value={
                 "total_genes": 1,
                 "significant_genes": 0,
                 "upregulated_genes": 0,
                 "downregulated_genes": 0,
                 "alpha": 0.05,
             }):
            run_deseq2_cli(args)

        load_counts_from_file.assert_called_once_with(str(precomputed_counts))
        load_counts_from_salmon.assert_not_called()

    def test_main_all_passes_quant_gene_counts_into_deseq(self) -> None:
        coldata_path = self.root / "coldata.csv"
        with coldata_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["sample", "id", "condition", "r1", "r2"])
            writer.writeheader()
            writer.writerow({
                "sample": "sample1",
                "id": "group1",
                "condition": "A",
                "r1": str(self.root / "sample1_R1.fq.gz"),
                "r2": str(self.root / "sample1_R2.fq.gz"),
            })

        args = argparse.Namespace(
            coldata=str(coldata_path),
            genome_fasta=str(self.root / "genome.fa"),
            gtf_file=str(self.root / "annotation.gtf"),
            transcript_fasta=str(self.root / "transcripts.fa"),
            output_dir=str(self.root / "results"),
            index_dir=None,
            tx2gene=None,
            threads=4,
            parallel=None,
            trim=False,
            force_index=False,
            skip_existing=False,
            design="~condition",
            contrast=None,
            alpha=0.05,
            lfc_threshold=2.0,
        )

        exported = {"gene_counts": str(self.root / "results" / "03_quant" / "gene_counts.csv")}

        with mock.patch("rskit.cli.build_index_if_needed"), \
             mock.patch("rskit.cli.prepare_samples", return_value={"sample1": {"fq1": "a", "fq2": "b"}}), \
             mock.patch("rskit.cli.run_quantification", return_value={"sample1": {}}), \
             mock.patch("rskit.cli.export_quant_expression_tables", return_value=exported), \
             mock.patch("rskit.cli.run_deseq2_cli") as run_deseq2:
            cli.main_all(args)

        passed_args = run_deseq2.call_args.args[0]
        self.assertEqual(passed_args.gene_counts, exported["gene_counts"])
        self.assertIsNone(passed_args.salmon_dir)


if __name__ == "__main__":
    unittest.main()
