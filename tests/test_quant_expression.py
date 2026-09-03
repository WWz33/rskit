import argparse
import csv
import json
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from rskit import cli
from rskit.config import DESeq2Config, PipelineConfig, SalmonConfig, StarConfig
from rskit.core.deseq2 import (
    Deseq2Analyzer,
    _lfc_shrink_coefficient,
    parse_contrast,
    run_deseq2_cli,
)
from rskit.core.pipeline import RNAseqPipeline
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

    def _patch_pydeseq2_modules(self, captured: dict):
        pydeseq2_module = types.ModuleType("pydeseq2")
        dds_module = types.ModuleType("pydeseq2.dds")
        ds_module = types.ModuleType("pydeseq2.ds")
        inference_module = types.ModuleType("pydeseq2.default_inference")

        class FakeDefaultInference:
            def __init__(self, n_cpus=None):
                captured["n_cpus"] = n_cpus

        class FakeDeseqDataSet:
            def __init__(self, **kwargs):
                self.counts = kwargs["counts"].copy()
                self.metadata = kwargs["metadata"].copy()
                self.obsm = {"design_matrix": pd.DataFrame(index=self.metadata.index)}
                self.layers = {}
                captured["counts"] = self.counts.copy()
                captured["metadata"] = self.metadata.copy()

            def deseq2(self):
                captured["ran_deseq2"] = True

        class FakeDeseqStats:
            def __init__(self, dds, **kwargs):
                self.dds = dds
                self.results_df = pd.DataFrame(
                    {
                        "baseMean": [1.0] * len(dds.counts.columns),
                        "log2FoldChange": [0.0] * len(dds.counts.columns),
                        "pvalue": [1.0] * len(dds.counts.columns),
                        "padj": [1.0] * len(dds.counts.columns),
                    },
                    index=dds.counts.columns,
                )
                captured["contrast"] = kwargs["contrast"]

            def summary(self):
                captured["summary"] = True

            def lfc_shrink(self, coeff):
                captured["lfc_coeff"] = coeff

        dds_module.DeseqDataSet = FakeDeseqDataSet
        ds_module.DeseqStats = FakeDeseqStats
        inference_module.DefaultInference = FakeDefaultInference

        return mock.patch.dict(
            sys.modules,
            {
                "pydeseq2": pydeseq2_module,
                "pydeseq2.dds": dds_module,
                "pydeseq2.ds": ds_module,
                "pydeseq2.default_inference": inference_module,
            },
        )

    def test_analyze_prefilters_genes_with_total_count_below_default(self) -> None:
        counts = pd.DataFrame(
            {
                "gene_low": [1, 2],
                "gene_edge": [5, 5],
                "gene_high": [10, 0],
            },
            index=["sample1", "sample2"],
        )
        metadata = pd.DataFrame({"condition": ["A", "B"]}, index=["sample1", "sample2"])
        captured = {}

        analyzer = Deseq2Analyzer(DESeq2Config())
        with self._patch_pydeseq2_modules(captured):
            analyzer.analyze(
                counts_df=counts,
                metadata_df=metadata,
                contrast=["condition", "B", "A"],
            )

        self.assertEqual(list(captured["counts"].columns), ["gene_edge", "gene_high"])
        self.assertEqual(list(analyzer.counts_df.columns), ["gene_edge", "gene_high"])

    def test_analyze_raises_when_prefilter_removes_all_genes(self) -> None:
        counts = pd.DataFrame(
            {"gene_low": [1, 2], "gene_lower": [0, 3]},
            index=["sample1", "sample2"],
        )
        metadata = pd.DataFrame({"condition": ["A", "B"]}, index=["sample1", "sample2"])
        captured = {}

        analyzer = Deseq2Analyzer(DESeq2Config())
        with self._patch_pydeseq2_modules(captured):
            with self.assertRaisesRegex(ValueError, "No genes remain"):
                analyzer.analyze(
                    counts_df=counts,
                    metadata_df=metadata,
                    contrast=["condition", "B", "A"],
                )

        self.assertNotIn("counts", captured)

    def test_parse_samples_from_coldata_accepts_tsv_contract(self) -> None:
        coldata_path = self.root / "coldata.tsv"
        for name in ("sample1_R1.fq.gz", "sample1_R2.fq.gz"):
            (self.root / name).write_text("stub", encoding="utf-8")
        coldata_path.write_text(
            "sample\tcondition\tr1\tr2\n"
            "sample1\tA\tsample1_R1.fq.gz\tsample1_R2.fq.gz\n",
            encoding="utf-8",
        )

        samples = cli.parse_samples_from_coldata(str(coldata_path))

        self.assertEqual(samples[0][0], "sample1")
        self.assertEqual(samples[0][1], (self.root / "sample1_R1.fq.gz").resolve())
        self.assertEqual(samples[0][2], (self.root / "sample1_R2.fq.gz").resolve())

    def test_parse_samples_from_coldata_requires_sample_column(self) -> None:
        coldata_path = self.root / "coldata.csv"
        coldata_path.write_text(
            "id,r1,r2\n"
            "sample1,sample1_R1.fq.gz,sample1_R2.fq.gz\n",
            encoding="utf-8",
        )

        with self.assertRaisesRegex(ValueError, "sample"):
            cli.parse_samples_from_coldata(str(coldata_path))

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
        tx2gene_df = pd.read_csv(outputs["tx2gene"], sep="\t")

        self.assertEqual(outputs["gene_log2_tpm"], str(self.root / "03_quant" / "gene_log2_tpm.csv"))
        log_tpm_df = pd.read_csv(outputs["gene_log2_tpm"], index_col=0)

        self.assertEqual(list(counts_df.index), ["geneA", "geneB"])
        self.assertEqual(list(counts_df.columns), ["sample1", "sample2"])
        self.assertEqual(counts_df.loc["geneA", "sample1"], 10.0)
        self.assertEqual(tpm_df.loc["geneB", "sample2"], 15.0)
        self.assertAlmostEqual(log_tpm_df.loc["geneB", "sample1"], 2.321928, places=5)
        self.assertEqual(tx2gene_df.to_dict("records"), [
            {"transcript_id": "tx1", "gene_id": "geneA"},
            {"transcript_id": "tx2", "gene_id": "geneB"},
        ])

    def test_quant_export_uses_current_command_samples_by_default(self) -> None:
        with mock.patch("rskit.cli.SalmonExpressionExporter.export_gene_tables", return_value={}) as export_gene_tables, \
             mock.patch("rskit.cli.merge_salmon_quant_tables") as merge_salmon_quant_tables:
            cli.export_quant_expression_tables(
                quant_dir=self.root / "03_quant",
                gtf_file=str(self.root / "annotation.gtf"),
                tx2gene=None,
                sample_names=["sample1"],
                merge_sf=False,
            )

        export_gene_tables.assert_called_once_with(
            salmon_dir=str(self.root / "03_quant"),
            output_dir=str(self.root / "03_quant"),
            gtf_file=str(self.root / "annotation.gtf"),
            tx2gene=None,
            sample_names=["sample1"],
        )
        merge_salmon_quant_tables.assert_not_called()

    def test_quant_export_scans_all_quant_files_when_merge_sf_requested(self) -> None:
        with mock.patch("rskit.cli.merge_salmon_quant_tables", return_value={}) as merge_salmon_quant_tables, \
             mock.patch("rskit.cli.SalmonExpressionExporter.export_gene_tables") as export_gene_tables:
            cli.export_quant_expression_tables(
                quant_dir=self.root / "03_quant",
                gtf_file=str(self.root / "annotation.gtf"),
                tx2gene=None,
                sample_names=["sample1"],
                merge_sf=True,
            )

        merge_salmon_quant_tables.assert_called_once_with(
            salmon_dir=str(self.root / "03_quant"),
            output_dir=str(self.root / "03_quant"),
            gtf_file=str(self.root / "annotation.gtf"),
            tx2gene=None,
        )
        export_gene_tables.assert_not_called()

    def test_run_deseq2_cli_prefers_quant_gene_counts_file(self) -> None:
        quant_dir = self.root / "03_quant"
        quant_dir.mkdir(parents=True, exist_ok=True)
        precomputed_counts = quant_dir / "gene_counts.csv"
        pd.DataFrame({"sample1": [10], "sample2": [12]}, index=["geneA"]).to_csv(precomputed_counts)

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
            prefilter_min_count=10,
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
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.analyze", return_value=fake_results), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.save_results", return_value={}), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.plot_volcano"), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.plot_pca"), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.plot_ma"), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.get_summary", return_value={
                 "total_genes": 1,
                 "significant_genes": 0,
                 "upregulated_genes": 0,
                 "downregulated_genes": 0,
                 "alpha": 0.05,
                 "lfc_threshold": 2.0,
             }):
            run_deseq2_cli(args)

        load_counts_from_file.assert_called_once()
        self.assertEqual(load_counts_from_file.call_args.args[0], str(precomputed_counts))
        self.assertEqual(list(load_counts_from_file.call_args.kwargs["metadata_df"].index), ["sample1", "sample2"])

    def test_run_deseq2_cli_exports_gene_tables_before_deseq(self) -> None:
        quant_dir = self.root / "03_quant"
        quant_dir.mkdir(parents=True, exist_ok=True)
        exported_counts = quant_dir / "gene_counts.csv"
        gtf_path = self._write_gtf()

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
            gtf=str(gtf_path),
            tx2gene=None,
            output_dir=str(output_dir),
            design="~condition",
            alpha=0.05,
            lfc_threshold=2.0,
            prefilter_min_count=10,
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

        with mock.patch("rskit.core.deseq2.merge_salmon_quant_tables", return_value={
                "gene_counts": str(exported_counts),
                "tx2gene": str(quant_dir / "tx2gene.tsv"),
             }) as merge_salmon_quant_tables, \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.load_counts_from_file", return_value=pd.DataFrame({"geneA": [10, 12]}, index=["sample1", "sample2"])) as load_counts_from_file, \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.analyze", return_value=fake_results), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.save_results", return_value={}), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.plot_volcano"), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.plot_pca"), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.plot_ma"), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.get_summary", return_value={
                 "total_genes": 1,
                 "significant_genes": 0,
                 "upregulated_genes": 0,
                 "downregulated_genes": 0,
                 "alpha": 0.05,
                 "lfc_threshold": 2.0,
             }):
            run_deseq2_cli(args)

        merge_salmon_quant_tables.assert_called_once_with(
            salmon_dir=str(quant_dir),
            output_dir=str(quant_dir),
            gtf_file=str(gtf_path),
            tx2gene=None,
        )
        load_counts_from_file.assert_called_once()
        self.assertEqual(load_counts_from_file.call_args.args[0], str(exported_counts))
        self.assertEqual(list(load_counts_from_file.call_args.kwargs["metadata_df"].index), ["sample1", "sample2"])

    def test_parse_contrast_validates_factor_and_levels(self) -> None:
        metadata = pd.DataFrame({"condition": ["A", "B"]}, index=["sample1", "sample2"])

        self.assertIsNone(parse_contrast(None, metadata))
        self.assertEqual(parse_contrast("condition,B,A", metadata), ["condition", "B", "A"])

        with self.assertRaisesRegex(ValueError, "factor 'batch'"):
            parse_contrast("batch,B,A", metadata)
        with self.assertRaisesRegex(ValueError, "Available levels: A, B"):
            parse_contrast("condition,C,A", metadata)
        with self.assertRaisesRegex(ValueError, "factor,level1,level2"):
            parse_contrast("condition,B", metadata)

    def test_run_deseq2_cli_validates_contrast_before_loading_counts(self) -> None:
        counts_path = self.root / "counts.csv"
        pd.DataFrame({"sample1": [10], "sample2": [12]}, index=["geneA"]).to_csv(counts_path)
        coldata_path = self.root / "coldata.csv"
        coldata_path.write_text(
            "sample,condition\n"
            "sample1,A\n"
            "sample2,B\n",
            encoding="utf-8",
        )
        args = argparse.Namespace(
            salmon_dir=None,
            gene_counts=str(counts_path),
            coldata=str(coldata_path),
            gtf=None,
            tx2gene=None,
            output_dir=str(self.root / "04_deseq2"),
            design="~condition",
            alpha=0.05,
            lfc_threshold=2.0,
            prefilter_min_count=10,
            threads=None,
            contrast="condition,C,A",
        )

        with mock.patch("rskit.core.deseq2.Deseq2Analyzer.load_counts_from_file") as load_counts:
            with self.assertRaisesRegex(ValueError, "Available levels: A, B"):
                run_deseq2_cli(args)

        load_counts.assert_not_called()

    def test_run_deseq2_cli_passes_validated_contrast_to_analysis(self) -> None:
        counts_path = self.root / "counts.csv"
        pd.DataFrame({"sample1": [10], "sample2": [12]}, index=["geneA"]).to_csv(counts_path)
        coldata_path = self.root / "coldata.csv"
        coldata_path.write_text(
            "sample,condition\n"
            "sample1,A\n"
            "sample2,B\n",
            encoding="utf-8",
        )
        args = argparse.Namespace(
            salmon_dir=None,
            gene_counts=str(counts_path),
            coldata=str(coldata_path),
            gtf=None,
            tx2gene=None,
            output_dir=str(self.root / "04_deseq2"),
            design="~condition",
            alpha=0.05,
            lfc_threshold=2.0,
            prefilter_min_count=10,
            threads=None,
            contrast="condition,B,A",
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

        with mock.patch("rskit.core.deseq2.Deseq2Analyzer.load_counts_from_file", return_value=pd.DataFrame({"geneA": [10, 12]}, index=["sample1", "sample2"])), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.analyze", return_value=fake_results) as analyze, \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.save_results", return_value={}), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.plot_volcano"), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.plot_pca"), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.plot_ma"), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.get_summary", return_value={
                 "total_genes": 1,
                 "significant_genes": 0,
                 "upregulated_genes": 0,
                 "downregulated_genes": 0,
                 "alpha": 0.05,
                 "lfc_threshold": 2.0,
             }):
            run_deseq2_cli(args)

        analyze.assert_called_once_with(contrast=["condition", "B", "A"])

    def test_run_deseq2_cli_passes_prefilter_threshold_to_config(self) -> None:
        counts_path = self.root / "counts.csv"
        pd.DataFrame({"sample1": [10], "sample2": [12]}, index=["geneA"]).to_csv(counts_path)
        coldata_path = self.root / "coldata.csv"
        coldata_path.write_text(
            "sample,condition\n"
            "sample1,A\n"
            "sample2,B\n",
            encoding="utf-8",
        )
        output_dir = self.root / "04_deseq2"
        args = argparse.Namespace(
            salmon_dir=None,
            gene_counts=str(counts_path),
            coldata=str(coldata_path),
            gtf=None,
            tx2gene=None,
            output_dir=str(output_dir),
            design="~condition",
            alpha=0.05,
            lfc_threshold=2.0,
            prefilter_min_count=25,
            threads=None,
            contrast=None,
        )
        metadata = pd.DataFrame({"condition": ["A", "B"]}, index=["sample1", "sample2"])

        with mock.patch("rskit.core.deseq2.Deseq2Analyzer") as analyzer_class:
            analyzer = analyzer_class.return_value
            analyzer.load_metadata.return_value = metadata
            analyzer.load_counts_from_file.return_value = pd.DataFrame(
                {"geneA": [10, 12]},
                index=["sample1", "sample2"],
            )
            analyzer.analyze.return_value = pd.DataFrame(
                {
                    "baseMean": [1.0],
                    "log2FoldChange": [0.5],
                    "pvalue": [0.1],
                    "padj": [0.2],
                },
                index=["geneA"],
            )
            analyzer.save_results.return_value = {}
            analyzer.get_summary.return_value = {
                "total_genes": 1,
                "significant_genes": 0,
                "upregulated_genes": 0,
                "downregulated_genes": 0,
                "alpha": 0.05,
                "lfc_threshold": 2.0,
            }

            run_deseq2_cli(args)

        config = analyzer_class.call_args.args[0]
        self.assertEqual(config.prefilter_min_count, 25)

    def test_run_deseq2_cli_writes_manifest(self) -> None:
        counts_path = self.root / "counts.csv"
        pd.DataFrame({"sample1": [10], "sample2": [12]}, index=["geneA"]).to_csv(counts_path)
        coldata_path = self.root / "coldata.csv"
        coldata_path.write_text(
            "sample,condition\n"
            "sample1,A\n"
            "sample2,B\n",
            encoding="utf-8",
        )
        output_dir = self.root / "04_deseq2"
        args = argparse.Namespace(
            salmon_dir=None,
            gene_counts=str(counts_path),
            coldata=str(coldata_path),
            gtf=None,
            tx2gene=None,
            output_dir=str(output_dir),
            design="~condition",
            alpha=0.05,
            lfc_threshold=2.0,
            prefilter_min_count=10,
            threads=None,
            contrast="condition,B,A",
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

        with mock.patch("rskit.core.deseq2.Deseq2Analyzer.load_counts_from_file", return_value=pd.DataFrame({"geneA": [10, 12]}, index=["sample1", "sample2"])), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.analyze", return_value=fake_results), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.save_results", return_value={"results": str(output_dir / "deseq2_results.csv")}), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.plot_volcano"), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.plot_pca"), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.plot_ma"), \
             mock.patch("rskit.core.deseq2.Deseq2Analyzer.get_summary", return_value={
                 "total_genes": 1,
                 "significant_genes": 0,
                 "upregulated_genes": 0,
                 "downregulated_genes": 0,
                 "alpha": 0.05,
                 "lfc_threshold": 2.0,
             }):
            run_deseq2_cli(args)

        manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["command"], "deseq2")
        self.assertEqual(manifest["samples"], ["sample1", "sample2"])
        self.assertEqual(manifest["design"], "~condition")
        self.assertEqual(manifest["contrast"], ["condition", "B", "A"])
        self.assertEqual(manifest["inputs"]["counts_file"], str(counts_path))
        self.assertEqual(manifest["outputs"]["results"], str(output_dir / "deseq2_results.csv"))
        self.assertEqual(manifest["summary"]["total_genes"], 1)

    def test_main_all_passes_quant_gene_counts_into_deseq(self) -> None:
        coldata_path = self.root / "coldata.csv"
        for name in ("sample1_R1.fq.gz", "sample1_R2.fq.gz"):
            (self.root / name).write_text("stub", encoding="utf-8")
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
            jobs=1,
            trim=False,
            force_index=False,
            skip_existing=False,
            merge_sf=True,
            design="~condition",
            contrast=None,
            alpha=0.05,
            lfc_threshold=2.0,
            prefilter_min_count=25,
        )

        exported = {"gene_counts": str(self.root / "results" / "03_quant" / "gene_counts.csv")}

        with mock.patch("rskit.cli.build_index_if_needed"), \
             mock.patch("rskit.cli.prepare_samples", return_value={"sample1": {"fq1": "a", "fq2": "b"}}), \
             mock.patch("rskit.cli.run_quantification", return_value={"sample1": {}}), \
             mock.patch("rskit.cli.merge_salmon_quant_tables", return_value=exported), \
             mock.patch("rskit.cli.os.cpu_count", return_value=2), \
             mock.patch("rskit.cli.run_deseq2_cli") as run_deseq2:
            cli.main_all(args)

        passed_args = run_deseq2.call_args.args[0]
        self.assertEqual(passed_args.gene_counts, exported["gene_counts"])
        self.assertIsNone(passed_args.salmon_dir)
        self.assertEqual(passed_args.threads, 2)
        self.assertEqual(passed_args.prefilter_min_count, 25)

    def test_main_quant_exports_current_samples_in_coldata_order(self) -> None:
        coldata_path = self.root / "coldata.csv"
        for name in ("sample1_R1.fq.gz", "sample1_R2.fq.gz", "sample2_R1.fq.gz", "sample2_R2.fq.gz"):
            (self.root / name).write_text("stub", encoding="utf-8")
        with coldata_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["sample", "r1", "r2"])
            writer.writeheader()
            writer.writerow({
                "sample": "sample1",
                "r1": str(self.root / "sample1_R1.fq.gz"),
                "r2": str(self.root / "sample1_R2.fq.gz"),
            })
            writer.writerow({
                "sample": "sample2",
                "r1": str(self.root / "sample2_R1.fq.gz"),
                "r2": str(self.root / "sample2_R2.fq.gz"),
            })

        args = argparse.Namespace(
            sample=None,
            r1=None,
            r2=None,
            coldata=str(coldata_path),
            genome_fasta=str(self.root / "genome.fa"),
            gtf_file=str(self.root / "annotation.gtf"),
            transcript_fasta=str(self.root / "transcripts.fa"),
            output_dir=str(self.root / "results"),
            index_dir=None,
            tx2gene=None,
            threads=8,
            jobs=2,
            trim=True,
            force_index=False,
            skip_existing=False,
            merge_sf=False,
        )

        with mock.patch("rskit.cli.build_index_if_needed"), \
             mock.patch("rskit.cli.prepare_samples", return_value={
                 "sample2": {"fq1": "b1", "fq2": "b2"},
                 "sample1": {"fq1": "a1", "fq2": "a2"},
             }), \
             mock.patch("rskit.cli.run_quantification", return_value={"sample1": {}, "sample2": {}}), \
             mock.patch("rskit.cli.export_quant_expression_tables", return_value={}) as export_tables:
            cli.main_quant(args)

        self.assertEqual(
            export_tables.call_args.kwargs["sample_names"],
            ["sample1", "sample2"],
        )

    def test_lfc_shrink_coefficient_matches_design_column_naming(self) -> None:
        # pydeseq2 names LFC columns after the design matrix (formulaic treatment coding)
        columns = ["Intercept", "batch[T.y]", "condition[T.B]"]
        self.assertEqual(
            _lfc_shrink_coefficient(["condition", "B", "A"], columns), "condition[T.B]"
        )
        self.assertEqual(_lfc_shrink_coefficient(["batch", "y", "x"], columns), "batch[T.y]")
        # legacy naming and continuous factors still resolve
        self.assertEqual(
            _lfc_shrink_coefficient(["condition", "B", "A"], ["Intercept", "condition_B_vs_A"]),
            "condition_B_vs_A",
        )
        self.assertEqual(
            _lfc_shrink_coefficient(["age", "old", "young"], ["Intercept", "age"]), "age"
        )
        self.assertIsNone(_lfc_shrink_coefficient(["missing", "B", "A"], ["Intercept"]))

    def test_pipeline_force_index_rebuilds_valid_index(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            pipeline = RNAseqPipeline(PipelineConfig(output_dir=tempdir))
            pipeline.aligner.align = mock.Mock(
                return_value={"bam": "b", "transcriptome_bam": "tb", "log": "l"}
            )
            pipeline.quantifier.quantify = mock.Mock(return_value={"quant": "q"})
            pipeline.indexer.build_index = mock.Mock(return_value=True)

            with mock.patch("rskit.core.pipeline.check_star_index", return_value=True):
                pipeline.run(
                    samples={"sample1": {"fq1": "r1", "fq2": "r2"}},
                    genome_fasta="genome.fa",
                    gtf_file="genes.gtf",
                    transcript_fasta="transcripts.fa",
                    index_dir=str(Path(tempdir) / "index"),
                    output_dir=tempdir,
                    quant_output_dir=str(Path(tempdir) / "03_quant"),
                    force_index=True,
                )

            pipeline.indexer.build_index.assert_called_once()
            self.assertTrue(pipeline.indexer.build_index.call_args.kwargs.get("force", False))

    def test_run_with_deseq2_uses_gene_level_counts(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            gene_counts = pd.DataFrame(
                {"sample1": [10, 20], "sample2": [30, 40]},
                index=pd.Index(["ENSG1", "ENSG2"], name="gene_id"),
            )
            pipeline = RNAseqPipeline(PipelineConfig(output_dir=tempdir))
            pipeline.run = mock.Mock(return_value={})
            pipeline.deseq2_analyzer = mock.Mock()
            pipeline.deseq2_analyzer.analyze.return_value = pd.DataFrame()

            with mock.patch.object(
                SalmonExpressionExporter,
                "build_gene_tables",
                return_value={"counts": gene_counts},
            ) as build_gene_tables:
                pipeline.run_with_deseq2(
                    samples={"sample1": {"fq1": "r1", "fq2": "r2"}, "sample2": {"fq1": "r1", "fq2": "r2"}},
                    genome_fasta="genome.fa",
                    gtf_file="genes.gtf",
                    transcript_fasta="transcripts.fa",
                    index_dir="index_dir",
                    output_dir=tempdir,
                    quant_output_dir=str(Path(tempdir) / "03_quant"),
                    metadata={"sample1": "ctrl", "sample2": "treat"},
                )

            build_gene_tables.assert_called_once()
            passed_counts = pipeline.deseq2_analyzer.analyze.call_args.args[0]
            # counts must be genes x features transposed to samples x genes with gene-level IDs
            self.assertEqual(list(passed_counts.index), ["sample1", "sample2"])
            self.assertEqual(list(passed_counts.columns), ["ENSG1", "ENSG2"])
            self.assertEqual(list(passed_counts.dtypes), ["int64", "int64"])

    def test_plot_ma_uses_fitted_stat_res(self) -> None:
        analyzer = Deseq2Analyzer(DESeq2Config())
        with self.assertRaisesRegex(ValueError, "analyze"):
            analyzer.plot_ma()

        analyzer.stat_res = mock.Mock()
        analyzer.plot_ma(save_path="out/ma.pdf")
        analyzer.stat_res.plot_MA.assert_called_once_with(save_path="out/ma.pdf")

    def test_parse_samples_from_coldata_rejects_unsafe_sample_names(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            coldata = Path(tempdir) / "coldata.csv"
            coldata.write_text(
                "sample,r1,r2\n../escape,r1.fq,r2.fq\n", encoding="utf-8"
            )

            with self.assertRaisesRegex(ValueError, "Invalid sample name"):
                cli.parse_samples_from_coldata(str(coldata))

    def test_parse_samples_from_coldata_rejects_missing_reads(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            coldata = Path(tempdir) / "coldata.csv"
            coldata.write_text(
                "sample,r1,r2\nsample1,missing_R1.fq,missing_R2.fq\n", encoding="utf-8"
            )

            with self.assertRaisesRegex(FileNotFoundError, "missing_R1.fq"):
                cli.parse_samples_from_coldata(str(coldata))

    def test_pipeline_run_skips_existing_quant_before_align(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            quant_dir = Path(tempdir) / "03_quant"
            (quant_dir / "sample1").mkdir(parents=True)
            (quant_dir / "sample1" / "quant.sf").write_text("stub", encoding="utf-8")

            pipeline = RNAseqPipeline(PipelineConfig(output_dir=tempdir))
            pipeline.aligner.align = mock.Mock()
            pipeline.quantifier.quantify = mock.Mock()

            pipeline.indexer.build_index = mock.Mock(return_value=True)

            pipeline.run(
                samples={"sample1": {"fq1": "r1", "fq2": "r2"}},
                genome_fasta="genome.fa",
                gtf_file="genes.gtf",
                transcript_fasta="transcripts.fa",
                index_dir=str(Path(tempdir) / "index"),
                output_dir=tempdir,
                quant_output_dir=str(quant_dir),
                skip_existing=True,
            )

        pipeline.aligner.align.assert_not_called()

    def test_trim_reads_reports_fastp_failure_with_stderr(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            workdirs = {
                "clean_data": root,
                "clean_data_json": root,
                "clean_data_html": root,
            }
            with mock.patch("rskit.core.base.Tool._run_command", side_effect=RuntimeError(
                "fastp failed with exit code 1: fastp\nfastp: bad input"
            )):
                with self.assertRaisesRegex(RuntimeError, "bad input"):
                    cli.trim_reads(
                        read1=root / "r1.fq",
                        read2=root / "r2.fq",
                        sample="sample1",
                        workdirs=workdirs,
                    )


if __name__ == "__main__":
    unittest.main()
