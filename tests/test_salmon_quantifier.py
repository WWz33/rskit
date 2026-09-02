import gzip
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from rskit.config import SalmonConfig
from rskit.core.salmon import SalmonExpressionExporter, SalmonQuantifier


class SalmonQuantifierTests(unittest.TestCase):
    def test_alignment_quant_command_omits_mapping_only_validate_mappings(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            transcript_fasta = root / "transcripts.fa"
            bam_file = root / "sample.bam"
            output_dir = root / "quant"
            transcript_fasta.write_text(">tx1\nACGT\n", encoding="utf-8")
            bam_file.write_text("stub", encoding="utf-8")

            quantifier = SalmonQuantifier(SalmonConfig(threads=8))

            with mock.patch("rskit.core.base.Tool._run_command", return_value=True) as run_command:
                result = quantifier.quantify(
                    transcript_fasta=str(transcript_fasta),
                    bam_file=str(bam_file),
                    output_dir=str(output_dir),
                    sample_name="sample1",
                    skip_if_exists=False,
                )

        command = run_command.call_args.args[0]
        self.assertEqual(command[:2], ["salmon", "quant"])
        self.assertIn("-t", command)
        self.assertIn("-a", command)
        self.assertIn("--seqBias", command)
        self.assertIn("--gcBias", command)
        self.assertNotIn("--posBias", command)  # experimental upstream; off by default
        self.assertNotIn("--validateMappings", command)
        self.assertEqual(result["quant"], str(output_dir / "quant.sf"))

    def test_salmon_args_enable_positional_bias(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            transcript_fasta = root / "transcripts.fa"
            bam_file = root / "sample.bam"
            transcript_fasta.write_text(">tx1\nACGT\n", encoding="utf-8")
            bam_file.write_text("stub", encoding="utf-8")

            quantifier = SalmonQuantifier(SalmonConfig(threads=8, pos_bias=True))

            with mock.patch("rskit.core.base.Tool._run_command", return_value=True) as run_command:
                quantifier.quantify(
                    transcript_fasta=str(transcript_fasta),
                    bam_file=str(bam_file),
                    output_dir=str(root / "quant"),
                    sample_name="sample1",
                    skip_if_exists=False,
                )

        self.assertIn("--posBias", run_command.call_args.args[0])

    def test_salmon_args_append_advanced_options(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            transcript_fasta = root / "transcripts.fa"
            bam_file = root / "sample.bam"
            output_dir = root / "quant"
            transcript_fasta.write_text(">tx1\nACGT\n", encoding="utf-8")
            bam_file.write_text("stub", encoding="utf-8")

            quantifier = SalmonQuantifier(
                SalmonConfig(
                    threads=8,
                    extra_args="--validateMappings --minScoreFraction 0.95",
                )
            )

            with mock.patch("rskit.core.base.Tool._run_command", return_value=True) as run_command:
                quantifier.quantify(
                    transcript_fasta=str(transcript_fasta),
                    bam_file=str(bam_file),
                    output_dir=str(output_dir),
                    sample_name="sample1",
                    skip_if_exists=False,
                )

        command = run_command.call_args.args[0]
        self.assertIn("--validateMappings", command)
        self.assertEqual(command[command.index("--minScoreFraction") + 1], "0.95")

    def test_salmon_args_replace_default_flags_without_duplicates(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            transcript_fasta = root / "transcripts.fa"
            bam_file = root / "sample.bam"
            output_dir = root / "quant"
            transcript_fasta.write_text(">tx1\nACGT\n", encoding="utf-8")
            bam_file.write_text("stub", encoding="utf-8")

            quantifier = SalmonQuantifier(
                SalmonConfig(seq_bias=True, extra_args="--seqBias")
            )

            with mock.patch("rskit.core.base.Tool._run_command", return_value=True) as run_command:
                quantifier.quantify(
                    transcript_fasta=str(transcript_fasta),
                    bam_file=str(bam_file),
                    output_dir=str(output_dir),
                    sample_name="sample1",
                    skip_if_exists=False,
                )

        command = run_command.call_args.args[0]
        self.assertEqual(command.count("--seqBias"), 1)

    def test_salmon_args_reject_protected_options(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            transcript_fasta = root / "transcripts.fa"
            bam_file = root / "sample.bam"
            output_dir = root / "quant"
            transcript_fasta.write_text(">tx1\nACGT\n", encoding="utf-8")
            bam_file.write_text("stub", encoding="utf-8")

            quantifier = SalmonQuantifier(
                SalmonConfig(extra_args="-o other_quant")
            )

            with self.assertRaisesRegex(ValueError, "-o"):
                quantifier.quantify(
                    transcript_fasta=str(transcript_fasta),
                    bam_file=str(bam_file),
                    output_dir=str(output_dir),
                    sample_name="sample1",
                    skip_if_exists=False,
                )

    def test_skip_if_exists_reruns_on_truncated_quant(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            transcript_fasta = root / "transcripts.fa"
            bam_file = root / "sample.bam"
            output_dir = root / "quant"
            transcript_fasta.write_text(">tx1\nACGT\n", encoding="utf-8")
            bam_file.write_text("stub", encoding="utf-8")
            output_dir.mkdir()
            (output_dir / "quant.sf").write_text("", encoding="utf-8")  # crashed run

            quantifier = SalmonQuantifier(SalmonConfig(threads=8))

            with mock.patch("rskit.core.base.Tool._run_command", return_value=True) as run_command:
                quantifier.quantify(
                    transcript_fasta=str(transcript_fasta),
                    bam_file=str(bam_file),
                    output_dir=str(output_dir),
                    sample_name="sample1",
                    skip_if_exists=True,
                )

        self.assertTrue(run_command.called)  # empty quant.sf must not be trusted

    def test_load_tx2gene_recovers_gene_first_column_order(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            tx2gene = Path(tempdir) / "tx2gene.tsv"
            tx2gene.write_text("ENSG1\tENST1\nENSG2\tENST2\n", encoding="utf-8")

            mapping, _ = SalmonExpressionExporter()._load_tx2gene_map(tx2gene=str(tx2gene))

        self.assertEqual(list(mapping.columns), ["transcript_id", "gene_id"])
        self.assertEqual(list(mapping["gene_id"]), ["ENSG1", "ENSG2"])

    def test_tx2gene_from_gff3_and_gzip_and_empty_fails_loud(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            gff3_content = (
                "##gff-version 3\n"
                "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene:ENSG1\n"
                "chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tID=transcript:ENST1;Parent=gene:ENSG1\n"
                "chr1\tsrc\ttranscript\t101\t200\t.\t-\t.\tID=transcript:ENST2;Parent=gene:ENSG2\n"
            )
            gff3_path = root / "annotation.gff3"
            gff3_path.write_text(gff3_content, encoding="utf-8")

            mapping = SalmonExpressionExporter()._create_tx2gene_from_gtf(str(gff3_path))
            self.assertEqual(
                list(mapping.itertuples(index=False)),
                [("ENST1", "ENSG1"), ("ENST2", "ENSG2")],
            )

            gz_path = root / "annotation.gtf.gz"
            gz_path.write_bytes(
                gzip.compress(
                    'chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST3";\n'.encode()
                )
            )
            gz_mapping = SalmonExpressionExporter()._create_tx2gene_from_gtf(str(gz_path))
            self.assertEqual(list(gz_mapping["gene_id"]), ["ENSG1"])

            empty_path = root / "genes_only.gtf"
            empty_path.write_text(
                'chr1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG1";\n', encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "No transcript-to-gene mappings"):
                SalmonExpressionExporter()._create_tx2gene_from_gtf(str(empty_path))


if __name__ == "__main__":
    unittest.main()
