import tempfile
import unittest
from pathlib import Path
from unittest import mock

from rskit.config import SalmonConfig
from rskit.core.salmon import SalmonQuantifier


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
        self.assertIn("--posBias", command)
        self.assertNotIn("--validateMappings", command)
        self.assertEqual(result["quant"], str(output_dir / "quant.sf"))

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


if __name__ == "__main__":
    unittest.main()
