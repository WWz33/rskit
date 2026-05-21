import tempfile
import unittest
from pathlib import Path
from unittest import mock

from rskit.config import StarConfig
from rskit.core.star import StarAligner


class StarArgsTests(unittest.TestCase):
    def test_star_args_replace_default_alignment_options(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            index_dir = root / "index"
            self._write_star_index(index_dir)
            r1 = root / "sample_R1.fq"
            r2 = root / "sample_R2.fq"
            r1.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
            r2.write_text("@r2\nACGT\n+\n!!!!\n", encoding="utf-8")

            aligner = StarAligner(
                StarConfig(
                    threads=4,
                    out_filter_multimap_nmax=20,
                    extra_args="--outFilterMultimapNmax 8 --alignIntronMax 500000",
                )
            )

            with mock.patch("rskit.core.base.Tool._run_command", return_value=True) as run_command:
                aligner.align(
                    index_dir=str(index_dir),
                    fq1=str(r1),
                    fq2=str(r2),
                    output_prefix=str(root / "bam" / "sample_"),
                    sample_name="sample",
                )

        command = run_command.call_args.args[0]
        self.assertEqual(command.count("--outFilterMultimapNmax"), 1)
        self.assertEqual(command[command.index("--outFilterMultimapNmax") + 1], "8")
        self.assertEqual(command[command.index("--alignIntronMax") + 1], "500000")

    def test_star_args_reject_protected_alignment_options(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            index_dir = root / "index"
            self._write_star_index(index_dir)
            r1 = root / "sample_R1.fq"
            r2 = root / "sample_R2.fq"
            r1.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
            r2.write_text("@r2\nACGT\n+\n!!!!\n", encoding="utf-8")

            aligner = StarAligner(
                StarConfig(extra_args="--outFileNamePrefix other_dir/sample_")
            )

            with self.assertRaisesRegex(ValueError, "--outFileNamePrefix"):
                aligner.align(
                    index_dir=str(index_dir),
                    fq1=str(r1),
                    fq2=str(r2),
                    output_prefix=str(root / "bam" / "sample_"),
                    sample_name="sample",
                )

    @staticmethod
    def _write_star_index(index_dir: Path) -> None:
        index_dir.mkdir(parents=True)
        for name in ("SA", "SAindex", "Genome", "genomeParameters.txt"):
            (index_dir / name).write_text("stub", encoding="utf-8")


if __name__ == "__main__":
    unittest.main()
