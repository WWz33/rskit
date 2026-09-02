import tempfile
import unittest
from pathlib import Path
from unittest import mock

from rskit.config import StarConfig
from rskit.core.star import StarAligner, StarIndexer
from rskit.utils.validators import check_star_index


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

    def test_build_index_clears_existing_directory_before_rebuild(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            index_dir = root / "index"
            self._write_star_index(index_dir)  # stale index: STAR refuses non-empty --genomeDir
            genome = root / "genome.fa"
            gtf = root / "genes.gtf"
            genome.write_text(">chr1\nACGT\n", encoding="utf-8")
            gtf.write_text("", encoding="utf-8")

            indexer = StarIndexer(StarConfig(threads=2))
            dir_state_at_run = {}

            def record_dir_state(cmd):
                dir_state_at_run["contents"] = [p.name for p in index_dir.iterdir()]
                return True

            with mock.patch("rskit.core.base.Tool._run_command", side_effect=record_dir_state):
                indexer.build_index(str(genome), str(gtf), str(index_dir), force=True)

        self.assertEqual(dir_state_at_run["contents"], [])

    def test_check_star_index_rejects_empty_files(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            index_dir = Path(tempdir) / "index"
            self._write_star_index(index_dir)
            self.assertTrue(check_star_index(str(index_dir)))

            (index_dir / "SA").write_text("", encoding="utf-8")  # interrupted build

            self.assertFalse(check_star_index(str(index_dir)))

    @staticmethod
    def _write_star_index(index_dir: Path) -> None:
        index_dir.mkdir(parents=True)
        for name in ("SA", "SAindex", "Genome", "chrNameLength", "genomeParameters.txt"):
            (index_dir / name).write_text("stub", encoding="utf-8")


if __name__ == "__main__":
    unittest.main()
