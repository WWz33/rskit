import json
import os
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

from rskit.config import StarConfig
from rskit.core.star import StarAligner, StarIndexer, index_input_changes
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

    def test_star_index_rejects_protected_output_prefix(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            genome = root / "genome.fa"
            gtf = root / "genes.gtf"
            genome.write_text(">chr1\nACGT\n", encoding="utf-8")
            gtf.write_text("", encoding="utf-8")

            indexer = StarIndexer(StarConfig(extra_args="--outFileNamePrefix other/"))

            with self.assertRaisesRegex(ValueError, "--outFileNamePrefix"):
                indexer.build_index(str(genome), str(gtf), str(root / "index"), force=True)

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

    def test_build_index_writes_and_checks_input_fingerprint(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            genome = root / "genome.fa"
            gtf = root / "genes.gtf"
            genome.write_text(">chr1\nACGT\n", encoding="utf-8")
            gtf.write_text(
                'chr1\tsrc\texon\t1\t4\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n',
                encoding="utf-8",
            )
            index_dir = root / "index"

            indexer = StarIndexer(StarConfig(threads=2))
            with mock.patch("rskit.core.base.Tool._run_command", return_value=True), \
                 mock.patch("rskit.core.star._star_version", return_value="2.7.11b"):
                indexer.build_index(str(genome), str(gtf), str(index_dir))

            fingerprint = json.loads(
                (index_dir / ".rskit_index.json").read_text(encoding="utf-8")
            )
            self.assertEqual(fingerprint["genome_fasta"]["path"], str(genome.resolve()))
            self.assertEqual(fingerprint["genome_fasta"]["size"], genome.stat().st_size)
            self.assertEqual(fingerprint["gtf_file"]["path"], str(gtf.resolve()))
            self.assertEqual(fingerprint["star_version"], "2.7.11b")

            self.assertEqual(
                index_input_changes(str(index_dir), str(genome), str(gtf)), []
            )

            # same size, new mtime: content may have changed, must be flagged
            touched = time.time() + 100
            os.utime(genome, (touched, touched))
            self.assertEqual(
                index_input_changes(str(index_dir), str(genome), str(gtf)), ["genome_fasta"]
            )

            # an index without a fingerprint cannot be verified
            self.assertIsNone(index_input_changes(str(root / "other"), str(genome), str(gtf)))

    def test_check_star_index_rejects_empty_files(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            index_dir = Path(tempdir) / "index"
            self._write_star_index(index_dir)
            self.assertTrue(check_star_index(str(index_dir)))

            (index_dir / "SA").write_text("", encoding="utf-8")  # interrupted build

            self.assertFalse(check_star_index(str(index_dir)))

    @staticmethod
    def _write_star_index(index_dir: Path) -> None:
        # use the on-disk names real STAR produces (chrNameLength.txt, not
        # chrNameLength) so the stub cannot mask check_star_index regressions
        index_dir.mkdir(parents=True)
        for name in ("SA", "SAindex", "Genome", "chrNameLength.txt", "genomeParameters.txt"):
            (index_dir / name).write_text("stub", encoding="utf-8")


if __name__ == "__main__":
    unittest.main()
