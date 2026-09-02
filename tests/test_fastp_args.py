import tempfile
import unittest
from pathlib import Path
from unittest import mock

from rskit.cli import trim_reads


class FastpArgsTests(unittest.TestCase):
    def test_fastp_args_append_trimming_options(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            read1 = root / "sample_R1.fq.gz"
            read2 = root / "sample_R2.fq.gz"
            read1.write_text("stub", encoding="utf-8")
            read2.write_text("stub", encoding="utf-8")
            workdirs = {
                "clean_data": root / "clean",
                "clean_data_json": root / "clean" / "json",
                "clean_data_html": root / "clean" / "html",
            }
            for path in workdirs.values():
                path.mkdir(parents=True, exist_ok=True)

            with mock.patch("subprocess.run") as run_command:
                trim_reads(
                    read1=read1,
                    read2=read2,
                    sample="sample1",
                    workdirs=workdirs,
                    threads=4,
                    fastp_args="--length_required 30 --cut_front",
                )

        command = run_command.call_args.args[0]
        self.assertIn("--length_required", command)
        self.assertEqual(command[command.index("--length_required") + 1], "30")
        self.assertIn("--cut_front", command)

    def test_fastp_args_reject_protected_options(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            read1 = root / "sample_R1.fq.gz"
            read2 = root / "sample_R2.fq.gz"
            read1.write_text("stub", encoding="utf-8")
            read2.write_text("stub", encoding="utf-8")
            workdirs = {
                "clean_data": root / "clean",
                "clean_data_json": root / "clean" / "json",
                "clean_data_html": root / "clean" / "html",
            }
            for path in workdirs.values():
                path.mkdir(parents=True, exist_ok=True)

            with self.assertRaisesRegex(ValueError, "--out1"):
                trim_reads(
                    read1=read1,
                    read2=read2,
                    sample="sample1",
                    workdirs=workdirs,
                    threads=4,
                    fastp_args="--out1 other_R1.fq.gz",
                )

    def test_fastp_args_reject_protected_equals_options(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            read1 = root / "sample_R1.fq.gz"
            read2 = root / "sample_R2.fq.gz"
            read1.write_text("stub", encoding="utf-8")
            read2.write_text("stub", encoding="utf-8")
            workdirs = {
                "clean_data": root / "clean",
                "clean_data_json": root / "clean" / "json",
                "clean_data_html": root / "clean" / "html",
            }
            for path in workdirs.values():
                path.mkdir(parents=True, exist_ok=True)

            with self.assertRaisesRegex(ValueError, "--out1"):
                trim_reads(
                    read1=read1,
                    read2=read2,
                    sample="sample1",
                    workdirs=workdirs,
                    threads=4,
                    fastp_args="--out1=other_R1.fq.gz",
                )

from rskit.cli import trim_reads
from rskit.cli_args import merge_extra_args


class PassthroughGuardTests(unittest.TestCase):
    def test_protected_options_block_prefix_and_attached_forms(self) -> None:
        # salmon/boost accepts unambiguous prefixes; short flags accept attached values
        for sneaky in ("--thre 4", "-p8", "--thread 4"):
            with self.assertRaisesRegex(ValueError, "Protected options"):
                merge_extra_args(
                    ["salmon", "quant", "-p", "8"], sneaky, {"-p", "--threads"}
                )

    def test_unprotected_passthrough_args_still_merge(self) -> None:
        merged = merge_extra_args(
            ["salmon", "quant", "-p", "8"], "--minScoreFraction 0.95", {"-p", "--threads"}
        )
        self.assertEqual(merged, ["salmon", "quant", "-p", "8", "--minScoreFraction", "0.95"])


if __name__ == "__main__":
    unittest.main()
