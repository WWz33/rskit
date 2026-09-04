import tempfile
import unittest
from pathlib import Path
from unittest import mock

from rskit.utils.parallel import calculate_sample_plan, process_single_sample, run_samples_parallel


class QuantSchedulingTests(unittest.TestCase):
    def test_calculate_sample_plan_caps_jobs_by_samples(self) -> None:
        plan = calculate_sample_plan(total_threads=100, requested_jobs=30, num_samples=20)

        self.assertEqual(plan.active_jobs, 20)
        self.assertEqual(plan.threads_per_sample, 5)
        self.assertEqual(plan.cap_reasons, ["sample count"])

    def test_calculate_sample_plan_caps_jobs_by_threads(self) -> None:
        plan = calculate_sample_plan(total_threads=8, requested_jobs=30, num_samples=20)

        self.assertEqual(plan.active_jobs, 8)
        self.assertEqual(plan.threads_per_sample, 1)
        self.assertEqual(plan.cap_reasons, ["total threads"])

    def test_calculate_sample_plan_rejects_non_positive_values(self) -> None:
        with self.assertRaisesRegex(ValueError, "threads"):
            calculate_sample_plan(total_threads=0, requested_jobs=1, num_samples=2)

        with self.assertRaisesRegex(ValueError, "jobs"):
            calculate_sample_plan(total_threads=8, requested_jobs=0, num_samples=2)

    def test_run_samples_parallel_uses_active_jobs_as_worker_cap(self) -> None:
        samples = {
            "sample1": {"fq1": "r1", "fq2": "r2"},
            "sample2": {"fq1": "r1", "fq2": "r2"},
            "sample3": {"fq1": "r1", "fq2": "r2"},
        }

        class ImmediateFuture:
            def __init__(self, result):
                self._result = result

            def result(self):
                return self._result

        class FakeExecutor:
            max_workers = None

            def __init__(self, max_workers):
                FakeExecutor.max_workers = max_workers

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

            def submit(self, fn, args):
                sample_name = args[0]
                return ImmediateFuture((sample_name, {"quantification": {}}))

        with mock.patch("rskit.utils.parallel.ThreadPoolExecutor", FakeExecutor), \
             mock.patch("rskit.utils.parallel.as_completed", side_effect=lambda futures: list(futures)):
            results = run_samples_parallel(
                samples=samples,
                index_dir="index",
                transcript_fasta="transcripts.fa",
                workdirs={"bam": "bam", "quant": "quant"},
                threads_per_sample=4,
                jobs=2,
            )

        self.assertEqual(FakeExecutor.max_workers, 2)
        self.assertEqual(sorted(results), ["sample1", "sample2", "sample3"])

    def test_run_samples_parallel_names_failed_sample_and_keeps_completed(self) -> None:
        samples = {
            "sample1": {"fq1": "r1", "fq2": "r2"},
            "bad": {"fq1": "r1", "fq2": "r2"},
            "sample3": {"fq1": "r1", "fq2": "r2"},
        }

        def fake_process(args):
            if args[0] == "bad":
                raise RuntimeError("STAR crashed")
            return args[0], {"quantification": {}}

        with mock.patch("rskit.utils.parallel.process_single_sample", side_effect=fake_process):
            with self.assertRaises(RuntimeError) as ctx:
                run_samples_parallel(
                    samples=samples,
                    index_dir="index",
                    transcript_fasta="transcripts.fa",
                    workdirs={"bam": "bam", "quant": "quant"},
                    threads_per_sample=4,
                    jobs=2,
                )

        self.assertIn("bad", str(ctx.exception))
        self.assertIn("sample", str(ctx.exception))  # failed sample(s) are named in the error

    def test_run_samples_parallel_handles_empty_samples(self) -> None:
        self.assertEqual(
            run_samples_parallel(
                samples={},
                index_dir="index",
                transcript_fasta="transcripts.fa",
                workdirs={"bam": "bam", "quant": "quant"},
                threads_per_sample=4,
                jobs=2,
            ),
            {},
        )

    def _sample_args(self, root: Path, skip_existing: bool = True):
        workdirs = {"bam": root / "bam", "quant": root / "quant"}
        (workdirs["bam"] / "sample1").mkdir(parents=True)
        (workdirs["quant"] / "sample1").mkdir(parents=True)
        return (
            "sample1",
            {"fq1": "missing_r1.fq", "fq2": "missing_r2.fq"},
            "index",
            "transcripts.fa",
            workdirs,
            4,
            skip_existing,
            "",
            "",
        ), workdirs

    def test_skip_existing_accepts_non_empty_quant(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            args, workdirs = self._sample_args(Path(tempdir))
            (workdirs["quant"] / "sample1" / "quant.sf").write_text("stub", encoding="utf-8")

            name, result = process_single_sample(args)

        self.assertEqual(name, "sample1")
        self.assertIn("quantification", result)

    def test_skip_existing_reruns_truncated_quant(self) -> None:
        # a crashed run leaves an empty quant.sf; it must not count as done
        with tempfile.TemporaryDirectory() as tempdir:
            args, workdirs = self._sample_args(Path(tempdir))
            (workdirs["quant"] / "sample1" / "quant.sf").write_text("", encoding="utf-8")

            # falls through to alignment, which fails on the missing reads
            with self.assertRaises(FileNotFoundError):
                process_single_sample(args)


if __name__ == "__main__":
    unittest.main()
