import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from rskit.utils.qc_summary import write_qc_summary


class QcSummaryTests(unittest.TestCase):
    def test_write_qc_summary_aggregates_all_three_tools(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            workdirs = {
                "clean_data": root / "01_clean_data",
                "clean_data_json": root / "01_clean_data" / "json",
                "bam": root / "02_bam",
                "quant": root / "03_quant",
                "summary": root / "00_summary",
            }
            workdirs["clean_data_json"].mkdir(parents=True)
            (workdirs["clean_data_json"] / "sample1.json").write_text(json.dumps({
                "summary": {
                    "before_filtering": {"total_reads": 10000},
                    "after_filtering": {"total_reads": 9000, "total_bases": 810000, "q30_rate": 0.95},
                }
            }), encoding="utf-8")
            bam1 = workdirs["bam"] / "sample1"
            bam1.mkdir(parents=True)
            (bam1 / "sample1_Log.final.out").write_text(
                "                      Number of input reads | 9000\n"
                "                  Uniquely mapped reads number | 7000\n"
                "                       Uniquely mapped reads % | 77.78%\n"
                " Number of reads mapped to multiple loci | 1200\n",
                encoding="utf-8",
            )
            q1 = workdirs["quant"] / "sample1"
            (q1 / "aux_info").mkdir(parents=True)
            (q1 / "aux_info" / "meta_info.json").write_text(
                json.dumps({"num_processed": 7000, "num_mapped": 6900}), encoding="utf-8"
            )
            # sample2: quant only (no trimming) -> still a row, salmon metrics only
            q2 = workdirs["quant"] / "sample2"
            (q2 / "aux_info").mkdir(parents=True)
            (q2 / "aux_info" / "meta_info.json").write_text(
                json.dumps({"num_processed": 5000, "num_mapped": 4000}), encoding="utf-8"
            )

            summary_path = write_qc_summary(workdirs)
            table = pd.read_csv(summary_path)

        self.assertEqual(list(table["sample"]), ["sample1", "sample2"])
        row = table[table["sample"] == "sample1"].iloc[0]
        self.assertEqual(int(row["input_reads"]), 10000)
        self.assertEqual(int(row["clean_reads"]), 9000)
        self.assertAlmostEqual(float(row["q30_rate"]), 0.95)
        self.assertEqual(int(row["star_unique_mapped"]), 7000)
        self.assertAlmostEqual(float(row["star_unique_rate"]), 0.7778, places=3)
        self.assertEqual(int(row["salmon_mapped_reads"]), 6900)
        self.assertAlmostEqual(float(row["salmon_mapping_rate"]), 6900 / 7000, places=4)
        row2 = table[table["sample"] == "sample2"].iloc[0]
        self.assertEqual(int(row2["salmon_mapped_reads"]), 4000)
        self.assertTrue(pd.isna(row2["input_reads"]))

    def test_write_qc_summary_returns_none_without_samples(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            (root / "03_quant").mkdir()
            self.assertIsNone(
                write_qc_summary({
                    "clean_data": root / "01_clean_data",
                    "clean_data_json": root / "01_clean_data" / "json",
                    "bam": root / "02_bam",
                    "quant": root / "03_quant",
                    "summary": root / "00_summary",
                })
            )


if __name__ == "__main__":
    unittest.main()
