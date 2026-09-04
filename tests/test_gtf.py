import argparse
import tempfile
import unittest
from pathlib import Path

from rskit.utils import gtf


class GtfTests(unittest.TestCase):
    def test_iter_gtf_iterates_gtf_records(self) -> None:
        lines = [
            'chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id "geneA"; transcript_id "tx1";\n'
        ]

        records = list(gtf.iter_gtf(lines, "ensembl"))

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].gene_id, "geneA")
        self.assertEqual(records[0].transcript_id, "tx1")

    def test_iter_records_tolerates_malformed_lines_and_attribute_styles(self) -> None:
        lines = [
            "\n",  # blank line must not crash
            "only\tthree\tfields\n",  # short line must be skipped
            'chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id "geneA"; transcript_id "tx1";\n',
            'chr1\tsrc\ttranscript\t101\t200\t.\t+\t.\tgene_id "geneB";transcript_id "tx2";\n',  # no space after ;
            'chr1\tsrc\ttranscript\t201\t300\t.\t+\t.\tgene_id geneC;transcript_id tx3\n',  # unquoted, no trailing ;
            'chr1\tsrc\ttranscript\t301\t400\t.\t+\t.\tID=transcript:tx4;Parent=gene:geneD\n',  # GFF3 style
            'chr1\tsrc\ttranscript\t401\t500\t.\t+\t.\tgene_id "geneE"; transcript_id "tx5"; extra\tfield "v";\n',  # tab in attrs
        ]

        records = list(gtf.iter_gtf(lines, "ensembl"))

        self.assertEqual(
            [(r.transcript_id, r.gene_id) for r in records[:3]],
            [("tx1", "geneA"), ("tx2", "geneB"), ("tx3", "geneC")],
        )
        # GFF3-style attributes parse without crashing; ID/Parent mapping is the caller's job
        self.assertEqual(records[3].meta, {"ID": "transcript:tx4", "Parent": "gene:geneD"})
        self.assertEqual((records[4].transcript_id, records[4].gene_id), ("tx5", "geneE"))

    def test_parse_line_treats_empty_numeric_fields_as_missing(self) -> None:
        rec = gtf.parse_line(
            'chr1\tsrc\tgene\t1\t100\t\t+\t\tgene_id "geneA"\n', "ensembl"
        )

        self.assertIsNone(rec.score)
        self.assertIsNone(rec.frame)
        self.assertEqual(rec.gene_id, "geneA")

    def test_gtf_tx2gene_writes_transcript_gene_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            input_path = root / "annotation.gtf"
            output_path = root / "tx2gene.csv"
            input_path.write_text(
                "\n".join(
                    [
                        'chr1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id "geneA";',
                        'chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id "geneA"; transcript_id "tx1";',
                        'chr1\tsrc\ttranscript\t101\t200\t.\t+\t.\tgene_id "geneB"; transcript_id "tx2";',
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            summary = gtf.gtf_tx2gene(
                argparse.Namespace(input=str(input_path), output=str(output_path))
            )

            self.assertEqual(summary, {"records": 3, "written": 2})
            self.assertEqual(
                output_path.read_text(encoding="utf-8"),
                "transcript_id,gene_id\n"
                "tx1,geneA\n"
                "tx2,geneB\n",
            )


if __name__ == "__main__":
    unittest.main()
