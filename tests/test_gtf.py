import argparse
import tempfile
import unittest
from pathlib import Path

from rskit.utils import gtf


class GtfTests(unittest.TestCase):
    def test_open_alias_iterates_gtf_records(self) -> None:
        lines = [
            'chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id "geneA"; transcript_id "tx1";\n'
        ]

        records = list(gtf.open(lines, "ensembl"))

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].gene_id, "geneA")
        self.assertEqual(records[0].transcript_id, "tx1")

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
