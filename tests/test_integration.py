"""End-to-end test against the real STAR and fastp binaries.

The unit tests mock every subprocess; this module exercises the actual
command wiring (index build -> fingerprint -> fastp trim -> STAR align with
gzipped mates) on a tiny genome. It is skipped automatically when STAR or
fastp is not on PATH, so CI (which installs neither) is unaffected.
"""

import gzip
import random
import shutil
import tempfile
import unittest
from pathlib import Path

from rskit import cli
from rskit.config import StarConfig
from rskit.core.star import StarAligner, StarIndexer, index_input_changes
from rskit.utils.validators import check_star_index

# a seeded, non-repetitive genome: a periodic ACGT sequence would make every
# read a multimapper
_RNG = random.Random(42)
GENOME = "".join(_RNG.choice("ACGT") for _ in range(2000))
GTF_LINES = [
    'chr1\tsrc\tgene\t1\t2000\t.\t+\t.\tgene_id "g1";',
    'chr1\tsrc\ttranscript\t1\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
    'chr1\tsrc\texon\t1\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
]

_COMPLEMENT = str.maketrans("ACGT", "TGCA")

# tiny genomes also need a small SA index; real genomes use the default
SMALL_GENOME_STAR_ARGS = "--genomeSAindexNbases 2"


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _write_paired_reads(reads_dir: Path, gzipped: bool) -> tuple:
    """Six FR read pairs (50 bp, Phred 40) cut out of the genome."""
    r1_records, r2_records = [], []
    for i, start in enumerate(range(100, 1600, 300)):
        r1_records += [f"@read{i}_1", GENOME[start:start + 50], "+", "I" * 50]
        r2_records += [f"@read{i}_2", _revcomp(GENOME[start + 100:start + 150]), "+", "I" * 50]
    text = "\n".join(r1_records) + "\n", "\n".join(r2_records) + "\n"

    suffix = ".fq.gz" if gzipped else ".fq"
    paths = []
    for name, content in zip(("sample1_R1", "sample1_R2"), text):
        path = reads_dir / (name + suffix)
        if gzipped:
            path.write_bytes(gzip.compress(content.encode("utf-8")))
        else:
            path.write_text(content, encoding="utf-8")
        paths.append(path)
    return tuple(paths)


@unittest.skipUnless(shutil.which("STAR") and shutil.which("fastp"),
                     "STAR and fastp must be on PATH for the integration test")
class StarFastpIntegration(unittest.TestCase):
    def test_index_fingerprint_trim_and_align_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            genome = root / "genome.fa"
            genome.write_text(">chr1\n" + GENOME + "\n", encoding="utf-8")
            gtf = root / "genes.gtf"
            gtf.write_text("\n".join(GTF_LINES) + "\n", encoding="utf-8")
            index_dir = root / "00_index"

            indexer = StarIndexer(
                StarConfig(threads=1, sjdb_overhang=49, extra_args=SMALL_GENOME_STAR_ARGS)
            )
            indexer.build_index(str(genome), str(gtf), str(index_dir))

            self.assertTrue(check_star_index(str(index_dir)))
            self.assertEqual(
                index_input_changes(str(index_dir), str(genome), str(gtf)), []
            )

            # an edited annotation must be flagged against the fingerprint
            gtf.write_text("\n".join(GTF_LINES) + "\n" + GTF_LINES[2] + "\n", encoding="utf-8")
            self.assertEqual(
                index_input_changes(str(index_dir), str(genome), str(gtf)), ["gtf_file"]
            )
            gtf.write_text("\n".join(GTF_LINES) + "\n", encoding="utf-8")

            # an index without a fingerprint cannot be verified
            empty_index = root / "handmade_index"
            empty_index.mkdir()
            self.assertIsNone(
                index_input_changes(str(empty_index), str(genome), str(gtf))
            )

            reads_dir = root / "reads"
            reads_dir.mkdir()
            r1, r2 = _write_paired_reads(reads_dir, gzipped=True)

            workdirs = {
                "clean_data": root / "01_clean_data",
                "clean_data_json": root / "01_clean_data" / "json",
                "clean_data_html": root / "01_clean_data" / "html",
            }
            for path in workdirs.values():
                path.mkdir(parents=True)

            r1_clean, r2_clean = cli.trim_reads(r1, r2, "sample1", workdirs, threads=1)
            self.assertTrue(r1_clean.endswith(".fq.gz"))
            self.assertTrue(Path(r1_clean).stat().st_size > 0)
            self.assertTrue(Path(r2_clean).stat().st_size > 0)
            self.assertTrue((workdirs["clean_data_json"] / "sample1.json").exists())

            aligner = StarAligner(StarConfig(threads=1, sjdb_overhang=49))
            results = aligner.align(
                str(index_dir), r1_clean, r2_clean,
                str(root / "02_bam" / "sample1_"), sample_name="sample1",
            )

            transcriptome_bam = Path(results["transcriptome_bam"])
            self.assertTrue(transcriptome_bam.exists())
            self.assertTrue(transcriptome_bam.stat().st_size > 0)

            # reads were cut verbatim from the genome: they must actually map
            unique_mapped = 0
            for line in Path(results["log"]).read_text(encoding="utf-8").splitlines():
                if "Uniquely mapped reads number" in line:
                    unique_mapped = int(line.partition("|")[2].strip())
            self.assertGreater(unique_mapped, 0)


if __name__ == "__main__":
    unittest.main()
