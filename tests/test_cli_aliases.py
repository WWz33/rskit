import sys
import unittest
from unittest import mock

from rskit import cli


class CliAliasTests(unittest.TestCase):
    def test_quant_short_aliases_parse_to_existing_fields(self) -> None:
        argv = [
            "rskit",
            "quant",
            "-s",
            "sample1",
            "-1",
            "r1.fq.gz",
            "-2",
            "r2.fq.gz",
            "-g",
            "genome.fa",
            "-gtf",
            "annotation.gtf",
            "-gf",
            "transcripts.fa",
            "-o",
            "results",
            "-tr",
            "-fi",
            "-se",
            "-j",
            "4",
            "-ms",
        ]

        with mock.patch.object(sys, "argv", argv), \
             mock.patch("rskit.cli.main_quant") as main_quant:
            cli.main()

        args = main_quant.call_args.args[0]
        self.assertTrue(args.trim)
        self.assertTrue(args.force_index)
        self.assertTrue(args.skip_existing)
        self.assertEqual(args.jobs, 4)
        self.assertTrue(args.merge_sf)

    def test_quant_parallel_argument_is_removed(self) -> None:
        argv = [
            "rskit",
            "quant",
            "-s",
            "sample1",
            "-1",
            "r1.fq.gz",
            "-2",
            "r2.fq.gz",
            "-g",
            "genome.fa",
            "-gtf",
            "annotation.gtf",
            "-gf",
            "transcripts.fa",
            "-o",
            "results",
            "-p",
            "8",
        ]

        with mock.patch.object(sys, "argv", argv):
            with self.assertRaises(SystemExit):
                cli.main()

    def test_all_parallel_argument_is_removed(self) -> None:
        argv = [
            "rskit",
            "all",
            "-S",
            "coldata.csv",
            "-g",
            "genome.fa",
            "-gtf",
            "annotation.gtf",
            "-gf",
            "transcripts.fa",
            "-o",
            "results",
            "-p",
            "8",
        ]

        with mock.patch.object(sys, "argv", argv):
            with self.assertRaises(SystemExit):
                cli.main()

    def test_deseq2_short_aliases_parse_to_existing_fields(self) -> None:
        argv = [
            "rskit",
            "deseq2",
            "-gc",
            "counts.csv",
            "-S",
            "coldata.csv",
            "-d",
            "~batch + condition",
            "-c",
            "condition,treatment,control",
            "-a",
            "0.01",
            "-l",
            "1.5",
            "-F",
            "25",
        ]

        with mock.patch.object(sys, "argv", argv), \
             mock.patch("rskit.cli.main_deseq2") as main_deseq2:
            cli.main()

        args = main_deseq2.call_args.args[0]
        self.assertEqual(args.design, "~batch + condition")
        self.assertEqual(args.contrast, "condition,treatment,control")
        self.assertEqual(args.alpha, 0.01)
        self.assertEqual(args.lfc_threshold, 1.5)
        self.assertEqual(args.prefilter_min_count, 25)

    def test_validate_and_template_short_aliases_parse_to_existing_fields(self) -> None:
        with mock.patch.object(
            sys,
            "argv",
            ["rskit", "validate", "-S", "coldata.csv", "-d", "~condition", "-r"],
        ), mock.patch("rskit.cli.main_validate") as main_validate:
            cli.main()

        validate_args = main_validate.call_args.args[0]
        self.assertEqual(validate_args.design, "~condition")
        self.assertTrue(validate_args.check_reads)

        with mock.patch.object(
            sys,
            "argv",
            ["rskit", "template", "coldata", "-o", "coldata.csv", "-f"],
        ), mock.patch("rskit.cli.main_template") as main_template:
            cli.main()

        template_args = main_template.call_args.args[0]
        self.assertTrue(template_args.force)

    def test_wgcna_short_aliases_parse_to_existing_fields(self) -> None:
        argv = [
            "rskit",
            "wgcna",
            "-e",
            "expression.tsv",
            "-o",
            "wgcna_out",
            "-S",
            "coldata.tsv",
            "-G",
            "gene_info.tsv",
            "-sp",
            "\t",
            "-nw",
            "signed",
            "-tt",
            "signed",
            "-ms",
            "30",
            "-rs",
            "0.8",
            "-mc",
            "50",
            "-md",
            "0.25",
            "-tc",
            "2",
        ]

        with mock.patch.object(sys, "argv", argv), \
             mock.patch("rskit.cli.main_wgcna") as main_wgcna:
            cli.main()

        args = main_wgcna.call_args.args[0]
        self.assertEqual(args.sep, "\t")
        self.assertEqual(args.network_type, "signed")
        self.assertEqual(args.tom_type, "signed")
        self.assertEqual(args.min_module_size, 30)
        self.assertEqual(args.rsquared_cut, 0.8)
        self.assertEqual(args.mean_cut, 50)
        self.assertEqual(args.mediss_thresh, 0.25)
        self.assertEqual(args.tpm_cutoff, 2)


if __name__ == "__main__":
    unittest.main()
