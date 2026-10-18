"""Optional assertions on real runs; set LONGREAD_TEST_OUTDIR after bin/run succeeds."""
import csv
import gzip
import os
from pathlib import Path
import unittest


@unittest.skipUnless(os.environ.get("LONGREAD_TEST_OUTDIR"), "Requires an actual completed Nextflow run, not stub outputs")
class RealLongReadOutputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.outdir = Path(os.environ["LONGREAD_TEST_OUTDIR"])
        with (cls.outdir / "calls.csv").open() as handle:
            cls.calls = list(csv.DictReader(handle))

    def test_manifest_and_indexed_calls(self):
        self.assertGreater(len(self.calls), 0)
        self.assertEqual({c["variant_type"] for c in self.calls}, {"small", "sv"})
        for call in self.calls:
            with self.subTest(call=call):
                vcf = Path(call["vcf"])
                self.assertTrue(vcf.is_absolute())
                self.assertTrue(vcf.is_file())
                self.assertTrue(Path(str(vcf) + ".tbi").is_file())
                with gzip.open(vcf, "rt") as handle:
                    columns = next(line for line in handle if line.startswith("#CHROM")).strip().split("\t")
                self.assertEqual(columns[9:], [call["sample"]])

    def test_records_use_requested_autosomes_and_sv_size_threshold(self):
        regions = {}
        with (self.outdir / "reference" / "calling-regions.bed").open() as handle:
            for line in handle:
                chrom, start, end, *_ = line.split()
                regions.setdefault(chrom, []).append((int(start), int(end)))
        for call in self.calls:
            with gzip.open(call["vcf"], "rt") as handle:
                for line in handle:
                    if line.startswith("#"):
                        continue
                    fields = line.split("\t")
                    chrom = fields[0]
                    # Boundary-spanning indels can be represented a few bases outside
                    # a supplied window. The actual truth benchmark intersects its own
                    # high-confidence regions; require autosomal contigs here.
                    self.assertIn(chrom, regions)
                    self.assertIn(chrom.removeprefix("chr"), map(str, range(1, 23)))
                    if call["variant_type"] == "sv":
                        info = dict(value.split("=", 1) for value in fields[7].split(";") if "=" in value)
                        if info.get("SVTYPE") in ("INS", "DEL"):
                            self.assertGreaterEqual(abs(int(info["SVLEN"])), 50)

    def test_phasing_and_qc_are_published(self):
        for sample in {call["sample"] for call in self.calls}:
            for suffix in ("phased.vcf.gz", "phased.vcf.gz.tbi", "haplotagged.bam", "haplotagged.bam.bai", "phasing.tsv"):
                self.assertTrue((self.outdir / "phasing" / sample / f"{sample}.{suffix}").is_file())
            self.assertTrue((self.outdir / "qc" / sample / "nanoplot" / "NanoPlot-report.html").is_file())
            self.assertTrue((self.outdir / "qc" / sample / "samtools" / f"{sample}.flagstat.txt").is_file())
            self.assertTrue((self.outdir / "qc" / sample / "mosdepth" / f"{sample}.regions.bed.gz").is_file())

    def test_real_processes_succeeded(self):
        with (self.outdir / "provenance" / "trace.tsv").open() as handle:
            trace = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(trace), 9)
        self.assertTrue(all(row["status"] in ("COMPLETED", "CACHED") and row["exit"] == "0" for row in trace))
        self.assertTrue(all(row["container"] for row in trace))


if __name__ == "__main__":
    unittest.main()
