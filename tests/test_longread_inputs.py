"""Exercise reference/sample/model safeguards without large data or external tools."""
import csv
import gzip
import hashlib
import importlib.util
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
import struct
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

SPEC = importlib.util.spec_from_file_location("longreads", Path(__file__).resolve().parents[1] / "scripts" / "run_long_reads.py")
lr = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(lr)


class LongReadValidationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.contigs = {"chr1": 20, "chr2": 30}
        self.header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\n@SQ\tSN:chr2\tLN:30\n@RG\tID:A\tSM:HG003\tPL:PACBIO\n"
        self.reads = self.root / "reads.fastq"
        self.reads.write_text("@read1\nACGT\n+\nIIII\n")

    def write_samples(self, rows):
        path = self.root / "samples.csv"
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["sample", "platform", "reads", "chemistry", "basecaller", "model"])
            writer.writerows(rows)
        return path

    def test_valid_bam_header_preserves_platform_and_sample(self):
        lr.validate_bam_header(self.header, "HG003", "hifi", self.contigs)
        bam = self.root / "header.bam"
        header = self.header.encode()
        with gzip.open(bam, "wb") as handle:
            handle.write(b"BAM\x01" + struct.pack("<i", len(header)) + header)
        self.assertEqual(lr.read_bam_header(bam), self.header)

    def test_bam_reference_mismatch_is_rejected(self):
        with self.assertRaisesRegex(lr.InputError, "reference mismatch"):
            lr.validate_bam_header(self.header.replace("LN:20", "LN:19"), "HG003", "hifi", self.contigs)

    def test_bam_missing_contig_is_rejected(self):
        with self.assertRaisesRegex(lr.InputError, "same contig"):
            lr.validate_bam_header(self.header.replace("@SQ\tSN:chr2\tLN:30\n", ""), "HG003", "hifi", self.contigs)

    def test_bam_other_sample_platform_and_sort_order_rejected(self):
        for header in (self.header.replace("SM:HG003", "SM:HG002"), self.header.replace("PL:PACBIO", "PL:ONT"), self.header.replace("SO:coordinate", "SO:queryname"), self.header.replace("@RG\tID:A\tSM:HG003\tPL:PACBIO\n", "")):
            with self.subTest(header=header), self.assertRaises(lr.InputError):
                lr.validate_bam_header(header, "HG003", "hifi", self.contigs)

    def test_fastq_checks_sequence_quality_length(self):
        lr.validate_fastq(self.reads)
        self.reads.write_text("@r\nACGT\n+\nII\n")
        with self.assertRaisesRegex(lr.InputError, "FASTQ record"):
            lr.validate_fastq(self.reads)

    def test_samples_resolve_relative_paths_and_no_implicit_models(self):
        path = self.write_samples([["HG003", "hifi", "reads.fastq", "", "", ""]])
        samples = lr.read_samples(path, self.contigs)
        self.assertEqual(samples[0]["reads"], str(self.reads))
        self.assertEqual(samples[0]["kind"], "fastq")
        self.assertIsNone(samples[0]["model"])

    def test_duplicate_sample_is_rejected(self):
        row = ["HG003", "hifi", "reads.fastq", "", "", ""]
        with self.assertRaisesRegex(lr.InputError, "Duplicate sample"):
            lr.read_samples(self.write_samples([row, row]), self.contigs)

    def test_empty_samples_and_unsupported_inputs_are_rejected(self):
        with self.assertRaisesRegex(lr.InputError, "no samples"):
            lr.read_samples(self.write_samples([]), self.contigs)
        with self.assertRaisesRegex(lr.InputError, "Invalid sample"):
            lr.read_samples(self.write_samples([["bad;command", "hifi", "reads.fastq", "", "", ""]]), self.contigs)

    def test_ont_requires_chemistry_basecaller_and_model(self):
        with self.assertRaisesRegex(lr.InputError, "ONT chemistry"):
            lr.read_samples(self.write_samples([["HG003", "ont", "reads.fastq", "", "", ""]]), self.contigs)

    def test_model_compatibility_and_checksum(self):
        model = self.root / "model"
        model.mkdir()
        checksums = {}
        for name in ("pileup.pt", "full_alignment.pt"):
            (model / name).write_bytes(b"unit-test-placeholder-not-real-weights")
            checksums[name] = lr.digest(model / name)
        chemistry = "r10.4.1_e8.2_400bps"
        basecaller = "dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
        metadata = {"caller": "clair3", "caller_version": "2.0.2", "format": "pytorch", "chemistry": chemistry, "basecaller": basecaller, "uses_move_table": False, "sha256": checksums}
        (model / "model.json").write_text(json.dumps(metadata))
        self.assertEqual(lr.validate_model(model, chemistry, basecaller), metadata)
        with self.assertRaisesRegex(lr.InputError, "basecaller"):
            lr.validate_model(model, chemistry, basecaller.replace("5.0.0", "4.3.0"))
        (model / "pileup.pt").write_bytes(b"changed")
        with self.assertRaisesRegex(lr.InputError, "SHA256 mismatch"):
            lr.validate_model(model, chemistry, basecaller)

    def test_move_table_model_is_rejected(self):
        model = self.root / "model"
        model.mkdir()
        (model / "model.json").write_text(json.dumps({"caller": "clair3", "caller_version": "2.0.2", "format": "pytorch", "chemistry": "r10.4.1_e8.2_400bps", "basecaller": "dna_r10.4.1_e8.2_400bps_sup@v5.0.0", "uses_move_table": True, "sha256": {}}))
        with self.assertRaisesRegex(lr.InputError, "uses_move_table"):
            lr.validate_model(model, "r10.4.1_e8.2_400bps", "dna_r10.4.1_e8.2_400bps_sup@v5.0.0")

    def test_regions_allow_valid_autosomes(self):
        bed = self.root / "regions.bed"
        bed.write_text("# bounded analysis\nchr1\t0\t10\nchr1\t10\t20\n")
        lr.validate_regions(bed, self.contigs)

    def test_regions_reject_overlap_outside_and_sex_chromosomes(self):
        bed = self.root / "regions.bed"
        for content in ("chr1\t0\t21\n", "chr1\t0\t10\nchr1\t9\t20\n", "chrX\t0\t10\n", "chr1\t10\t10\n", "# empty\n"):
            bed.write_text(content)
            with self.subTest(content=content), self.assertRaises(lr.InputError):
                lr.validate_regions(bed, self.contigs)

    def test_reference_digest_and_build_are_checked(self):
        fasta = self.root / "ref.fa"
        fasta.write_text(">chr1\nACGT\n")
        sha = lr.digest(fasta)
        with mock.patch.dict(lr.GRCH38_LENGTHS, {"1": 4}, clear=True):
            self.assertEqual(lr.read_reference(fasta, sha), {"chr1": 4})
        with self.assertRaisesRegex(lr.InputError, "GRCh38 length mismatch"):
            lr.read_reference(fasta, sha)
        with self.assertRaisesRegex(lr.InputError, "reference_sha256"):
            lr.read_reference(fasta, "0" * 64)

    def test_dry_run_does_not_create_output(self):
        fasta = self.root / "ref.fa"
        fasta.write_text(">chr1\nACGT\n")
        resources = self.root / "resources.json"
        resources.write_text(json.dumps({"schema_version": 1, "genome_build": "GRCh38", "reference_fasta": str(fasta), "reference_sha256": lr.digest(fasta)}))
        samples = self.write_samples([["HG003", "hifi", "reads.fastq", "", "", ""]])
        output = self.root / "new-output"
        with mock.patch.dict(lr.GRCH38_LENGTHS, {"1": 4}, clear=True), mock.patch("builtins.print"):
            self.assertEqual(lr.main(["--input", str(samples), "--resources", str(resources), "--outdir", str(output), "--profile", "apptainer", "--dry-run"]), 0)
        self.assertFalse(output.exists())

    def test_resume_rejects_changed_input(self):
        fasta = self.root / "ref.fa"
        fasta.write_text(">chr1\nACGT\n")
        resources = self.root / "resources.json"
        resources.write_text(json.dumps({"schema_version": 1, "genome_build": "GRCh38", "reference_fasta": str(fasta), "reference_sha256": lr.digest(fasta)}))
        samples = self.write_samples([["HG003", "hifi", "reads.fastq", "", "", ""]])
        output = self.root / "results"
        argv = ["--input", str(samples), "--resources", str(resources), "--outdir", str(output), "--profile", "apptainer"]
        with mock.patch.dict(lr.GRCH38_LENGTHS, {"1": 4}, clear=True), mock.patch.object(lr.shutil, "which", return_value="/bin/true"), mock.patch.object(lr, "get_nextflow_version", return_value="25.10.4"), mock.patch.object(lr, "workflow_identity", return_value={"sha256": "test"}), mock.patch.object(lr, "run_nextflow", return_value=0) as run:
            self.assertEqual(lr.main(argv), 0)
            self.assertEqual(lr.main(argv + ["--resume"]), 0)
            self.reads.write_text("@read2\nACGT\n+\nIIII\n")
            with mock.patch("builtins.print"):
                self.assertEqual(lr.main(argv + ["--resume"]), 2)
            self.assertEqual(run.call_count, 2)


class LongReadLaunchSafetyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.outdir = self.root / "output"
        self.plan = {"outdir": str(self.outdir), "workdir": str(self.outdir / "work"), "samples": [{"id": "HG003"}], "reference_sha256": "0" * 64}
        self.args = SimpleNamespace(resume=False, profile="apptainer")

    def test_unrelated_output_is_rejected_without_changes(self):
        self.outdir.mkdir()
        (self.outdir / "important.txt").write_text("existing user content")
        with self.assertRaisesRegex(lr.InputError, "nonempty"):
            lr.validate_run_destination(self.plan, False)
        self.assertEqual(list(self.outdir.iterdir()), [self.outdir / "important.txt"])

    def test_resume_requires_prior_run_without_creating_directory(self):
        with self.assertRaisesRegex(lr.InputError, "prior matching run"):
            lr.validate_run_destination(self.plan, True)
        self.assertFalse(self.outdir.exists())

    def test_concurrent_launch_is_blocked(self):
        with lr.output_lock(self.outdir):
            with self.assertRaisesRegex(lr.InputError, "Another launcher"):
                with lr.output_lock(self.outdir):
                    self.fail("Concurrent writer acquired the lock")
        with lr.output_lock(self.outdir):
            pass

    def test_version_drift_is_rejected(self):
        for version in ("25.10.3", "26.10.2", "25.10.4-edge"):
            with mock.patch.object(lr.subprocess, "run", return_value=mock.Mock(returncode=0, stdout=f"nextflow version {version}")), self.assertRaisesRegex(lr.InputError, "25.10.4 is required"):
                lr.get_nextflow_version()
        with mock.patch.object(lr.subprocess, "run", return_value=mock.Mock(returncode=0, stdout="\n      N E X T F L O W\n      version 25.10.4 build 11177\n")):
            self.assertEqual(lr.get_nextflow_version(), "25.10.4")

    def test_failure_status_is_recorded_and_lock_is_released(self):
        with mock.patch.object(lr, "workflow_identity", return_value={"sha256": "test"}), mock.patch.object(lr, "run_nextflow", return_value=17):
            self.assertEqual(lr.execute_plan(self.plan, self.args, "25.10.4"), 17)
        status = json.loads((self.outdir / "provenance/run-status.json").read_text())
        self.assertEqual(status["status"], "failed")
        self.assertEqual(status["exit_code"], 17)
        self.assertEqual(status["nextflow_version"], "25.10.4")
        self.assertIsNotNone(status["elapsed_seconds"])
        with lr.output_lock(self.outdir):
            pass

    def test_resume_preserves_previous_metrics_and_identity(self):
        with mock.patch.object(lr, "workflow_identity", return_value={"sha256": "test"}), mock.patch.object(lr, "run_nextflow", return_value=0):
            self.assertEqual(lr.execute_plan(self.plan, self.args, "25.10.4"), 0)
            trace = self.outdir / "provenance/trace.tsv"
            trace.write_text("original measured trace\n")
            original_status = (self.outdir / "provenance/run-status.json").read_text()
            self.args.resume = True
            self.assertEqual(lr.execute_plan(self.plan, self.args, "25.10.4"), 0)
        history = list((self.outdir / "provenance/history").iterdir())
        self.assertEqual(len(history), 1)
        self.assertEqual((history[0] / "trace.tsv").read_text(), "original measured trace\n")
        self.assertEqual((history[0] / "run-status.json").read_text(), original_status)
        status = json.loads((self.outdir / "provenance/run-status.json").read_text())
        self.assertTrue(status["resumed"])
        self.assertEqual(status["status"], "completed")

    def test_interrupted_status_is_written_before_output_unlock(self):
        original_write = lr.write_json
        statuses_checked = []
        def checked_write(path, value):
            if path.name == "run-status.json" and value.get("status") == "interrupted":
                with self.assertRaisesRegex(lr.InputError, "Another launcher"):
                    with lr.output_lock(self.outdir):
                        self.fail("Output lock released before final status")
                statuses_checked.append(value["exit_code"])
            original_write(path, value)
        with mock.patch.object(lr, "workflow_identity", return_value={"sha256": "test"}), mock.patch.object(lr, "run_nextflow", side_effect=KeyboardInterrupt), mock.patch.object(lr, "write_json", side_effect=checked_write):
            self.assertEqual(lr.execute_plan(self.plan, self.args, "25.10.4"), 130)
        self.assertEqual(statuses_checked, [130])
        self.assertEqual(json.loads((self.outdir / "provenance/run-status.json").read_text())["status"], "interrupted")

    def test_interruption_terminates_real_process_and_stubborn_descendant(self):
        ready = self.root / "child-pids.json"
        child_script = "import os,signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); open('child-ready','w').write(str(os.getpid())); time.sleep(30)"
        parent_script = "import json,os,subprocess,sys,time; child=subprocess.Popen([sys.executable,'-c',sys.argv[1]]); " + \
                        "time.sleep(0.1); open(sys.argv[2],'w').write(json.dumps([os.getpid(),child.pid])); time.sleep(30)"
        original_wait = subprocess.Popen.wait
        interrupted = False
        processes = []
        def interrupt_wait(process, *args, **kwargs):
            nonlocal interrupted
            if not interrupted and not args and not kwargs:
                interrupted = True
                processes.append(process)
                deadline = time.monotonic() + 5
                while (not ready.exists() or not (self.root / "child-ready").exists()) and time.monotonic() < deadline:
                    time.sleep(0.01)
                if not ready.exists() or not (self.root / "child-ready").exists():
                    raise AssertionError("Test child did not start")
                raise KeyboardInterrupt
            return original_wait(process, *args, **kwargs)
        try:
            with mock.patch.object(subprocess.Popen, "wait", new=interrupt_wait):
                with self.assertRaises(KeyboardInterrupt):
                    lr.run_nextflow([sys.executable, '-c', parent_script, child_script, str(ready)], self.root)
            parent_pid, child_pid = json.loads(ready.read_text())
            self.assertIsNotNone(processes[0].returncode)
            with self.assertRaises(ProcessLookupError):
                os.kill(parent_pid, 0)
            # An orphan killed by SIGKILL may briefly remain as a zombie until
            # init reaps it; it must never remain executing after cleanup.
            deadline = time.monotonic() + 2
            while time.monotonic() < deadline:
                stat = Path(f"/proc/{child_pid}/stat")
                if not stat.exists() or stat.read_text().split()[2] == "Z":
                    break
                time.sleep(0.01)
            else:
                self.fail("Nextflow descendant remained running after interruption")
        finally:
            for process in processes:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                original_wait(process)


if __name__ == "__main__":
    unittest.main()
