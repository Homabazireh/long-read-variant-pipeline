#!/usr/bin/env python3
"""Validate human long-read inputs and launch the version-pinned Nextflow workflow."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import signal
import struct
import subprocess
import sys
import time
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GRCH38_LENGTHS = dict(zip(map(str, range(1, 23)), [248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468]))
SAMPLE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$")
SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")
CLAIR3_VERSION = "2.0.2"
NEXTFLOW_VERSION = "25.10.4"


class InputError(ValueError):
    pass


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            sha.update(chunk)
    return sha.hexdigest()


def require_file(value: str, base: Path, label: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    path = path.resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise InputError(f"{label}: missing or empty file: {path}")
    # Nextflow stages with controlled basenames, but avoid unsupported shell paths.
    if any(c in str(path) for c in "\n\r\x00"):
        raise InputError(f"{label}: path contains a control character")
    return path


def read_reference(path: Path, expected_sha: str) -> dict[str, int]:
    """Hash the actual FASTA and read its contigs in one streaming pass."""
    sha = hashlib.sha256()
    contigs: dict[str, int] = {}
    name = None
    with path.open("rb") as handle:
        for raw in handle:
            sha.update(raw)
            line = raw.strip()
            if not line:
                continue
            if line.startswith(b">"):
                name = line[1:].split()[0].decode("ascii")
                if name in contigs:
                    raise InputError(f"Duplicate FASTA contig: {name}")
                if not re.fullmatch(r"[A-Za-z0-9_.:-]+", name):
                    raise InputError(f"Unsupported FASTA contig name: {name}")
                contigs[name] = 0
            elif name is None:
                raise InputError("reference_fasta must be an uncompressed FASTA")
            else:
                contigs[name] += len(line)
    if sha.hexdigest() != expected_sha.lower():
        raise InputError("reference_sha256 does not match the actual reference FASTA")
    autosomes = {k: v for k, v in contigs.items() if k.removeprefix("chr") in GRCH38_LENGTHS}
    if not autosomes:
        raise InputError("Reference has no GRCh38 autosomes (1..22 or chr1..chr22)")
    if any("chr" + x in autosomes and x in autosomes for x in GRCH38_LENGTHS):
        raise InputError("Reference mixes duplicate chr-prefixed and unprefixed autosomes")
    for name, length in autosomes.items():
        if length != GRCH38_LENGTHS[name.removeprefix("chr")]:
            raise InputError(f"GRCh38 length mismatch for {name}: {length}")
    return contigs


def read_bam_header(path: Path) -> str:
    try:
        with gzip.open(path, "rb") as handle:
            if handle.read(4) != b"BAM\x01":
                raise InputError(f"Not a BAM file: {path}")
            size_raw = handle.read(4)
            if len(size_raw) != 4:
                raise InputError(f"Truncated BAM header: {path}")
            size = struct.unpack("<i", size_raw)[0]
            if not 0 < size <= 64 * 1024 * 1024:
                raise InputError(f"Invalid or unsupported BAM header size: {path}")
            header = handle.read(size)
            if len(header) != size:
                raise InputError(f"Truncated BAM header: {path}")
            return header.decode("utf-8")
    except (OSError, EOFError, UnicodeDecodeError) as exc:
        raise InputError(f"Cannot read BAM header: {path}: {exc}") from exc


def validate_bam_header(header: str, sample: str, platform: str, contigs: dict[str, int]) -> None:
    rows: dict[str, list[dict[str, str]]] = {}
    for line in header.splitlines():
        fields = line.split("\t")
        if not fields:
            continue
        tags = dict(x.split(":", 1) for x in fields[1:] if ":" in x)
        rows.setdefault(fields[0], []).append(tags)
    if len(rows.get("@HD", [])) != 1 or rows["@HD"][0].get("SO") != "coordinate":
        raise InputError(f"{sample}: BAM must declare coordinate sort order; unsorted/unmapped BAM is unsupported")
    sq = rows.get("@SQ", [])
    if not sq:
        raise InputError(f"{sample}: BAM must contain aligned-reference @SQ records")
    seen = set()
    for record in sq:
        name = record.get("SN", "")
        if name in seen:
            raise InputError(f"{sample}: duplicate BAM @SQ {name}")
        seen.add(name)
        if name not in contigs or record.get("LN") != str(contigs[name]):
            raise InputError(f"{sample}: BAM reference mismatch for {name}")
    # Same reference dictionary is required; avoiding accidental subset/reference mixing.
    if seen != set(contigs):
        raise InputError(f"{sample}: BAM and reference must have the same contig dictionary")
    groups = rows.get("@RG", [])
    if not groups:
        raise InputError(f"{sample}: BAM requires read groups with SM and PL")
    expected_pl = "PACBIO" if platform == "hifi" else "ONT"
    ids = set()
    for group in groups:
        if not group.get("ID") or group["ID"] in ids:
            raise InputError(f"{sample}: missing or duplicated BAM read-group ID")
        ids.add(group["ID"])
        if group.get("SM") != sample or group.get("PL", "").upper() != expected_pl:
            raise InputError(f"{sample}: BAM read groups must have SM={sample} and PL={expected_pl}")


def validate_fastq(path: Path) -> None:
    opener = gzip.open if path.name.endswith(".gz") else open
    try:
        with opener(path, "rt") as handle:
            name, bases, plus, qualities = [handle.readline().rstrip("\r\n") for _ in range(4)]
        if not name.startswith("@") or not plus.startswith("+") or not bases or len(bases) != len(qualities):
            raise InputError(f"Invalid first FASTQ record: {path}")
    except (OSError, EOFError, UnicodeDecodeError) as exc:
        raise InputError(f"Cannot read FASTQ: {path}: {exc}") from exc


def validate_model(path: Path, chemistry: str, basecaller: str) -> dict:
    metadata = require_file(str(path / "model.json"), path, "Clair3 model metadata")
    try:
        model = json.loads(metadata.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise InputError(f"Invalid Clair3 model.json: {exc}") from exc
    if not isinstance(model, dict) or not isinstance(model.get("sha256"), dict):
        raise InputError("Clair3 model.json must be an object with a sha256 object")
    expected = {"caller": "clair3", "caller_version": CLAIR3_VERSION, "format": "pytorch", "chemistry": chemistry, "basecaller": basecaller, "uses_move_table": False}
    for key, value in expected.items():
        if model.get(key) != value:
            raise InputError(f"Clair3 model.json {key} must equal {value!r}")
    modern = chemistry == "r10.4.1_e8.2_400bps" and re.fullmatch(r"dna_r10\.4\.1_e8\.2_400bps_(sup|hac)@v[0-9]+\.[0-9]+\.[0-9]+", basecaller)
    legacy = chemistry == "r9.4.1_450bps" and re.fullmatch(r"guppy_(sup|hac)@v[0-9]+\.[0-9]+\.[0-9]+", basecaller)
    if not (modern or legacy):
        raise InputError("Declare supported ONT chemistry and exact basecaller: R10.4.1/E8.2/400bps Dorado or R9.4.1/450bps Guppy SUP/HAC")
    for name in ("pileup.pt", "full_alignment.pt"):
        weights = require_file(str(path / name), path, "Clair3 PyTorch checkpoint")
        expected_sha = model.get("sha256", {}).get(name, "")
        if not SHA256.fullmatch(expected_sha) or digest(weights) != expected_sha.lower():
            raise InputError(f"Clair3 checkpoint SHA256 mismatch or missing: {name}")
    return model


def validate_regions(path: Path, contigs: dict[str, int]) -> None:
    previous: dict[str, int] = {}
    count = 0
    with path.open() as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise InputError(f"BED line {lineno}: expected at least three columns")
            chrom = fields[0]
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise InputError(f"BED line {lineno}: coordinates must be integers") from exc
            if chrom not in contigs or chrom.removeprefix("chr") not in GRCH38_LENGTHS:
                raise InputError(f"BED line {lineno}: unsupported or non-autosomal contig {chrom}")
            if not 0 <= start < end <= contigs[chrom]:
                raise InputError(f"BED line {lineno}: interval is outside the reference")
            if start < previous.get(chrom, 0):
                raise InputError(f"BED line {lineno}: intervals must be sorted and nonoverlapping within each contig")
            previous[chrom] = end
            count += 1
    if not count:
        raise InputError("Regions BED contains no intervals")


def read_samples(path: Path, contigs: dict[str, int]) -> list[dict]:
    samples = []
    seen = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        expected = {"sample", "platform", "reads", "chemistry", "basecaller", "model"}
        if not reader.fieldnames or not expected.issubset(reader.fieldnames):
            raise InputError(f"Samples CSV requires columns: {','.join(sorted(expected))}")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise InputError("Malformed samples CSV row")
            row = {key: value.strip() for key, value in row.items()}
            sample, platform = row["sample"], row["platform"]
            if not SAMPLE.fullmatch(sample):
                raise InputError(f"Invalid sample ID: {sample!r}")
            if sample in seen:
                raise InputError(f"Duplicate sample: {sample}")
            seen.add(sample)
            if platform not in ("hifi", "ont"):
                raise InputError(f"{sample}: platform must be hifi or ont")
            reads = require_file(row["reads"], path.parent, f"{sample} reads")
            if reads.name.endswith(".bam"):
                kind = "bam"
                validate_bam_header(read_bam_header(reads), sample, platform, contigs)
            elif reads.name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
                kind = "fastq"
                validate_fastq(reads)
            else:
                raise InputError(f"{sample}: reads must be basecalled FASTQ[.gz] or aligned BAM")
            model_path, model_metadata = None, None
            if platform == "ont":
                if not all(row[x] for x in ("chemistry", "basecaller", "model")):
                    raise InputError(f"{sample}: ONT chemistry, basecaller and model directory are required")
                model_path = Path(row["model"]).expanduser()
                if not model_path.is_absolute():
                    model_path = path.parent / model_path
                model_path = model_path.resolve()
                model_metadata = validate_model(model_path, row["chemistry"], row["basecaller"])
            elif row["model"]:
                raise InputError(f"{sample}: HiFi uses DeepVariant PACBIO bundled weights; leave model blank")
            samples.append({"id": sample, "platform": platform, "reads": str(reads), "kind": kind, "chemistry": row["chemistry"], "basecaller": row["basecaller"], "model": str(model_path) if model_path else None, "model_metadata": model_metadata})
    if not samples:
        raise InputError("Samples CSV contains no samples")
    return samples


def build_plan(args: argparse.Namespace) -> dict:
    sample_file = require_file(args.input, Path.cwd(), "samples CSV")
    resources_file = require_file(args.resources, Path.cwd(), "resources JSON")
    try:
        resources = json.loads(resources_file.read_text())
    except json.JSONDecodeError as exc:
        raise InputError(f"Invalid resources JSON: {exc}") from exc
    if not isinstance(resources, dict):
        raise InputError("Resources JSON must be an object")
    if resources.get("schema_version") != 1 or resources.get("genome_build") != "GRCh38":
        raise InputError("Resources require schema_version=1 and genome_build=GRCh38")
    fasta_value = resources.get("reference_fasta", "")
    if not fasta_value or not Path(fasta_value).is_absolute():
        raise InputError("reference_fasta must be an absolute path")
    fasta = require_file(fasta_value, resources_file.parent, "reference FASTA")
    sha = resources.get("reference_sha256", "")
    if not isinstance(sha, str) or not SHA256.fullmatch(sha):
        raise InputError("reference_sha256 must be a 64-character SHA256 digest")
    contigs = read_reference(fasta, sha)
    if args.cpus < 2 or args.memory_gb < 4:
        raise InputError("At least --cpus 2 and --memory-gb 4 are required; whole genomes need substantially more")
    regions = None
    if args.regions:
        regions = require_file(args.regions, Path.cwd(), "regions BED")
        validate_regions(regions, contigs)
    samples = read_samples(sample_file, contigs)
    outdir = Path(args.outdir).expanduser().resolve()
    workdir = Path(args.work_dir).expanduser().resolve() if args.work_dir else outdir / "work"
    return {"schema_version": 1, "samples": samples, "reference_fasta": str(fasta), "reference_sha256": sha.lower(), "reference_contigs": contigs, "genome_build": "GRCh38", "regions": str(regions) if regions else None, "outdir": str(outdir), "workdir": str(workdir), "cpus": args.cpus, "memory_gb": args.memory_gb, "executor_name": args.executor, "queue": args.queue or "", "profile": args.profile, "samples_sha256": digest(sample_file), "resources_sha256": digest(resources_file), "regions_sha256": digest(regions) if regions else None, "read_provenance": [{"sample": s["id"], "path": s["reads"], "size_bytes": Path(s["reads"]).stat().st_size, "mtime_ns": Path(s["reads"]).stat().st_mtime_ns} for s in samples]}


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--input", required=True, help="Samples CSV")
    result.add_argument("--resources", required=True, help="Reference/resource JSON")
    result.add_argument("--outdir", required=True)
    result.add_argument("--profile", choices=("apptainer", "docker"), required=True)
    result.add_argument("--executor", choices=("local", "slurm"), default="local")
    result.add_argument("--queue", help="Optional Slurm partition")
    result.add_argument("--cpus", type=int, default=8, help="Per-task maximum; local scheduling also respects this total")
    result.add_argument("--memory-gb", type=int, default=32)
    result.add_argument("--regions", help="Optional 0-based half-open autosomal BED, for scoped analysis")
    result.add_argument("--work-dir", help="Persistent Nextflow work directory (default OUTDIR/work)")
    result.add_argument("--resume", action="store_true")
    result.add_argument("--dry-run", action="store_true", help="Validate inputs and print plan without executing or writing")
    return result


def validate_run_destination(plan: dict, resume: bool) -> None:
    outdir = Path(plan["outdir"])
    params_file = outdir / "provenance" / "run-inputs.json"
    if outdir.exists() and not outdir.is_dir():
        raise InputError("--outdir is not a directory")
    if not params_file.is_file():
        if resume:
            raise InputError("--resume requires an existing run-inputs.json from a prior matching run")
        if outdir.exists() and any(p.name != ".run.lock" for p in outdir.iterdir()):
            raise InputError("Output directory is nonempty and does not contain a matching pipeline run")
        return
    if not resume:
        raise InputError("Output directory already contains a run; use --resume or a new --outdir")
    try:
        previous = json.loads(params_file.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise InputError(f"Cannot read previous run manifest: {exc}") from exc
    keys = ("samples", "reference_fasta", "reference_sha256", "regions", "regions_sha256", "read_provenance", "samples_sha256", "resources_sha256")
    if not isinstance(previous, dict) or any(previous.get(key) != plan.get(key) for key in keys):
        raise InputError("--resume inputs differ from the existing run; choose a new --outdir")


@contextmanager
def output_lock(outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / ".run.lock").open("a+") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise InputError(f"Another launcher is writing to {outdir}") from exc
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def get_nextflow_version() -> str:
    result = subprocess.run(["nextflow", "-version"], text=True, capture_output=True, check=False)
    match = re.search(r"(?m)^\s*(?:nextflow\s+)?version\s+([0-9]+\.[0-9]+\.[0-9]+(?:-[a-z0-9]+)?)", result.stdout, re.IGNORECASE)
    if result.returncode != 0 or not match:
        raise InputError("Could not determine the installed Nextflow version")
    version = match.group(1)
    if version != NEXTFLOW_VERSION:
        raise InputError(f"Nextflow {NEXTFLOW_VERSION} is required; found {version}")
    return version


def workflow_identity() -> dict:
    paths = [ROOT / "main.nf", ROOT / "nextflow.config", ROOT / "scripts" / "run_long_reads.py", ROOT / "bin" / "run"]
    paths += sorted((ROOT / "modules" / "longread").glob("*.nf"))
    paths += sorted((ROOT / "configs").glob("longread*.config"))
    files = {str(p.relative_to(ROOT)): digest(p) for p in paths}
    identity = {"sha256": hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest(), "files": files, "git_commit": None, "git_dirty": None}
    if shutil.which("git"):
        revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=False)
        if revision.returncode == 0:
            identity["git_commit"] = revision.stdout.strip()
            status = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, capture_output=True, text=True, check=False)
            identity["git_dirty"] = bool(status.stdout.strip())
    return identity


def write_json(path: Path, value) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def terminate_process_group(process) -> None:
    """Stop Nextflow and any remaining local descendants before unlocking outputs."""
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        pass
    # A child can outlive a parent that exits promptly after SIGTERM. Always
    # stop any remaining members, even if the group leader was already reaped.
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait()


def run_nextflow(command: list[str], cwd: Path) -> int:
    process = subprocess.Popen(command, cwd=cwd, start_new_session=True)
    try:
        code = process.wait()
        if code:
            terminate_process_group(process)
        return code
    except BaseException:
        terminate_process_group(process)
        raise


def execute_plan(plan: dict, args: argparse.Namespace, version: str) -> int:
    outdir = Path(plan["outdir"])
    validate_run_destination(plan, args.resume)
    with output_lock(outdir):
        # Recheck under the lock to prevent a second launcher racing the first.
        validate_run_destination(plan, args.resume)
        provenance = outdir / "provenance"
        provenance.mkdir(parents=True, exist_ok=True)
        attempt_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:8]
        if args.resume:
            archive = provenance / "history" / attempt_id
            archive.mkdir(parents=True, exist_ok=True)
            for name in ("run-inputs.json", "run-status.json", "command.json", "nextflow.log", "trace.tsv", "execution-report.html", "timeline.html", "workflow-dag.html"):
                source = provenance / name
                if source.is_file():
                    shutil.copy2(source, archive / name)
        params_file = provenance / "run-inputs.json"
        write_json(params_file, plan)
        command = ["nextflow", "-log", str(provenance / "nextflow.log"), "run", str(ROOT / "main.nf"), "-profile", args.profile, "-params-file", str(params_file), "-work-dir", plan["workdir"]]
        if args.resume:
            command += ["-resume"]
        write_json(provenance / "command.json", command)
        started = time.monotonic()
        status = {"attempt_id": attempt_id, "status": "running", "started_at": datetime.now(timezone.utc).isoformat(), "finished_at": None, "elapsed_seconds": None, "exit_code": None, "resumed": args.resume, "nextflow_version": version, "workflow": workflow_identity()}
        write_json(provenance / "run-status.json", status)
        exit_code = 1
        try:
            exit_code = run_nextflow(command, provenance)
            return exit_code
        except KeyboardInterrupt:
            exit_code = 130
            return exit_code
        finally:
            status.update({"status": "completed" if exit_code == 0 else "interrupted" if exit_code == 130 else "failed", "finished_at": datetime.now(timezone.utc).isoformat(), "elapsed_seconds": round(time.monotonic() - started, 3), "exit_code": exit_code})
            write_json(provenance / "run-status.json", status)


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    try:
        plan = build_plan(args)
        validate_run_destination(plan, args.resume)
        if args.dry_run:
            print(json.dumps(plan, indent=2))
            return 0
        required = ["nextflow", args.profile] + (["sbatch"] if args.executor == "slurm" else [])
        for executable in required:
            if shutil.which(executable) is None:
                raise InputError(f"Required executable is unavailable: {executable}")
        return execute_plan(plan, args, get_nextflow_version())
    except (InputError, OSError, TypeError, KeyError, UnicodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    def interrupt_on_termination(signum, frame):
        raise KeyboardInterrupt("Termination requested")

    signal.signal(signal.SIGTERM, interrupt_on_termination)
    raise SystemExit(main())
