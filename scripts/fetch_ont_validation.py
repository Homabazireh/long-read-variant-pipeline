#!/usr/bin/env python3
"""Prepare the documented Clair3 HG003 Guppy 3.6.0 ONT integration example.

This legacy ONT fixture is not evidence of modern ONT chemistry performance.
Source: https://github.com/HKU-BAL/Clair3/blob/v2.0.2/docs/quick_demo/ont_quick_demo.md
"""
from __future__ import annotations

import argparse
import csv
import io
from pathlib import Path
import shutil
import subprocess
import tempfile
from fetch_public_validation import (acquisition_lock, digest, fetch, immutable_json,
                                     immutable_text, record_derived, validate_derived)


def acquire(out):
    base = 'https://www.bio8.cs.hku.hk/clair3/demo/quick_demo/ont/'
    bam = out / 'HG003_chr20_demo.bam'
    sources = [fetch(base + bam.name, bam)]
    sources.append(fetch(base + bam.name + '.bai', out / (bam.name + '.bai')))
    model_name = 'r941_prom_hac_g360+g422'
    model = out / model_name
    model.mkdir(exist_ok=True)
    model_base = 'https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/' + model_name + '/'
    checksums = {}
    for name in ('pileup.pt', 'full_alignment.pt'):
        sources.append(fetch(model_base + name, model / name))
        checksums[name] = digest(model / name)
    immutable_json(model / 'model.json', {
        'caller': 'clair3', 'caller_version': '2.0.2', 'format': 'pytorch',
        'chemistry': 'r9.4.1_450bps', 'basecaller': 'guppy_hac@v3.6.0',
        'uses_move_table': False, 'sha256': checksums, 'source': model_base,
    })
    # Publisher-documented identity is added to a new BAM; original bytes remain.
    rg_bam = out / 'HG003.ont.rg.bam'
    index = Path(str(rg_bam) + '.bai')
    marker = out / 'ont-preparation.json'
    command = ['samtools', 'addreplacerg', '-r', 'ID:HG003_ONT', '-r', 'SM:HG003',
               '-r', 'PL:ONT', '-m', 'overwrite_all', '-o', str(rg_bam), str(bam)]
    identity = {'original_sha256': digest(bam), 'command': command}
    if not validate_derived(marker, identity, [rg_bam, index]):
        with tempfile.TemporaryDirectory(prefix='.prepare-ont-', dir=out) as directory:
            temporary = Path(directory) / rg_bam.name
            temporary_command = list(command)
            temporary_command[temporary_command.index('-o') + 1] = str(temporary)
            subprocess.run(temporary_command, check=True)
            subprocess.run(['samtools', 'index', str(temporary)], check=True)
            subprocess.run(['samtools', 'quickcheck', '-v', str(temporary)], check=True)
            temporary.replace(rg_bam)
            Path(str(temporary) + '.bai').replace(index)
        record_derived(marker, identity, [rg_bam, index], read_group_command=command)
    immutable_text(out / 'evaluation.bed', 'chr20\t100000\t300000\n')
    text = io.StringIO(newline='')
    writer = csv.writer(text)
    writer.writerow(['sample', 'platform', 'reads', 'chemistry', 'basecaller', 'model'])
    writer.writerow(['HG003', 'ont', rg_bam, 'r9.4.1_450bps', 'guppy_hac@v3.6.0', model])
    immutable_text(out / 'ont.samples.csv', text.getvalue())
    immutable_json(out / 'provenance.json', {
        'purpose': 'legacy ONT bounded integration test, not modern ONT benchmark',
        'documentation': 'https://github.com/HKU-BAL/Clair3/blob/v2.0.2/docs/quick_demo/ont_quick_demo.md',
        'sources': sources, 'read_group_command': command,
        'derived_bam_sha256': digest(rg_bam),
    })
    print(out / 'ont.samples.csv')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    if args.outdir.is_symlink():
        parser.error('acquisition output directory cannot be a symbolic link')
    if not shutil.which('samtools'):
        parser.error('samtools is required for the audited read-group preparation')
    out = args.outdir.resolve()
    spec = {'sample': 'HG003', 'fixture': 'clair3-2.0.2-ont-quick-demo',
            'model': 'r941_prom_hac_g360+g422', 'region': 'chr20:100001-300000'}
    with acquisition_lock(out, 'ont-acquisition.json', spec):
        acquire(out)


if __name__ == '__main__':
    main()
