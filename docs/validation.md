# Measured validation

Validation date: **2026-10-18**. This initial implementation passed real public
HG003 regional integration runs on Linux x86-64 with Apptainer, Nextflow 25.10.4,
Java 21, CPU execution, and the full GRCh38 no-alt analysis reference. The host
has an AMD EPYC 9454P and approximately 251 GiB RAM. Calling runs requested
8 CPUs and 32 GiB; annotation and benchmarking requested 4 CPUs and 16 GiB.

These measurements establish that the tested workflow paths execute and produce
assessable results. They do not establish whole-genome performance, clinical
validity, or a technology ranking. Model-training overlap with HG003 was not
audited. No DRAGEN installation or DRAGEN callset was benchmarked.

## Executed workflows

| Input | Completed real stages | Original small-variant records | SV records |
| --- | --- | ---: | ---: |
| PacBio HiFi BAM | Input validation, QC, coverage, DeepVariant, Sniffles2, phasing, haplotagging | 603 | 2 |
| PacBio HiFi reconstructed FASTQ | Whole-reference minimap2 alignment and all downstream stages | 603 | 2 |
| Legacy ONT BAM | Input validation, QC, coverage, Clair3, Sniffles2, phasing, haplotagging | 553 | 3 |

Each run completed all nine Nextflow tasks. HiFi and ONT resumes reused all nine
cached tasks. The accepted HiFi BAM and its haplotagged output both contain 807
alignment records; haplotagging also preserved all 807 records with an empty
small-variant input in a separate real-tool test. The FASTQ test used 804 primary
reads reconstructed from that public BAM; it is not validation of raw sequencing
acquisition or upstream basecalling.

The [per-task execution measurements](validation/long-read-execution.csv)
retain original non-cached timings and peak resident memory. For example, HiFi
DeepVariant took 36.5 seconds of task runtime with a sampled 6.4 GB peak RSS;
these regional timings exclude acquisition and are not whole-genome estimates.
[Image artifact hashes](validation/long-read-images.json) identify the exact
Apptainer SIF files used alongside the software version pins.

The HiFi source has 23 additional unused decoy dictionary entries and a published
sample alias (`HG3-bc20`). Fixture preparation explicitly checks every read and
mate contig, re-encodes through SAM to update BAM reference IDs safely, removes
only unused extra dictionary entries, and standardizes the documented sample
to HG003. Original files and transformation provenance are retained. This
fixture helper is not applied automatically to study data.

The ONT fixture is the official Clair3 HG003 quick demonstration using
`r9.4.1_450bps`, Guppy HAC 3.6.0, and the publisher's converted PyTorch
`r941_prom_hac_g360+g422` model. The original BAM lacks read groups; an explicitly
recorded step adds the publisher-documented sample/platform to a separate BAM.
Modern R10.4.1/Dorado configurations are supported by the input contract but
have not been executed in this validation.

## Small-variant benchmarks

Both evaluations use **GIAB HG003 GRCh38 v4.2.1**, hap.py 0.3.12 with RTG vcfeval,
the complete truth VCF, and truth confidence regions intersected with the
specified evaluation BED. PASS results are shown below. The HiFi interval is
`chr20:10000001-10200000` (194,051 input confidence BED bases); the ONT interval
is `chr20:100001-300000` (184,613 input confidence BED bases). hap.py's standard
insertion-boundary adjustment was enabled; its native adjusted confidence
sizes are 194,597 bases for HiFi and 185,165 for ONT. Both input and adjusted
sizes are recorded. Intervals are 1-based inclusive here; BED files are
0-based half-open.

**The inputs and intervals differ, so the rows cannot support a platform ranking.**

| Input | Type | Truth TP | Query TP | FN | FP | Precision | Recall | F1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| HiFi / DeepVariant | SNP | 388 | 390 | 1 | 0 | 1.000000 | 0.997429 | 0.998713 |
| HiFi / DeepVariant | INDEL | 62 | 67 | 0 | 0 | 1.000000 | 1.000000 | 1.000000 |
| Legacy ONT / Clair3 | SNP | 402 | 403 | 0 | 1 | 0.997525 | 1.000000 | 0.998761 |
| Legacy ONT / Clair3 | INDEL | 44 | 44 | 15 | 5 | 0.897959 | 0.745763 | 0.814815 |

Truth and query TP counts can differ because equivalent haplotypes may have
different variant representations. Precision uses query TP; recall uses truth
TP. One ONT indel genotype error was reported. Calls outside assessed regions
are counted separately, not treated as false positives. Native hap.py results
were retained; the initial report was regenerated after correcting extraction
of query TP from hap.py's extended aggregate table. That correction and its
script hash are recorded in the exported provenance.

See [machine-readable metrics](validation/long-read-metrics.csv) and
[benchmark provenance](validation/long-read-provenance.json) for full-precision
numbers, resource hashes, matching policy, container digests, and scope.

## Annotation and SV benchmarking

Real offline VEP 116 annotation completed for all four original HiFi/ONT small
and structural callsets using the complete Ensembl 116 GRCh38 cache. The
ClinVar 2026-10-11 and gnomAD 4.1 genomes resources were explicitly restricted to
the test intervals. Resource checksums, allele matches, HGVS, MANE fields, and
all transcript consequences were checked. These regional resource files must
not be reused as whole-genome annotation databases.

| Callset | Records after small-variant normalization | Consequence rows |
| --- | ---: | ---: |
| HiFi DeepVariant | 618 | 13,272 |
| HiFi Sniffles2 | 2 | 44 |
| ONT Clair3 | 561 | 1,010 |
| ONT Sniffles2 | 3 | 12 |

Normalization can split multiallelic records, so counts differ from original
calls. ClinVar/gnomAD match counts in the
[annotation evidence](validation/long-read-annotation.json) count consequence
rows, not unique variants. SV annotations explicitly mark small-variant
resources as not applicable.

The actual Truvari workflow passed deliberately perturbed synthetic tests
(2 SV TP, 1 FN, 1 FP, and 1 genotype error), empty-query tests, and interval/size
bin checks. These test report semantics; **no real GIAB SV accuracy result is
claimed**. GIAB SV evaluation requires a matching release, its documented
assessment regions, and an explicit checksum-bound matching policy.

## Reproduction and remaining validation

Install prerequisites with the README instructions and fetch bounded fixtures:

```bash
python3 scripts/fetch_public_validation.py --outdir /analysis/hg003 --platforms hifi
python3 scripts/fetch_ont_validation.py --outdir /analysis/hg003-ont

bin/run --input /analysis/hg003/hifi.samples.csv \
  --resources /analysis/hg003/resources.json --regions /analysis/hg003/evaluation.bed \
  --outdir /analysis/hifi-run --profile apptainer --cpus 8 --memory-gb 32

bin/run --input /analysis/hg003-ont/ont.samples.csv \
  --resources /analysis/hg003/resources.json --regions /analysis/hg003-ont/evaluation.bed \
  --outdir /analysis/ont-run --profile apptainer --cpus 8 --memory-gb 32
```

Use each run's `calls.csv` with the [annotation](annotation.md) interface. Add
the complete truth VCF and confidence BED to a benchmark CSV, and supply the
corresponding `--regions` to [benchmarking](benchmarking.md). The full VEP cache
and optional custom annotation resources require the explicit preparation in
[resources](resources.md); they are not distributed in Git.

CI runs Python validation/report-parser tests and compiles the actual Nextflow
DSL2 workflows. It does not download human genomes or replace real integration
execution with mocked accuracy results. Public fixture hashes and actual
tool outputs remain outside the checkout; compact reports are committed here.

Outstanding work is tracked in [the implementation roadmap](roadmap.md): whole
genomes with an independently audited training holdout, real GIAB SV truth,
modern ONT models and FASTQ, Docker/Slurm execution, and dedicated performance
studies. Sex chromosomes, mitochondrial calling, somatic analysis, dedicated
CNVs, repeat expansions, and cohort joint genotyping are outside this release.
