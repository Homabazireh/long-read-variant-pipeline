# Long-read calling workflow

`bin/run` launches the real Nextflow DSL2 workflow in this repository. It accepts human germline HiFi and ONT data, preserves accepted BAM alignments and tags, and produces small variants, structural variants, phasing, QC, and a callset manifest for the annotation and benchmarking commands. Software execution tests and scientific accuracy benchmarks are separate claims; consult the repository's validation report for measured results.

## Requirements and invocation

Install Python 3.10+, Java 17–25, exactly Nextflow 25.10.4, and Docker or Apptainer. Use Linux x86-64 with AVX for the CPU DeepVariant image. Whole genomes require substantially more disk and runtime than interval tests. The default request is 8 CPUs and 32 GiB **per task**; the local executor also caps concurrent tasks at those totals. Slurm requests these resources per job. `--queue` selects an optional Slurm partition.

```bash
bin/run \
  --input samples.csv \
  --resources resources.json \
  --outdir /results/long-read-run \
  --profile apptainer \
  --cpus 8 --memory-gb 32
```

Add `--dry-run` to validate reference, reads, and model resources and print the launch manifest without writing files or executing tools. Add `--resume` to reuse an existing matching Nextflow session and work directory; it cannot create a fresh run. The launcher rejects unrelated nonempty output directories and locks the destination against concurrent writers. Previous trace, reports, command, input manifest, and run status are archived under `provenance/history/` before each resume. Biological inputs, reference, region file, and recorded read size/modification time must be unchanged; changed inputs require a new output directory. `--work-dir` selects an alternative persistent work directory. Avoid deleting the work directory or `provenance/.nextflow` until the run is complete and no further resume is needed.

For a bounded, real analysis, add `--regions regions.bed`. Regions are 0-based, half-open, nonoverlapping autosomal intervals. Without it, intervals are derived from the reference FAI for chromosomes 1–22. Alignment and BAM QC retain all reads; small-variant and SV callers run in selected intervals. Coverage reports include global BAM statistics and interval-specific depth. This option does not download or synthesize reads and does not turn a small region into a whole-genome validation.

## Input contract

The CSV must have `sample,platform,reads,chemistry,basecaller,model` columns. Paths to reads and model directories can be absolute or relative to the CSV. IDs must be unique and use letters, digits, dots, underscores, and dashes, starting with a letter/digit.

- `platform=hifi`: basecalled HiFi FASTQ/FASTQ.gz or coordinate-sorted aligned BAM. The bundled DeepVariant `PACBIO` checkpoint is selected; leave `model` empty. Chemistry and basecaller metadata are recorded when supplied.
- `platform=ont`: basecalled FASTQ/FASTQ.gz or coordinate-sorted aligned BAM, with all chemistry/basecaller/model fields required. Modern chemistry is `r10.4.1_e8.2_400bps` with a complete Dorado model such as `dna_r10.4.1_e8.2_400bps_sup@v5.0.0`. Legacy `r9.4.1_450bps` with `guppy_hac@v3.6.0` is also supported when a matching PyTorch model is supplied. The model must declare the exact supplied basecaller, not merely a broadly compatible flow cell.
- BAMs must have the same reference contig dictionary and lengths as the FASTA, coordinate sort order, and read groups with `SM` exactly matching the CSV sample. `PL` must be `PACBIO` or `ONT`. The BAM header is inspected before launch, then `samtools quickcheck` and indexing validate the file during execution. BAMs are never silently converted to FASTQ, realigned, or stripped of platform tags. Unaligned BAM/CRAM/POD5/FAST5 inputs are unsupported in this release.

The resource JSON requires `schema_version: 1`, `genome_build: "GRCh38"`, an absolute `reference_fasta`, and the SHA256 of that exact uncompressed FASTA. The launcher hashes the reference and checks primary-autosome lengths against GRCh38. Use the complete reference used for alignment. Mismatched BAM dictionaries are an error, including extra decoy/custom contigs. A chromosome-specific reference may be used only with a correspondingly aligned chromosome-specific input; results must be labeled accordingly.

A trusted ONT model directory contains `pileup.pt`, `full_alignment.pt`, and `model.json` (see `examples/longread-model.json`). Both checkpoint SHA256 values are verified before launch. The metadata must declare Clair3 `2.0.2`, PyTorch format, matching chemistry/basecaller, and `uses_move_table: false`. Use the [official Clair3 models and downloads](https://github.com/HKU-BAL/Clair3/tree/v2.0.2), including explicitly converted PyTorch models for legacy training configurations. TensorFlow checkpoints and move-table-dependent checkpoints are rejected. Fill the metadata from the model publisher's actual training/basecaller compatibility information; declaring a model compatible does not make it so. Only load checkpoints from trusted sources.

## Processing and outputs

1. Generate the FASTA index and caller intervals. Align FASTQ with minimap2 `map-hifi`, `lr:hq` for modern ONT, or `map-ont` for legacy ONT; retain supplementary alignments and add sample/platform read groups. Sort and index with samtools. Compatible BAM input is preserved and reindexed.
2. Run NanoPlot on the aligned BAM, samtools statistics/flagstat/idxstats, and mosdepth with 1×/5×/10×/20×/30× thresholds. NanoPlot reports read/alignment metrics from retained BAM records, rather than claiming raw-signal or pre-basecalling QC.
3. Run DeepVariant 1.10.0 `PACBIO` for HiFi or Clair3 2.0.2 with the supplied ONT model. Small-variant VCFs retain caller filters and reference representations. HiFi also produces gVCF output.
4. Run Sniffles2 2.8.1 on the selected intervals with `--minsvlen 50 --minsupport auto`. Retain its VCF, index, and SNF. These SV outputs are not a dedicated CNV or repeat-expansion assay.
5. Phase small variants with WhatsHap 2.8, report phasing statistics, and haplotag a separate BAM. Complete the phase-input VCF contig header from the validated FASTA so haplotagging can retain reads on uncalled chromosomes without silently dropping them. Preserve original caller VCFs and alignments alongside the phased products.

Output layout:

```text
alignments/SAMPLE/SAMPLE.bam[.bai]
variants/SAMPLE/SAMPLE.deepvariant.vcf.gz[.tbi]  # or .clair3.vcf.gz
variants/SAMPLE/SAMPLE.deepvariant.g.vcf.gz[.tbi] # HiFi
variants/SAMPLE/SAMPLE.sniffles2.vcf.gz[.tbi]
variants/SAMPLE/SAMPLE.sniffles2.snf
phasing/SAMPLE/SAMPLE.phased.vcf.gz[.tbi]
phasing/SAMPLE/SAMPLE.haplotagged.bam[.bai]
phasing/SAMPLE/SAMPLE.phasing.tsv
qc/SAMPLE/{nanoplot,samtools,mosdepth}/
reference/{reference.fa.fai,calling-regions.bed}
provenance/{run-inputs.json,run-status.json,command.json,nextflow.log,trace.tsv,execution-report.html,timeline.html,workflow-dag.html,versions/,history/}
calls.csv
```

`calls.csv` has `sample,label,vcf,variant_type` columns and absolute paths. It includes original small-variant and Sniffles2 callsets, so annotation and benchmarks see caller output without phasing-specific transformations. Headers-only VCFs remain valid outputs. The manifest may be created before concurrent publishing finishes; consume it only after Nextflow exits successfully.

## Software pins and checks

The fixed image references are in `configs/longread-containers.config`. Runtime commands and tool versions are recorded in the output provenance. `run-status.json` records completion/failure, elapsed time, Nextflow version, Git revision/dirty state, and SHA256 hashes of the workflow, modules, configuration, and launcher for each attempt. Pins were obtained from the [nf-core module implementations](https://github.com/nf-core/modules/tree/master/modules/nf-core) for minimap2 2.30/samtools 1.23.1, samtools/htslib 1.24, NanoPlot 1.47.0, mosdepth 0.3.14, Sniffles2 2.8.1 and WhatsHap 2.8; the caller images follow [DeepVariant 1.10.0](https://github.com/google/deepvariant/blob/r1.10/docs/deepvariant-quick-start.md) and [Clair3 2.0.2](https://github.com/HKU-BAL/Clair3/tree/v2.0.2). Container tags are fixed; registry tags may still be mutable, so archive pulled images/digests with analysis records for long-term reproducibility.

```bash
python3 -m unittest discover -s tests -p 'test_longread*.py' -v
nextflow run main.nf --help -profile apptainer
# Assert real outputs after a completed run (unset by default in unit-test CI):
LONGREAD_TEST_OUTDIR=/results/long-read-run python3 -m unittest discover -s tests -p 'test_longread*.py' -v
```

Unit tests verify malformed inputs, reference/build mismatches, BAM sample/platform/dictionary checks, ONT metadata and checkpoint hashes, BED bounds and overlap, dry-run behavior, and safe resume. They do not invoke pretend callers or count mocked outputs as integration results. Real integration requires public reads, matched reference/model resources, and the actual pinned images. Whole-genome performance, ONT model configurations not exercised by those data, sex chromosomes, mitochondrial calls, somatic analysis, dedicated CNV calls, and joint genotyping require separate validation.
