nextflow.enable.dsl=2

include { PREPARE_REFERENCE; ALIGN_FASTQ; ACCEPT_BAM; ALIGNMENT_QC; NANOPLOT; MOSDEPTH; DEEPVARIANT; CLAIR3; SNIFFLES; PHASE; INDEX_PHASED } from './modules/longread/main'

def csvRow(values) {
    values.collect { value -> '"' + value.toString().replace('"', '""') + '"' }.join(',')
}

workflow {
    if (params.help) {
        log.info 'Launch with bin/run --input samples.csv --resources resources.json --outdir results --profile apptainer. See docs/long-read-workflow.md.'
        return
    }
    if (!params.samples || !params.reference_fasta || !params.reference_sha256) {
        error 'Use bin/run to validate inputs and generate the Nextflow parameter manifest.'
    }
    if (!(workflow.containerEngine in ['docker', 'apptainer'])) {
        error 'Use the docker or apptainer execution profile.'
    }
    input_samples = Channel.fromList(params.samples).map { sample ->
        tuple(sample, file(sample.reads, checkIfExists: true))
    }
    inputs = input_samples.branch {
        fastq: it[0].kind == 'fastq'
        bam: it[0].kind == 'bam'
    }
    PREPARE_REFERENCE(file(params.reference_fasta, checkIfExists: true), params.regions ? file(params.regions, checkIfExists: true) : [])
    reference = PREPARE_REFERENCE.out.reference
    ALIGN_FASTQ(inputs.fastq, reference)
    ACCEPT_BAM(inputs.bam)
    alignments = ALIGN_FASTQ.out.bam.mix(ACCEPT_BAM.out.bam)
    ALIGNMENT_QC(alignments, reference)
    NANOPLOT(alignments)
    MOSDEPTH(alignments, reference)
    platforms = alignments.branch {
        hifi: it[0].platform == 'hifi'
        ont: it[0].platform == 'ont'
    }
    DEEPVARIANT(platforms.hifi, reference)
    ont_models = platforms.ont.map { meta, bam, bai -> tuple(meta, bam, bai, file(meta.model, checkIfExists: true)) }
    CLAIR3(ont_models, reference)
    small_calls = DEEPVARIANT.out.calls.mix(CLAIR3.out.calls)
    SNIFFLES(alignments, reference)
    phase_inputs = alignments.join(small_calls)
    PHASE(phase_inputs, reference)
    INDEX_PHASED(PHASE.out.phased)
    call_rows = small_calls.map { meta, vcf, index ->
        def label = meta.platform == 'hifi' ? 'deepvariant' : 'clair3'
        csvRow([meta.id, label, "${file(params.outdir).toAbsolutePath()}/variants/${meta.id}/${vcf.name}", "small"])
    }.mix(SNIFFLES.out.calls.map { meta, vcf, index ->
        csvRow([meta.id, "sniffles2", "${file(params.outdir).toAbsolutePath()}/variants/${meta.id}/${vcf.name}", "sv"])
    })
    call_rows.collectFile(name: 'calls.csv', storeDir: params.outdir, seed: 'sample,label,vcf,variant_type', newLine: true, sort: true)
}
