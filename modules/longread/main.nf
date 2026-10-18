process PREPARE_REFERENCE {
    label 'light'
    container params.images.samtools
    publishDir "${params.outdir}/reference", mode: 'copy', pattern: '*.{fai,bed}'
    publishDir "${params.outdir}/provenance/versions", mode: 'copy', pattern: '*.versions.txt'
    input:
    path fasta, stageAs: 'reference.fa'
    path regions, stageAs: 'requested.bed'
    output:
    tuple path('reference.fa'), path('reference.fa.fai'), path('calling-regions.bed'), emit: reference
    path 'reference.versions.txt'
    script:
    def select_regions = regions ? "awk 'BEGIN{OFS=\"\\t\"} !/^#/ && !/^track/ && !/^browser/ && NF>=3 {print \$1,\$2,\$3}' requested.bed > calling-regions.bed" : "awk 'BEGIN{OFS=\"\\t\"} \$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])\$/ {print \$1,0,\$2}' reference.fa.fai > calling-regions.bed"
    """
    samtools faidx reference.fa
    ${select_regions}
    test -s calling-regions.bed
    samtools --version > reference.versions.txt
    """
}

process ALIGN_FASTQ {
    tag "${meta.id}"
    container params.images.align
    publishDir "${params.outdir}/alignments/${meta.id}", mode: 'copy', pattern: '*.bam*'
    publishDir "${params.outdir}/provenance/versions", mode: 'copy', pattern: '*.versions.txt'
    input:
    tuple val(meta), path(reads, stageAs: 'input.reads')
    tuple path(fasta), path(fai), path(regions)
    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    path "${meta.id}.alignment.versions.txt"
    script:
    def preset = meta.platform == 'hifi' ? 'map-hifi' : (meta.chemistry == 'r9.4.1_450bps' ? 'map-ont' : 'lr:hq')
    def platform = meta.platform == 'hifi' ? 'PACBIO' : 'ONT'
    def align_threads = Math.max(1, task.cpus.intdiv(2))
    def sort_threads = Math.max(0, task.cpus - align_threads - 1)
    def sort_memory_mb = Math.max(256, (task.memory.toMega() * 0.5 / (sort_threads + 1)) as int)
    """
    minimap2 -ax ${preset} --MD -t ${align_threads} \\
      -R '@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:${platform}' \\
      reference.fa input.reads \\
      | samtools sort -@ ${sort_threads} -m ${sort_memory_mb}M -o ${meta.id}.bam -
    samtools quickcheck -v ${meta.id}.bam
    samtools index -@ ${Math.max(0, task.cpus - 1)} ${meta.id}.bam
    minimap2 --version > ${meta.id}.alignment.versions.txt
    samtools --version >> ${meta.id}.alignment.versions.txt
    """
}

process ACCEPT_BAM {
    tag "${meta.id}"
    container params.images.samtools
    publishDir "${params.outdir}/alignments/${meta.id}", mode: 'copy', pattern: '*.bam*'
    publishDir "${params.outdir}/provenance/versions", mode: 'copy', pattern: '*.versions.txt'
    input:
    tuple val(meta), path(reads, stageAs: 'input.bam')
    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    path "${meta.id}.alignment.versions.txt"
    script:
    """
    samtools quickcheck -v input.bam
    ln -s input.bam ${meta.id}.bam
    samtools index -@ ${Math.max(0, task.cpus - 1)} ${meta.id}.bam
    samtools --version > ${meta.id}.alignment.versions.txt
    """
}

process ALIGNMENT_QC {
    tag "${meta.id}"
    container params.images.samtools
    publishDir "${params.outdir}/qc/${meta.id}/samtools", mode: 'copy', pattern: '*.txt'
    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(regions)
    output:
    tuple val(meta), path("${meta.id}.*.txt"), emit: qc
    script:
    """
    samtools stats -@ ${Math.max(0, task.cpus - 1)} --reference reference.fa ${bam} > ${meta.id}.stats.txt
    samtools flagstat -@ ${Math.max(0, task.cpus - 1)} ${bam} > ${meta.id}.flagstat.txt
    samtools idxstats ${bam} > ${meta.id}.idxstats.txt
    """
}

process NANOPLOT {
    tag "${meta.id}"
    container params.images.nanoplot
    publishDir "${params.outdir}/qc/${meta.id}", mode: 'copy', pattern: 'NanoPlot', saveAs: { 'nanoplot' }
    publishDir "${params.outdir}/provenance/versions", mode: 'copy', pattern: '*.versions.txt'
    input:
    tuple val(meta), path(bam), path(bai)
    output:
    tuple val(meta), path('NanoPlot'), emit: qc
    path "${meta.id}.nanoplot.versions.txt"
    script:
    """
    NanoPlot --bam ${bam} --threads ${task.cpus} --outdir NanoPlot --tsv_stats --no_static
    NanoPlot --version > ${meta.id}.nanoplot.versions.txt
    """
}

process MOSDEPTH {
    tag "${meta.id}"
    container params.images.mosdepth
    publishDir "${params.outdir}/qc/${meta.id}/mosdepth", mode: 'copy', pattern: '*.mosdepth*'
    publishDir "${params.outdir}/qc/${meta.id}/mosdepth", mode: 'copy', pattern: '*.regions.bed.gz*'
    publishDir "${params.outdir}/qc/${meta.id}/mosdepth", mode: 'copy', pattern: '*.thresholds.bed.gz*'
    publishDir "${params.outdir}/provenance/versions", mode: 'copy', pattern: '*.versions.txt'
    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(regions)
    output:
    tuple val(meta), path("${meta.id}.mosdepth*"), path("${meta.id}.regions.bed.gz*"), path("${meta.id}.thresholds.bed.gz*"), emit: qc
    path "${meta.id}.mosdepth.versions.txt"
    script:
    """
    mosdepth --threads ${Math.max(0, task.cpus - 1)} --no-per-base --by calling-regions.bed \\
      --thresholds 1,5,10,20,30 ${meta.id} ${bam}
    mosdepth --version > ${meta.id}.mosdepth.versions.txt
    """
}

process DEEPVARIANT {
    tag "${meta.id}"
    container params.images.deepvariant
    publishDir "${params.outdir}/variants/${meta.id}", mode: 'copy', pattern: '*.vcf.gz*'
    publishDir "${params.outdir}/variants/${meta.id}", mode: 'copy', pattern: '*.html'
    publishDir "${params.outdir}/provenance/versions", mode: 'copy', pattern: '*.versions.txt'
    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(regions)
    output:
    tuple val(meta), path("${meta.id}.deepvariant.vcf.gz"), path("${meta.id}.deepvariant.vcf.gz.tbi"), emit: calls
    path "${meta.id}.deepvariant.g.vcf.gz*"
    path "${meta.id}.deepvariant.visual_report.html", optional: true
    path "${meta.id}.deepvariant.versions.txt"
    script:
    """
    export OMP_NUM_THREADS=1
    export TF_NUM_INTRAOP_THREADS=${task.cpus}
    export TF_NUM_INTEROP_THREADS=1
    /opt/deepvariant/bin/run_deepvariant \\
      --model_type=PACBIO --ref=reference.fa --reads=${bam} \\
      --regions=calling-regions.bed --sample_name=${meta.id} \\
      --output_vcf=${meta.id}.deepvariant.vcf.gz \\
      --output_gvcf=${meta.id}.deepvariant.g.vcf.gz \\
      --num_shards=${task.cpus} --intermediate_results_dir=intermediate \\
      --vcf_stats_report=true
    printf 'DeepVariant 1.10.0\\nmodel_type=PACBIO (bundled CPU checkpoint)\\n' > ${meta.id}.deepvariant.versions.txt
    """
}

process CLAIR3 {
    tag "${meta.id}"
    container params.images.clair3
    publishDir "${params.outdir}/variants/${meta.id}", mode: 'copy', pattern: '*.vcf.gz*'
    publishDir "${params.outdir}/provenance/versions", mode: 'copy', pattern: '*.versions.txt'
    input:
    tuple val(meta), path(bam), path(bai), path(model, stageAs: 'clair3-model')
    tuple path(fasta), path(fai), path(regions)
    output:
    tuple val(meta), path("${meta.id}.clair3.vcf.gz"), path("${meta.id}.clair3.vcf.gz.tbi"), emit: calls
    path "${meta.id}.clair3.versions.txt"
    script:
    """
    export OMP_NUM_THREADS=1
    /opt/bin/run_clair3.sh --bam_fn="\$PWD/${bam}" --ref_fn="\$PWD/reference.fa" \\
      --threads=${task.cpus} --platform=ont --model_path="\$PWD/clair3-model" \\
      --output="\$PWD/clair3-output" --bed_fn="\$PWD/calling-regions.bed" \\
      --sample_name=${meta.id}
    cp clair3-output/merge_output.vcf.gz ${meta.id}.clair3.vcf.gz
    cp clair3-output/merge_output.vcf.gz.tbi ${meta.id}.clair3.vcf.gz.tbi
    printf 'Clair3 2.0.2\\n' > ${meta.id}.clair3.versions.txt
    cat clair3-model/model.json >> ${meta.id}.clair3.versions.txt
    """
}

process SNIFFLES {
    tag "${meta.id}"
    container params.images.sniffles
    publishDir "${params.outdir}/variants/${meta.id}", mode: 'copy', pattern: '*.vcf.gz*'
    publishDir "${params.outdir}/variants/${meta.id}", mode: 'copy', pattern: '*.snf'
    publishDir "${params.outdir}/provenance/versions", mode: 'copy', pattern: '*.versions.txt'
    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(regions)
    output:
    tuple val(meta), path("${meta.id}.sniffles2.vcf.gz"), path("${meta.id}.sniffles2.vcf.gz.tbi"), emit: calls
    path "${meta.id}.sniffles2.snf"
    path "${meta.id}.sniffles2.versions.txt"
    script:
    """
    sniffles --input ${bam} --reference reference.fa --regions calling-regions.bed \\
      --vcf ${meta.id}.sniffles2.vcf.gz --snf ${meta.id}.sniffles2.snf \\
      --sample-id ${meta.id} --threads ${task.cpus} --minsvlen 50 --minsupport auto
    tabix -f -p vcf ${meta.id}.sniffles2.vcf.gz
    sniffles --version > ${meta.id}.sniffles2.versions.txt
    """
}

process PHASE {
    tag "${meta.id}"
    container params.images.whatshap
    publishDir "${params.outdir}/phasing/${meta.id}", mode: 'copy', pattern: '*.vcf.gz*'
    publishDir "${params.outdir}/phasing/${meta.id}", mode: 'copy', pattern: '*.tsv'
    publishDir "${params.outdir}/provenance/versions", mode: 'copy', pattern: '*.versions.txt'
    input:
    tuple val(meta), path(bam), path(bai), path(vcf), path(tbi)
    tuple path(fasta), path(fai), path(regions)
    output:
    tuple val(meta), path("${meta.id}.haplotagged.bam"), emit: phased
    path "${meta.id}.phased.vcf.gz*"
    path "${meta.id}.phasing.tsv"
    path "${meta.id}.whatshap.versions.txt"
    script:
    """
    python - <<'PYTHON'
    import pysam
    with pysam.VariantFile('${vcf}') as source:
        header = source.header.copy()
        with open('reference.fa.fai') as fai:
            for line in fai:
                name, length, *_ = line.split()
                if name not in header.contigs:
                    header.contigs.add(name, length=int(length))
                elif header.contigs[name].length not in (None, int(length)):
                    raise ValueError('VCF/reference contig length mismatch: ' + name)
        with pysam.VariantFile('phase-input.vcf.gz', 'wz', header=header) as target:
            for record in source:
                record.translate(header)
                target.write(record)
    pysam.tabix_index('phase-input.vcf.gz', preset='vcf', force=True)
    PYTHON
    whatshap phase --reference reference.fa --sample ${meta.id} \\
      --output ${meta.id}.phased.vcf.gz phase-input.vcf.gz ${bam}
    tabix -p vcf ${meta.id}.phased.vcf.gz
    whatshap stats --tsv ${meta.id}.phasing.tsv ${meta.id}.phased.vcf.gz
    whatshap haplotag --reference reference.fa --sample ${meta.id} \\
      --output ${meta.id}.haplotagged.bam ${meta.id}.phased.vcf.gz ${bam}
    whatshap --version > ${meta.id}.whatshap.versions.txt
    """
}

process INDEX_PHASED {
    tag "${meta.id}"
    container params.images.samtools
    publishDir "${params.outdir}/phasing/${meta.id}", mode: 'copy', pattern: '*.bam*'
    input:
    tuple val(meta), path(bam)
    output:
    tuple val(meta), path("${meta.id}.haplotagged.bam"), path("${meta.id}.haplotagged.bam.bai"), emit: bam
    script:
    """
    samtools quickcheck -v ${bam}
    samtools index -@ ${Math.max(0, task.cpus - 1)} ${bam}
    """
}
