process FASTQ_CONCAT {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(reads, stageAs: 'fastq??/*')

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val("${task.process}"), val('cat'), eval("printf 'unknown'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_files = reads instanceof List ? reads : [reads]

    if (meta.single_end) {
        """
        cat ${read_files.join(' ')} > ${prefix}.fastq.gz
        """
    } else {
        if (read_files.size() % 2 != 0) {
            error "FASTQ_CONCAT expected an even number of FASTQ files for paired-end sample '${meta.id}', got ${read_files.size()}."
        }
        def r1_reads = []
        def r2_reads = []
        read_files.eachWithIndex { read, index ->
            if (index % 2 == 0) {
                r1_reads << read
            } else {
                r2_reads << read
            }
        }
        """
        cat ${r1_reads.join(' ')} > ${prefix}_R1.fastq.gz
        cat ${r2_reads.join(' ')} > ${prefix}_R2.fastq.gz
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        echo '' | gzip > ${prefix}.fastq.gz
        """
    } else {
        """
        echo '' | gzip > ${prefix}_R1.fastq.gz
        echo '' | gzip > ${prefix}_R2.fastq.gz
        """
    }
}
