process CONSENSUS_SUPER_ENHANCERS {
    tag "consensus"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    path tables, stageAs: "rose_super_enhancers/*"

    output:
    path "consensus_super_enhancers.bed", emit: bed
    path "consensus_super_enhancers.saf", emit: saf
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk '
        BEGIN { FS=OFS="\\t" }
        FNR == 1 {
            chrom_col = start_col = end_col = 0
            is_super_col = 0
        }
        /^#/ || NF == 0 {
            next
        }
        {
            upper_first = toupper(\$1)
        }
        upper_first == "REGION_ID" || upper_first == "CHROM" || upper_first == "CHR" || upper_first == "CHROMOSOME" {
            chrom_col = start_col = end_col = 0
            is_super_col = 0
            for (i = 1; i <= NF; i++) {
                col = toupper(\$i)
                if (col == "CHROM" || col == "CHR" || col == "CHROMOSOME") chrom_col = i
                if (col == "START") start_col = i
                if (col == "STOP" || col == "END") end_col = i
                if (col == "ISSUPER" || col == "IS_SUPER") is_super_col = i
            }
            next
        }
        chrom_col && start_col && end_col && \\
            (!is_super_col || \$is_super_col == 1 || toupper(\$is_super_col) == "TRUE") && \\
            \$start_col ~ "^[0-9]+\$" && \$end_col ~ "^[0-9]+\$" {
            print \$chrom_col, \$start_col, \$end_col
        }
    ' rose_super_enhancers/* \\
        | sort -k1,1 -k2,2n \\
        | bedtools merge -i - \\
        > consensus_super_enhancers.bed

    if [ ! -s consensus_super_enhancers.bed ]; then
        echo "No consensus super-enhancer intervals were extracted from ROSE SuperStitched tables." >&2
    fi

    awk 'BEGIN{OFS="\\t"; print "GeneID","Chr","Start","End","Strand"} \\
         {print \$1":"\$2"-"\$3, \$1, \$2+1, \$3, "."}' consensus_super_enhancers.bed \\
        > consensus_super_enhancers.saf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e 's/bedtools v//g')
    END_VERSIONS
    """

    stub:
    """
    touch consensus_super_enhancers.bed
    printf "GeneID\\tChr\\tStart\\tEnd\\tStrand\\n" > consensus_super_enhancers.saf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.31.1
    END_VERSIONS
    """
}
