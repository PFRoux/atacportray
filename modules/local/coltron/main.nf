process COLTRON {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(enhancers), path(bam), path(bai)
    val   genome_build
    path  subpeaks
    path  activity
    path  motifs
    val   enhancer_number
    val   exp_cutoff

    output:
    tuple val(meta), path("*coltron") , emit: outdir
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.coltron"
    def subpeaks_arg = subpeaks ? "--subpeaks ${subpeaks}" : ''
    def activity_arg = activity ? "--activity ${activity}" : ''
    def motifs_arg = motifs ? "--motifs ${motifs}" : ''
    def enhancer_arg = enhancer_number ? "--enhancer_number ${enhancer_number}" : ''
    def exp_arg = exp_cutoff ? "--expCutoff ${exp_cutoff}" : ''
    """
    coltron \\
        --enhancer_file ${enhancers} \\
        --bam ${bam} \\
        --genome ${genome_build} \\
        --output ${prefix} \\
        --name ${meta.id} \\
        ${subpeaks_arg} \\
        ${activity_arg} \\
        ${motifs_arg} \\
        ${enhancer_arg} \\
        ${exp_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coltron: \$(python - <<'PY'
import pkg_resources
print(pkg_resources.get_distribution("coltron").version)
PY
)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.coltron"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${meta.id}.coltron.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coltron: 1.0.2
    END_VERSIONS
    """
}
