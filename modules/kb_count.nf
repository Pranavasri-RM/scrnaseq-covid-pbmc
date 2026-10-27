process KB_COUNT {
    tag "${meta.id}"
    container 'ghcr.io/pranavasri-rm/kb-python:0.30.2'
    publishDir "${params.outdir}/kb_count", mode: 'copy'

    input:
    tuple val(meta), path(read1), path(read2)
    path index
    path t2g

    output:
    tuple val(meta), path("${meta.id}"), emit: counts
    tuple val(meta), path("${meta.id}/counts_unfiltered/adata.h5ad"), emit: h5ad

    script:
    """
    kb count \\
        -i ${index} \\
        -g ${t2g} \\
        -x 0,0,12:0,12,20:1,0,0 \\
        -o ${meta.id} \\
        -t ${task.cpus} \\
        -m ${Math.round(task.memory.toGiga() * 0.6)}G \\
        --h5ad \\
        --delete-bus \\
        ${read1} ${read2}
    """

    
    stub:
    """
    mkdir -p ${meta.id}/counts_unfiltered
    touch ${meta.id}/counts_unfiltered/adata.h5ad
    """
}