process QC_SAMPLE {
    tag "${meta.id}"
    label 'analysis'

    input:
    tuple val(meta), path(kb_dir)
    path cutoffs

    output:
    tuple val(meta), path("${meta.id}.h5ad"), emit: h5ad

    script:
    """
    qc_sample.py \\
        --counts-dir ${kb_dir}/counts_unfiltered \\
        --sample ${meta.id} \\
        --donor ${meta.donor} \\
        --condition ${meta.condition} \\
        --cutoffs ${cutoffs} \\
        --out ${meta.id}.h5ad
    """
}