process CLUSTER {
    label 'analysis'
    publishDir "${params.outdir}/clustering", mode: 'copy'

    input:
    path h5ad

    output:
    path "clustered.h5ad", emit: h5ad

    script:
    """
    cluster.py --input ${h5ad} --out clustered.h5ad
    """
}