process QC_MERGE {
    label 'analysis'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path h5ads
    path samplesheet

    output:
    path "combined_qc.h5ad", emit: h5ad

    script:
    """
    qc_merge.py --inputs ${h5ads} --samplesheet ${samplesheet} --out combined_qc.h5ad
    """

    stub:
    """
    touch combined_qc.h5ad
    """
}