process ANNOTATE {
    label 'analysis'
    publishDir "${params.outdir}/annotation", mode: 'copy'

    input:
    path h5ad
    path labels
    val model   // path inside the container, so a value rather than a staged file

    output:
    path "annotated.h5ad",              emit: h5ad
    path "cell_counts_per_sample.csv",  emit: counts

    script:
    """
    annotate.py --input ${h5ad} --labels ${labels} --model ${model} \\
        --out annotated.h5ad --counts-out cell_counts_per_sample.csv
    """
}