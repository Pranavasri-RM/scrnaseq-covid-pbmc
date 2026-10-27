process COMPOSITION_DE {
    label 'analysis'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path h5ad
    val gene_sets   // .gmt inside the container

    output:
    path "de", emit: tables

    script:
    """
    composition_de.py --input ${h5ad} --outdir de --gene-sets ${gene_sets}
    """

    stub:
    """
    mkdir -p de
    touch de/de_summary.csv
    """
}