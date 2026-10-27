process KB_REF {
    container 'ghcr.io/pranavasri-rm/kb-python:0.30.2'

    input:
    path genome_fasta
    path genome_gtf
    path virus_fasta
    path virus_gtf

    output:
    path "index.idx", emit: index
    path "t2g.txt",   emit: t2g

    script:
    """
    kb ref \\
        -i index.idx \\
        -g t2g.txt \\
        -f1 cdna.fa \\
        -t ${task.cpus} \\
        ${genome_fasta},${virus_fasta} \\
        ${genome_gtf},${virus_gtf}
    """

    stub:
    """
    touch index.idx t2g.txt
    """
}