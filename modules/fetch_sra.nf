process FETCH_SRA {
    tag "${meta.id}"
    container 'quay.io/biocontainers/sra-tools:3.1.1--h4304569_2'

    input:
    tuple val(meta), val(run_id)

    output:
    tuple val(meta), path("${run_id}_1.fastq.gz"), path("${run_id}_2.fastq.gz")

    script:
    """ 
    prefetch ${run_id}
    fasterq-dump --split-files --threads ${task.cpus} ${run_id}
    gzip ${run_id}_1.fastq &
    gzip ${run_id}_2.fastq &
    wait
    rm -r ${run_id}
    """
}