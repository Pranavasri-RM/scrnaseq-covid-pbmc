params.samplesheet = "${projectDir}/assets/samplesheet.csv"

workflow {
    samples_ch = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample, donor: row.donor, condition: row.condition]
            tuple(meta, row.run_accession) }
        
    samples_ch.view()
}