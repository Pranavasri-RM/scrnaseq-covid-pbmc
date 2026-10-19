include { FETCH_SRA } from './modules/fetch_sra'

params.samplesheet = "${projectDir}/assets/samplesheet.csv"

workflow {
    samples_ch = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample, donor: row.donor, condition: row.condition]
            tuple(meta, row.run_accession) 
        }
        
    FETCH_SRA(samples_ch)
    FETCH_SRA.out.view()
}