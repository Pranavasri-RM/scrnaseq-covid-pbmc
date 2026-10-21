include { FETCH_SRA } from './modules/fetch_sra'
include { KB_REF    } from './modules/kb_ref'
include { KB_COUNT  } from './modules/kb_count'

workflow {
    samples_ch = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample, donor: row.donor, condition: row.condition]
            tuple(meta, row.run_accession)
        }

    FETCH_SRA(samples_ch)

    KB_REF(
        file(params.genome_fasta),
        file(params.genome_gtf),
        file(params.virus_fasta),
        file(params.virus_gtf)
    )

    KB_COUNT(FETCH_SRA.out, KB_REF.out.index, KB_REF.out.t2g)
    KB_COUNT.out.counts.view()
}