include { FETCH_SRA } from './modules/fetch_sra'
include { KB_REF    } from './modules/kb_ref'
include { KB_COUNT  } from './modules/kb_count'
include { QC_SAMPLE } from './modules/qc_sample'
include { QC_MERGE  } from './modules/qc_merge'
include { CLUSTER   } from './modules/cluster'
include { ANNOTATE  } from './modules/annotate'
include { COMPOSITION_DE } from './modules/composition_de'

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
    
    QC_SAMPLE(KB_COUNT.out.counts, file(params.qc_cutoffs))

    QC_MERGE(
        QC_SAMPLE.out.h5ad.map { meta, h5ad -> h5ad }.collect(),
        file(params.samplesheet)
    )

    CLUSTER(QC_MERGE.out.h5ad)

    ANNOTATE(CLUSTER.out.h5ad, file(params.cluster_labels), params.celltypist_model)

    COMPOSITION_DE(ANNOTATE.out.h5ad, params.gene_sets)
}