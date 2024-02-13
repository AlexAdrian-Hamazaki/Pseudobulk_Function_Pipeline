
process simulateBulk { // Sample single cells and create simulated bulk samples

    publishDir "${params.publish}/non_collapsed/${tissue_ch.getSimpleName()}/sample_size${sample_size_ch}/variance${variance_ch}/", mode: 'copy'
    conda params.python3_9 
    cpus 4
    memory '32 GB'
    maxForks 12

    input:
        val sample_size_ch // int, number of cells to sample for each bulk simulated sample
        val num_simulations_ch // the number of bulk sample simulations per variance level 
        path tissue_ch // path to h5ad
        each variance_ch // list of variance levels
    output:
        tuple val("${tissue_ch}"), val(sample_size_ch), val(num_simulations_ch), val(variance_ch), path("n_sim_*.h5ad")
        path "n_sim_*.csv" // This is the metadata for each simulation
    shell:
    '''
    simBulk.py !{tissue_ch} !{tissue_ch.getSimpleName()} !{num_simulations_ch} !{sample_size_ch} !{variance_ch}
    '''

}

process combineProfiles {

    publishDir "${params.publish}/collapsed/${tissue_ch}/sample_size${sample_size_ch}/variance${variance_ch}/", mode: 'copy'
    conda params.python3_9

    input:
    
        tuple val(tissue_ch), val(sample_size_ch), val(num_simulations_ch), val(variance_ch), path("n_sim_*.h5ad")

    output:
        tuple val(tissue_ch), val(sample_size_ch), val(num_simulations_ch), val(variance_ch), path("n_sim_*.h5ad.csv")
    shell:
    '''
    combineProfiles.py n_sim_*.h5ad
    '''
}

process combineSimultedBulk {

    publishDir "${params.publish}/simulated_bulk_dataset/${tissue_ch}/numbsimulations${num_simulations_ch}/num_cells${sample_size_ch}", mode: 'copy'
    conda params.python3_9

    input: 
        tuple val(tissue_ch), val(sample_size_ch), val(num_simulations_ch), val(variance_ch), path("n_sim_*.h5ad.csv")

    output:
        tuple val(tissue_ch), val(sample_size_ch), path("simulated_ss${sample_size_ch}_var${variance_ch}_nsim${num_simulations_ch}.csv")
    script:
    """
    combineSimulatedBulk.py ${sample_size_ch} ${variance_ch} ${num_simulations_ch} n_sim_*.h5ad.csv
    """
}

process sim_bulk_EGAD {
    cpus 4
    memory '32 GB'
    maxForks 12

    publishDir "${params.publish}/EGAD/${tissue_ch}/numbsimulations${num_simulations_ch}/ss${sample_size_ch}", mode: 'copy'

    input:
        tuple val(tissue_ch), val(sample_size_ch), path(expression_matrix)
        each go_annotations_ch
        val go_annot_column
    output:
        path "*_EGAD.csv"
    shell:
    '''
    simBulkEgad.R !{expression_matrix} !{expression_matrix.getBaseName()} !{go_annotations_ch} !{go_annotations_ch.getSimpleName()} !{go_annot_column}
    '''
}

sample_size_ch = Channel.value(params.sample_size)
num_simulations_ch = Channel.value(params.num_simulations)
tissue_dir_ch = Channel.fromPath(params.tissue_dir)
variance_ch = Channel.value(params.compositional_variance)

// Go annotations is all of the GO terms and what genes are affiliated with them
go_annotations_ch = Channel.fromPath(params.go_annotations)
go_annot_gene_column_ch = Channel.value(params.gene_column)
sample_size_ch.view()
num_simulations_ch.view()
tissue_dir_ch.view()
variance_ch.view()
workflow{
    //  For each h5ad, make a h5ad with only our specific subsample. Make sure the file is named such that each sim is traceable
    simulateBulk(sample_size_ch, num_simulations_ch, tissue_dir_ch, variance_ch)
    //simulateBulk.out[0].view()

    // // Combine the subsamples into 1 profile to "simulate a bulk"
    // combineProfiles(simulateBulk.out[0])
    // //combineProfiles.out.view()

    // // Combine simulated bulk
    // combineSimultedBulk(combineProfiles.out)
    // //combineSimultedBulk.out.view()

    // // run EGAD
    // sim_bulk_EGAD(combineSimultedBulk.out, go_annotations_ch, go_annot_gene_column_ch)

}


