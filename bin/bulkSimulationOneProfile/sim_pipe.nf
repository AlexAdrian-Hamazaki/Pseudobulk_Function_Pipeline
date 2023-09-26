
process makeCTProfiles {
    publishDir "data/CTProfiles", mode: 'copy'
    conda params.python3_9
    
    input:
        path h5ad_path
    output:
        path "${h5ad_path.getSimpleName()}_cell_type_profiles.csv"
    shell:
    '''
    makeCellTypeProfiles.py !{h5ad_path} !{h5ad_path.getSimpleName()}
    '''
}

process simBulk {
    publishDir "data/simulations/${variance_ch}/${cell_type_profiles_csv}", mode: "copy"
    conda params.python3_9
    
    input:
        tuple path(cell_type_profiles_csv), path(cell_type_proportions)
        each variance_ch
        
    output:
        tuple val(variance_ch), path("${cell_type_profiles_csv.getSimpleName()}.csv.gz")
        path "${cell_type_profiles_csv.getSimpleName()}_n_sim_*_profiles.csv"

    shell:
    '''
    simBulk.py !{cell_type_profiles_csv} !{cell_type_profiles_csv.getSimpleName()} !{cell_type_proportions} !{variance_ch}
    '''
}


process sim_bulk_EGAD {
    cpus 4
    memory '32 GB'
    maxForks 12

    publishDir "data/EGAD/${variance_ch}/${expression_matrix.getBaseName()}", mode: 'copy'

    input:
        tuple val(variance_ch), path(expression_matrix), path(go_annotations_ch)
        val go_annot_column
    output:
        path "*_EGAD.csv"
    shell:
    '''
    simBulkEgad.R !{expression_matrix} !{expression_matrix.getBaseName()} !{go_annotations_ch} !{go_annotations_ch.getSimpleName()} !{go_annot_column}
    '''
}


tissue_ch = Channel.fromPath(params.tissue_dir)
cell_type_proportions = Channel.fromPath(params.cell_type_proportions)
variance_ch = Channel.value(params.compositional_variance)

// Go annotations is all of the GO terms and what genes are affiliated with them
go_annotations_ch = Channel.fromPath(params.go_annotations)
go_annot_gene_column_ch = Channel.value(params.gene_column)

workflow {
    // First Make cell type profiles for each tissue
    makeCTProfiles(tissue_ch)


    tuple_ch = makeCTProfiles.out.combine(cell_type_proportions)

    // Simulate a bunch of bulk profiles at different variances
    simBulk(tuple_ch, variance_ch)

    tuple_ch2 = simBulk.out[0].combine(go_annotations_ch)
    //tuple_ch2.view()


    // Run EGAD on the different samples
    sim_bulk_EGAD(tuple_ch2, go_annot_gene_column_ch)


}