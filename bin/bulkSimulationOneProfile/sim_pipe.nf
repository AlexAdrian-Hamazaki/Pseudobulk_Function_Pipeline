
process makeCTProfiles {
    publishDir "${params.publish}/CTProfiles", mode: 'copy'
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
    cpus 2
    memory '32 GB'
    executor "local"
    maxForks 1

    publishDir "${params.publish}/simulations/${variance_ch}/${cell_type_profiles_csv}", mode: "copy"
    conda params.python3_9
    
    input:
        tuple path(cell_type_profiles_csv), path(cell_type_proportions)
        each variance_ch
        val num_simulations_ch
        
    output:
        tuple val(variance_ch), path("${cell_type_profiles_csv.getSimpleName()}.csv.gz"), path("${cell_type_profiles_csv.getSimpleName()}_n_sim_*_profiles.csv") // This will output a list of these files
        tuple val(variance_ch), path("${cell_type_profiles_csv.getSimpleName()}_n_sim_*_profiles_uncollapsed.csv")
    shell:
    '''
    simBulk.py !{cell_type_profiles_csv} !{cell_type_profiles_csv.getSimpleName()} !{cell_type_proportions} !{variance_ch} !{num_simulations_ch}
    '''
}


process sim_bulk_EGAD {
    cpus 4
    memory '32 GB'
    maxForks 12

    publishDir "${params.publish}/EGAD/${variance_ch}/${expression_matrix.getBaseName()}", mode: 'copy'

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

process getSimulationStats {
    publishDir "${params.publish}/stats/avg_composition/${variance_ch}/${expression_matrix.getBaseName()}", mode: "copy"
    conda params.python3_9
    
    input:
        tuple val(variance_ch), path(expression_matrix), path(lo_simulation_compositions)
    output:
        tuple val(variance_ch), val("${expression_matrix.getSimpleName()}"), path("all_simulation_${variance_ch}_${expression_matrix.getSimpleName()}_proportions.csv"), path("average_simulation_${variance_ch}_${expression_matrix.getSimpleName()}_proportion.csv")
    script:
    """
    getSimulationStats.py ${variance_ch} ${expression_matrix.getSimpleName()} ${lo_simulation_compositions} 
    """
}

process graphSimulationComposition {
    // Get a graph that shows the average composition proportion across all variances for an organism part
    publishDir "${params.publish}/graphs/${variance_ch}/${op_name}", mode: "copy"
    conda params.python3_9
    
    input:
        tuple val(variance_ch), val(op_name), path(all_simulation_var_opname_proportions), path(average_simulation_var_opname_proportion)
    output:
        path "simulation_consistancy.png"

    script:
    """
    graphSimulationComposition.py ${all_simulation_var_opname_proportions}
    """
}

process graphAverageComposition {
    // Get a graph that shows the average composition proportion across all variances for an organism part
    publishDir "${params.publish}/graphs/average_composition/${organism_part_name}", mode: "copy"
    conda params.python3_9
    
    input:
        tuple val(lo_variances), val(organism_part_name), path(lo_all_simulation_proportions), path(lo_average_simulation_proportion)
    output:
        path "avg_composition_graph.png"

    script:
    """
    graphAverageComposition.py ${lo_average_simulation_proportion}
    """
}


tissue_ch = Channel.fromPath(params.tissue_dir)
cell_type_proportions = Channel.fromPath(params.cell_type_proportions)
variance_ch = Channel.value(params.compositional_variance)


// Go annotations is all of the GO terms and what genes are affiliated with them
go_annotations_ch = Channel.fromPath(params.go_annotations)
go_annot_gene_column_ch = Channel.value(params.gene_column)

// Num sims is how many simulations we want to run
num_simulations_ch = Channel.value(params.num_sims)

workflow {
    // First Make cell type profiles for each tissue
    makeCTProfiles(tissue_ch)

    tuple_ch = makeCTProfiles.out.combine(cell_type_proportions)

    // Simulate a bunch of bulk profiles at different variances
    simBulk(tuple_ch, variance_ch, num_simulations_ch)

    // Get a bunch of stats about the simulation
    getSimulationStats(simBulk.out[0]) // For one level of variance, get the average of the composition across simulations
    // Make one graph per organism part and variance level that shows cell sampling at each simulation
    graphSimulationComposition(getSimulationStats.out)

    // Create a new channel. Its a grouped tuple, keys are the organism part
    grouped_tuple_ch = getSimulationStats.out.groupTuple(by:1)
    // Make one graph per organism part that shows the average composition for each variance level
    graphAverageComposition(grouped_tuple_ch)

    // Create tuple to run EGAD with metadata
    tuple_ch2 = simBulk.out[0].combine(go_annotations_ch)
    //tuple_ch2.view()


    // Run EGAD on the different samples
    //sim_bulk_EGAD(tuple_ch2, go_annot_gene_column_ch)


}