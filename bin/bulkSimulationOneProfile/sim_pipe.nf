
process makeCTProfiles {
    publishDir "${params.publish}/CTProfiles", mode: 'copy'

    conda params.python3_9 
    
    input:
        path h5ad_path
    output:
        path "exp_${h5ad_path.getSimpleName()}_cell_type_profiles.csv"
        path "cntrl_${h5ad_path.getSimpleName()}_cell_type_profiles.csv"
    shell:
    '''
    makeCellTypeProfiles.py !{h5ad_path} !{h5ad_path.getSimpleName()}
    '''
}

process graphSparsity {
    publishDir "${params.publish}/graphs", mode: "copy"
    conda params.python3_9

    input:
        each path(ct_profile_csv)
    output:
        path "${ct_profile_csv.getSimpleName()}_sparsity_graph.png"
    script:
        """
        graphSparsity.py ${ct_profile_csv} ${ct_profile_csv.getSimpleName()}
        """
}

process simBulk {
    publishDir "${params.publish}/simulations/${variance_ch}/${cell_type_profiles_csv}", mode: "copy"
    conda params.python3_9
    
    input:
        tuple path(cell_type_profiles_csv), path(cell_type_proportions)
        each variance_ch
        val num_simulations_ch
        
    output:
        tuple val(variance_ch), path("${cell_type_profiles_csv.getSimpleName()}_.csv.gz"), path("${cell_type_profiles_csv.getSimpleName()}_n_sim_*_profiles.csv") // This will output a list of these files
        tuple val(variance_ch), path("${cell_type_profiles_csv.getSimpleName()}_n_sim_*_profiles_uncollapsed.csv")
    shell:
    '''
    simBulk.py !{cell_type_profiles_csv} !{cell_type_profiles_csv.getSimpleName()} !{cell_type_proportions} !{variance_ch} !{num_simulations_ch}
    '''
}


process sim_bulk_EGAD {
    memory '32 GB'
    executor "local"
    maxForks 8

    publishDir "${params.publish}/EGAD/${expression_matrix.getBaseName()}", mode: 'copy'

    input:
        tuple val(variance_ch), path(expression_matrix), path(lo_simulation_compositions), path(go_annotations_ch)
        val go_annot_column
    output:
        path "*_EGAD.csv"
    shell:
    '''
    simBulkEgad.R !{expression_matrix} !{expression_matrix.getSimpleName()} !{go_annotations_ch} !{go_annotations_ch.getSimpleName()} !{go_annot_column} !{variance_ch}
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

process makeMeltedDF {
    // For each organism part, make a melted DF containing all the EGAD results for graphing

    publishDir "${params.publish}/EGAD/"
    conda params.python3_9

    input:
        path lo_EGAD_expression_matrixes
    output:
        path "*_melted_df.csv" // Melted DF for graphing purpose
    script:
    """
    makeMeltedDf.py ${lo_EGAD_expression_matrixes}
    """
}

process graphEGADavg {
    // Graph the EGAD average AUC over variance levels

    publishDir "${params.publish}/graphs/EGAD/"
    conda params.python3_9

    input:
        each path(melted_df_csv)
    output:
        path "*_EGAD_avg_AUC_boxplot.png" // Graphs with avg for Control vs Experimental across variance levels
        path "*_EGAD_avg_AUC_lineplot.png" // Graphs with avg for Control vs Experimental across variance levels
      
    script:
    """
    graphEGAD_avg.py ${melted_df_csv}
    """
}

process graphEGAD_CTrelateness {
    // Graph the EGAD average AUC over variance levels and stratifying the AUCs by their GO CT Relatedness

    publishDir "${params.publish}/graphs/EGAD/CTAffiliation"
    conda params.python3_9

    input:
        each path(melted_df_csv)
        path ctaffiliation_csv
    output:
        path "*_EGAD_avg_AUC_boxplot_CTAffiliation.png" // Graphs with avg for Control vs Experimental across variance levels
        path "*_EGAD_avg_AUC_lineplot_CTAffiliation.png" // Graphs with avg for Control vs Experimental across variance levels
    
    script:
    """
    graphEGAD_splitCTAffiliation.py ${melted_df_csv} ${ctaffiliation_csv}
    """
}


process graphEGAD_BulkPerformers {
    // Graph the EGAD average AUC over variance levels and stratifying the AUCs by their GO CT Relatedness

    publishDir "${params.publish}/graphs/EGAD/Bulk"
    conda params.python3_9

    input:
        tuple val(organism_part), path(melted_egad), path(bulk_egad)
    output:
        path "*_EGAD_avg_AUC_boxplot_CTAffiliation.png" // Graphs with avg for Control vs Experimental across variance levels
        path "*_EGAD_avg_AUC_lineplot_CTAffiliation.png" // Graphs with avg for Control vs Experimental across variance levels
        path "*_EGAD_recovery_graph.png"
    
    script:
    """
    graphEGAD_splitCTAffiliation.py ${melted_df_csv} ${ctaffiliation_csv}
    """
}



tissue_ch = Channel.fromPath(params.tissue_dir)
cell_type_proportions = Channel.fromPath(params.cell_type_proportions)
variance_ch = Channel.value(params.compositional_variance)


// Go annotations is all of the GO terms and what genes are affiliated with them
go_annotations_ch = Channel.fromPath(params.go_annotations)
go_annot_gene_column_ch = Channel.value(params.gene_column)
ctaffiliation_csv_ch = Channel.fromPath(params.ctaffiliation_csv)

// Num sims is how many simulations we want to run
num_simulations_ch = Channel.value(params.num_sims)

//PAth to the Bulk performance GO Terms for BP
bulk_egad_ch = Channel.fromPath(params.bulk_performance)
bulk_egad_ch = bulk_egad_ch.map{[it.toString().split("/")[-1].split("_")[0], it]}

workflow {
    // First Make cell type profiles for each tissue
    makeCTProfiles(tissue_ch)
    // Combine all the CT Profiles into 1 channel
    CT_profiles_ch = makeCTProfiles.out[0].concat(makeCTProfiles.out[1])
    // Graph the sparsity of the CT Profiles
    graphSparsity(CT_profiles_ch)
    
    // Add metadata pertaining to the desired cell_type_proportions
    tuple_ch = CT_profiles_ch.combine(cell_type_proportions)
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
    // Graph the Correlations of each bulk dataset
    // Create tuple to run EGAD with metadata
    tuple_ch2 = simBulk.out[0].combine(go_annotations_ch)
    
    //Run EGAD on the different samples
    sim_bulk_EGAD(tuple_ch2, go_annot_gene_column_ch)

    // Make 1 channel for EGAD graphing that passes all of the EGAD results into 1 channel
    egad_melt_ch = sim_bulk_EGAD.out.collect()
    
    // Make melted dataframe for graphing for each organism part
    makeMeltedDF(egad_melt_ch)

    // Graph the AVG AUC of each organism part across the variance levels, stratify by experiemnt and control
    graphEGADavg(makeMeltedDF.out)

    // Graph the AVG AUC of each organism part across the variance levels, stratify by experiment and control AND by if the GO term is CT Related or not
    graphEGAD_CTrelateness(makeMeltedDF.out, ctaffiliation_csv_ch)

    // Graph the avg AUC for the GTEX BULK AUC GO term performers
    
    // // For every element of this channel, convert it to a string, split in pieces separated by --, get the second part, then split by _3p and get the first part. Return a list with this as the first value, and then the original element as the second value. This part has to be customized depending on what part of the String you want to get as matching key    
    // makeMeltedDF.out.flatten().branch{brain: it.toString().split("/")[-1].split("_")[0] == 'brain'
    //                         blood: it.toString().split("/")[-1].split("_")[0] == 'pbmc'
    //                         other: true
    //                         }.set{melted_ch}
    // melted_ch.blood.map{['Blood', it]}.set{blood_ch}
    // melted_ch.brain.map{['Brain', it]}.set{brain_ch}
    // ch_collected_melts = brain_ch.concat(blood_ch)

    // ch_bulk_sc_egads = ch_collected_melts.join(bulk_egad_ch, by:0).view()
    

    // Combine according to a key that is the first value of every first element, which is a list according to what we did above
    // For every element of this channel, which consists of three values now, the matching key (id), the first element of the first channel, and the second, keep only the second and the third.
// https://nextflow-io.github.io/patterns/create-key-to-combine-channels/
}