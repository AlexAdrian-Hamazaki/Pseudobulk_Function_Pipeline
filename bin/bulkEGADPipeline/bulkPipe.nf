


process splitMerged  {
    publishDir "${params.publish}/splits", mode: 'copy'
    conda params.python3_9 

    input:
        path bulk_merged_ch
        each organism_parts_ch
        
    output:
        path "${organism_parts_ch}_split.csv.gz"
    
    script:
    """
    splitadata.py ${bulk_merged_ch} ${organism_parts_ch}
    """
}

process sim_bulk_EGAD {
    memory '8 GB'
    executor "local"
    maxForks 64


    publishDir "${params.publish}/EGAD", mode: 'copy'

    input:
        tuple val(bootstrap), val(organism_part), path(expression_matrix), path(go_annotations_ch)
        val go_annot_column
    output:
        tuple val(bootstrap), val(organism_part), path("${expression_matrix.getSimpleName()}_${go_annotations_ch.getSimpleName()}_${bootstrap}_EGAD.csv")
    shell:
    '''
    bulkEgad.R !{expression_matrix} !{expression_matrix.getSimpleName()} !{go_annotations_ch} !{go_annotations_ch.getSimpleName()} !{go_annot_column} !{bootstrap}
    '''
}

process convertGtexToCSV {
    // publishDir <path>, mode: 'copy'
    conda params.python3_9 
    
    input:
        path gtex_h5ad
    output:
        path 'gtex.csv.gz'
    shell:
    '''
    adata_to_csv.py !{gtex_h5ad}
    '''
}

process makeMeltedMergedDF {
    publishDir "${params.publish}/EGAD/melted_dfs", mode: 'copy'
    conda params.python3_9 

    input:
        tuple val(lo_bootstraps), val(organism_part), path(lo_EGADs)
    output:
        path "${organism_part}_melted_EGADs.csv.gz"
    script:
    """
    make_merged_EGAD_of_splits.py --lo_bootstraps "${lo_bootstraps}" --lo_EGAD_paths "${lo_EGADs}" --organism_part ${organism_part}
    """
}

process sampleWithReplacement {
    publishDir "${params.publish}/splits/subsamples", mode: 'copy'
    conda params.python3_9 

    input:
        val bootstrap
        each path(organism_part_samples)
        val sample_size 
    output:
        tuple val(bootstrap), val("${organism_part_samples.getBaseName()}"), path("${organism_part_samples.getBaseName()}_${bootstrap}.csv.gz")
    script:
        """
        subsampleBulk.py ${bootstrap} ${organism_part_samples} ${organism_part_samples.getBaseName()} ${sample_size}
        """
}

bulk_merged_ch = Channel.fromPath(params.bulk_merged)
organism_parts_ch = Channel.fromList(params.organism_parts)

// Go annotations is all of the GO terms and what genes are affiliated with them
go_annotations_ch = Channel.fromPath(params.go_annotations)
go_annot_gene_column_ch = Channel.value(params.gene_column)

// Number of bulk samples to take for each simulation
bulk_sim_size_ch = Channel.value(params.bulk_size)

workflow bulk_pipe {
    take:
        bulk_merged_ch
        organism_parts_ch
        go_annotations_ch
        go_annot_gene_column_ch
        bootstrap_ch
        bulk_sim_size_ch

    main:

        // First  split the merged gtex dataset into different organism parts of interest as csv.gz files
        splitMerged(bulk_merged_ch, organism_parts_ch)

        // splitMerged.out.view()

        // Sample N Bulk samples with replacement for each organism part. The output of this will be dataframes of N samples for each Organism Part
        sampleWithReplacement(bootstrap_ch, splitMerged.out, bulk_sim_size_ch)
        // Make tuple containing the Bulk datasets as well as the required GO Annotations
        split_gtex_tuple_ch = sampleWithReplacement.out.combine(go_annotations_ch)
        // split_gtex_tuple_ch.view()

        // Then run EGAD on the splits and the whole EGAD
        sim_bulk_EGAD(split_gtex_tuple_ch, go_annot_gene_column_ch)

        // sim_bulk_EGAD.out.view()

        grouped_by_op_ch = sim_bulk_EGAD.out.groupTuple(by:1)

        // grouped_by_op_ch.view()
        // Make a merged dataframe that contains the EGAD results for each of the EGAD splits
        makeMeltedMergedDF(grouped_by_op_ch)
}

workflow run_on_all_Gtex { // This does not work currently but it is used to run EGAD on all GTEX across OPs
    take:
        bulk_merged_ch
        organism_parts_ch
        go_annotations_ch
        go_annot_gene_column_ch
        bootstrap_ch
        bulk_sim_size_ch

    main:


        // Turn h5ad to csv
        convertGtexToCSV(bulk_merged_ch)

        // Sample N Bulk samples with replacement for each organism part. The output of this will be dataframes of N samples for each Organism Part
        sampleWithReplacement(bootstrap_ch, splitMerged.out, bulk_sim_size_ch)
        sampleWithReplacement.out.view()
        // Make tuple containing the Bulk datasets as well as the required GO Annotations
        // split_gtex_tuple_ch = sampleWithReplacement.out.combine(go_annotations_ch)
        // gtex_tuple_ch = convertGtexToCSV.out.combine(go_annotations_ch)


        // // Combine the split channel and gtex channel
        // all_bulk_for_EGAD_ch = Channel.of().concat(split_gtex_tuple_ch,gtex_tuple_ch)

        // // Then run EGAD on the splits and the whole EGAD
        // sim_bulk_EGAD(all_bulk_for_EGAD_ch, go_annot_gene_column_ch)

        // // Make a merged dataframe that contains the EGAD results for each of the EGAD splits
        // makeMeltedMergedDF(sim_bulk_EGAD.out.collect())
}

workflow bootstrap {
    
    main:
        def bootstrap_num = 0
        def bootstrap_lists = []    
        while (bootstrap_num++ < params.num_bootstrap) {
            bootstrap_lists << bootstrap_num
        }
        bootstrap_ch = Channel.fromList(bootstrap_lists)

        bulk_pipe(bulk_merged_ch,
            organism_parts_ch,
            go_annotations_ch,
            go_annot_gene_column_ch,
            bootstrap_ch,
            bulk_sim_size_ch)

}