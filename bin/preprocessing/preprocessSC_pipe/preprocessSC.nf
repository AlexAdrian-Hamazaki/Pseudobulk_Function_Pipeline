

process transformToH5AD {
    conda params.python3_9 

    publishDir "data/h5ad_datasets/raw", mode: 'copy'

    input: 
        each sc_datasets
        path sc_cluster_meta
        path genes_meta
    output:
        path "${sc_datasets.getSimpleName()}_sc_with_metadata.h5ad"
    shell:
    '''
    transformtoH5AD.py !{sc_datasets} !{sc_datasets.getSimpleName()} !{sc_cluster_meta} !{genes_meta}
    '''
}

process generateQCFigs {
    conda params.python3_9 

    publishDir "data/scqc/${h5ad_dataset.getSimpleName()}", mode: 'move'

    input:
        path h5ad_dataset
    output:
        path "figures/highest_expr_genes_${h5ad_dataset.getSimpleName()}.png"
        path "figures/violin_${h5ad_dataset.getSimpleName()}.png"
        path "figures/scatter_${h5ad_dataset.getSimpleName()}_mito.png"
        path "figures/scatter_${h5ad_dataset.getSimpleName()}_duplet.png"

    shell:
    '''
    performSCQC.py !{h5ad_dataset} !{h5ad_dataset.getSimpleName()}
    '''
}

process performSCQC {
    conda params.python3_9 

    publishDir "data/h5ad_datasets/processed/cpm", mode: 'copy'

    input:
        each h5ad_dataset
        path cutoffs
    output:
        path "${h5ad_dataset.getSimpleName()}_cpm.h5ad"
        path "${h5ad_dataset.getSimpleName()}.h5ad"


    shell:
    '''
    performSCQC2.py !{h5ad_dataset} !{h5ad_dataset.getSimpleName()} !{cutoffs}
    '''
}

process filterPCGenes {
    conda params.python3_9 

    publishDir "data/h5ad_datasets/processed/pc", mode: 'copy'

    input:
        each h5ad_dataset
        path pc_map
    output:
        path "${h5ad_dataset.getSimpleName()}_pc.h5ad"

    shell:
    '''
    filterPCGenes.py !{h5ad_dataset} !{h5ad_dataset.getSimpleName()} !{pc_map}
    '''
}

process cluster{
    conda params.python3_9 

    publishDir "data/clusterFigs/${h5ad_dataset.getSimpleName()}", mode: 'move'

    input:
        path h5ad_dataset
    output:
        path "figures/filter_genes_dispersion_${h5ad_dataset.getSimpleName()}.png"
        path "figures/pca_variance_ratio_${h5ad_dataset.getSimpleName()}.png"
        path "figures/umap_${h5ad_dataset.getSimpleName()}.png"

    shell:
    '''
    cluster.py !{h5ad_dataset} !{h5ad_dataset.getSimpleName()}
    '''
}

process makeStats {
    conda params.python3_9 

    publishDir "data/h5ad_datasets/processed/withStats/unmerged", mode: 'copy'

    input:
        path h5ad_dataset
    output:
        path "${h5ad_dataset.getSimpleName()}_stats.h5ad"

    shell:
    '''
    makeStats.py !{h5ad_dataset} !{h5ad_dataset.getSimpleName()}
    '''
}

process makeStats2 {
    conda params.python3_9 

    publishDir "data/h5ad_datasets/processed/withStats", mode: 'copy'

    input:
        path h5ad_dataset
    output:
        path "${h5ad_dataset.getSimpleName()}_stats.h5ad"

    shell:
    '''
    makeStats.py !{h5ad_dataset} !{h5ad_dataset.getSimpleName()}
    '''
}

process combineH5AD {
    conda params.python3_9 

    publishDir "data/h5ad_datasets/", mode: 'copy'

    input:
        path "*" // This will accept all of the file inputs as a list of inputs
    output:
        path "merged_adata.h5ad"
    shell:
    '''
    mergeH5AD.py *  # this will pass all of the file inputs into sys.argv
    '''
}


process cluster2{
    conda params.python3_9 

    publishDir "data/clusterFigs/${h5ad_dataset.getSimpleName()}", mode: 'move'

    input:
        path h5ad_dataset
    output:
        path "figures/filter_genes_dispersion_${h5ad_dataset.getSimpleName()}.png"
        path "figures/pca_variance_ratio_${h5ad_dataset.getSimpleName()}.png"
        path "figures/umap_${h5ad_dataset.getSimpleName()}.png"

    shell:
    '''
    cluster.py !{h5ad_dataset} !{h5ad_dataset.getSimpleName()}
    '''
}



sc_datasets = Channel.fromPath(params.tissue_datasets, type:'dir')
sc_cluster_meta = Channel.fromPath(params.sc_cluster_meta)
genes_meta = Channel.fromPath(params.genes_meta)
curated_cutoffs = Channel.fromPath(params.curated_cutoffs)
pc_map = Channel.fromPath(params.pc_map)

 workflow generateQCFigs_pipe {

        // Turn all SC in to .h5ad objects
        transformToH5AD(sc_datasets, sc_cluster_meta, genes_meta)
    
        // Perform QC on SC dataset and output .h5ad object with good cells. This will happen for each dataset in sc datasets
        generateQCFigs(transformToH5AD.out) // Publish graphs for QC measures to some sort of data/qc. This would happen for pre- and post-processing

    emit:
        transformToH5AD.out
 }

workflow qcCells_pipe1 {
    
    take:
        transformToH5AD_out

    main:
        // filter for PC genes
        filterPCGenes(transformToH5AD_out, pc_map)

        // Actually remove the cells from the sc data based on thresholds identified in  processSC_pipe1 
        performSCQC(filterPCGenes.out, curated_cutoffs)

        // Make stats. NOTE: the STATS DF has logCPM as the .X
        makeStats(performSCQC.out[0])

        // Perform clustering to make sure cell types are clustering together within each dataset
        cluster(makeStats.out)


        // Combine all SC dataframes into h5ad object. Note. This will not have the stats. this is just PC filtered and logCPM taken
        combineH5AD(performSCQC.out[0].collect()) 

        // Make stats of the combined dataframe
        makeStats2(combineH5AD.out)

        // Perform clustering to see how all the cells clsuter across datasets as well.
        cluster2(makeStats2.out)
}

workflow scQC_pipe{
    main:
        generateQCFigs_pipe()
        qcCells_pipe1(generateQCFigs_pipe.out)

}