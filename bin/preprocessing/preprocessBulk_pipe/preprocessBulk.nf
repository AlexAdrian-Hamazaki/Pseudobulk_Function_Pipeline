
// Import Other pipelines

//include {downloadPC; createPCMap} from "./preprocess_proteins.nf"

// Processes

GTEX_ch = Channel.fromPath(params.GTEX)
GTEX_meta_ch = Channel.fromPath(params.GTEX_meta)

process turnGTExtoAD {
    conda params.python3_9 
    //publishDir "${params.stagingDir}/Gtex",  mode: 'symlink'

  
    input:
        path GTEX_ch
        path GTEX_meta_ch
    output:
        path "${GTEX_ch.getSimpleName()}.h5ad.gz"
    shell:
    '''
    bulkScripts_gTEXToAnndata.py !{GTEX_ch} !{GTEX_meta_ch} !{GTEX_ch.getSimpleName()}
    '''

}

process filterPCBySymbol {
    conda params.python3_9 

    publishDir 'data/bulk', mode: 'copy'

    input:
        val pc
        path dataset
    output:
        path 'GtexMergedPC.h5ad'
    shell:
    '''
    bulkScripts_turnENSG.py !{pc} !{dataset}
    '''
}


process mergeDatasets2 {
    conda params.python3_9
    publishDir 'data/bulk', mode: 'symlink'

    input:
        path '*ENSGdf.h5ad' // list of paths to h5ad bulk datas 
    output:
        path 'bulkMerged.h5ad'

    script:
    """
    bulkScripts_mergeH5AD.py ${'*ENSGdf.h5ad'} 
    """
}

process downsampleBulk{
    conda '/home/aadrian/anaconda3/envs/PseudoPipelineEnv/'

    publishDir "data/bulk/downsampled"

    input:
        path mergedh5ad
        each k // value you want to downsample too
        val bootstrap // How many downsample bootstraps you want
    
    output:
        path "K${k}_boot*_bulk.h5ad"
        path "csvs/K${k}_boot*_bulk.csv"


    shell:
        """
        bulkScripts_downsampleBulk.py !{mergedh5ad} !{k} !{bootstrap}
        """
}


process makeTissueProfile {
    conda '/home/aadrian/anaconda3/envs/PseudoPipelineEnv/'

    publishDir "data/bulk/"

    input:
        path mergedh5ad
    output: 
        path "organismPartProfile.h5ad"
    shell:
    '''
    bulkScripts_makeTissueProfile.py !{mergedh5ad}
    '''
}


// Channels

//k_ch = Channel.fromList(params.k)  // Tuple of K's to subsample down to
bootstrap_ch = Channel.of(params.bootstrap)
pc = Channel.of(params.pc)

// Workflow
workflow processBulk_pipe {

    // take:
    //     pc // master PC list from downloadPCGenes.nf pipeline
    //     datasets  // bulk GTEX dataset


    main:
    // Make bulk matrix into h5ad
    turnGTExtoAD(GTEX_ch, GTEX_meta_ch)

    // Filter by PC looking at HGNC May have to add a parameter later saying if the genes are in ENSG or HGNC symbol
    filterPCBySymbol(pc, turnGTExtoAD.out)
    
    // Combine all Bulk Dataframes into 1 h5ad dataframe
    //mergeDatasets2(filterPCBySymbol.out.collect()) // publish this a "large"

    // Process stats about this
    // How many N we have for each tissue.


    /// This should go in its own script later #############

    // Downsample Bulk Samples to config input Ks
    //downsampleBulk(filterPCBySymbol.out, params.k, bootstrap_ch) // Publish this in data


    //TODO: Once we truly have different experiments. Perform batch correction/normalization before collapsing

    // collapse into row per tissue cell type // publish this
    makeTissueProfile(filterPCBySymbol.out)

}