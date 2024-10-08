
// Processes

process processGaf { // Turn gaf file into csv

    publishDir "${params.PublishDir}/data/processing", mode: 'copy'
    conda params.python3_9


    input:
        path gaf_file
    output: 
        path "*.csv.gz"
    shell:
        '''
        processGO.py !{gaf_file}
        '''
}



process addENSG { // add protein data
    publishDir "${params.PublishDir}/data/processing", mode: 'copy'
    conda params.python3_9


    input:
        each path(go_annot)
        path pc_map

    output:
        path "*_withGeneData.csv.gz"
    shell:
    '''
    addENSG.py !{go_annot} !{go_annot.getSimpleName()} !{pc_map}
    '''
}

process QCGO { //QC the GO terms by getting only the ones with 20 or more genes
    publishDir "${params.PublishDir}/data/processing", mode: 'copy'
    conda params.python3_9


    input:
        path go_annot
    output:
        path "${go_annot.simpleName}_qc_annotations.csv" // The final output
        path "${go_annot.simpleName}_qc_goTerms.csv"// The final output
        path "${go_annot.simpleName}_no_duplicates.csv.gz"// The final output


        path "processed_${go_annot.simpleName}_distribution.png" // Distribution of how many genes are in each GO TErm
        path "unprocessed_${go_annot.simpleName}_distribution.png"
        // path "summary_${gocsv.simpleName}.csv"

    shell:
     '''
     qcGO.py !{go_annot} !{go_annot.simpleName}
     '''
}

process parse_obo {
    publishDir "${params.PublishDir}/data/processing", mode: 'copy'
    input:
        path go_obo_ch
    output:
        path "go_metadata.tsv"
    script:
    """
    getGOBasicInfo.R ${go_obo_ch}
    """
}

process combine_gos {
    // Also output a list of only the good GO terms that we are using, as well as their description
    publishDir "${params.PublishDir}/data/final", mode: 'copy'
    conda params.python3_9

    input:
        path go_metadata
        each path(go_annot)
    output:
        path "${go_annot.getSimpleName()}_final_terms_with_metadata.csv"
    script:
    """
    get_final_go_terms.py ${go_metadata} ${go_annot} ${go_annot.getSimpleName()}
    """
}

process cellTypeRelatedness_pipe {
    publishDir "${params.PublishDir}/data/cellTypeRelatedness", mode: 'copy'

    input:
        path obo
        path curated_regulations
        path good_terms
        path bad_terms
        path functions
    output:
        path "*_ctRelatedness.csv"
    shell:
    '''
    curateRegexAll.R !{obo} !{curated_regulations} !{good_terms} !{bad_terms} !{functions}
    '''
}

process remove_dependance {

    publishDir "${params.PublishDir}/data/final", mode: 'copy'

    input:
        path go_with_annotations
    output:
        path "${go_with_annotations.getSimpleName()}_no_dependance.csv"
        path "${go_with_annotations.getSimpleName()}_dependent_terms.csv"
        path "${go_with_annotations.getSimpleName()}_kept_terms.txt"
    script:
    """
    removeGOOverlap.py ${go_with_annotations} ${go_with_annotations.getSimpleName()}
    """

}

// Channels
gaf_ch = Channel.fromPath(params.gaf)
pc_map_ch = Channel.fromPath(params.pc_map)
go_obo = Channel.fromPath(params.GO_obo)

workflow go_annotations_pipe {
    // Process GAF file to proper format, save BP and MF
    processGaf(gaf_ch)

    // Add ENSG to the gene symbols in the GO ontologies
    addENSG(processGaf.out.flatten(), pc_map_ch)
    
    // Remove GO terms that don't have at least 20 genes in them. Get summary stats as well
    // QC to get only GO Terms that have 20+ Genes
    QCGO(addENSG.out)

    // Parse OBO file to just get GO name and description
    parse_obo(go_obo)

    
    // Remove GO Terms that are dependent, keepiing only one of the dependent GO terms
    remove_dependance(QCGO.out[0])
    remove_dependance.out[0].view()
    // Get the NAmes and description of the GO terms we want to use
    combine_gos(parse_obo.out, remove_dependance.out[0])
}


bad_terms_ch = Channel.fromPath(params.bad_terms)
good_terms_ch = Channel.fromPath(params.good_terms)
curated_regulations_ch = Channel.fromPath(params.curated_regulations)
go_obo_ch = Channel.fromPath(params.goobo)
functions = Channel.fromPath("/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/preprocessing/preprocessGO_pipe/bin/functions.R")

workflow cell_type_relatedness_pipe {

    // Run the ontology propagation script
    onto_regex(go_obo_ch, curated_regulations_ch, good_terms_ch, bad_terms_ch, functions)
}