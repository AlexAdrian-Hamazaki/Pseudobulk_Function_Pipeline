
process calc_ct_profiles {
    publishDir "${params.publish}"
    conda params.python3_9
    queueSize = 1
    pollInterval = '30 sec'
    memory = '16 GB'

    input:
        path adata
    output:
        path "all_human_CT_profiles.csv"
    
    script:
    """
    makeCellTypeProfiles.py ${adata}
    """
}

process remove_CTs_iteratively {
    publishDir "${params.publish}/removals", mode: 'copy'
    conda params.python3_9
    queueSize = 1
    pollInterval = '30 sec'
    memory = '16 GB'


    input:
        path all_human_CT_profiles
    output:
        path "*_removed_*"
        
    
    script:
    """
    remove1ProfileSystematically.py ${all_human_CT_profiles}
    """
}

process EGAD {
    executor "local"
    // maxForks 20

    publishDir "${params.publish}/EGAD", mode: 'copy'

    input:
        path expression_matrix
        each path(go_annotations_ch)
    output:
        path "*_EGAD.csv"
    script:
    """
    Egad.R ${expression_matrix} ${expression_matrix.getBaseName()} ${go_annotations_ch} ${go_annotations_ch.getSimpleName()}  ${expression_matrix.getSimpleName()}
    """
}

adata_ch = Channel.fromPath(params.path_to_sc_adata)
go_annotations_ch = Channel.fromPath(params.go_annotations)

workflow {
	// Purpose: Perform the CT removal experiment where I systematically remove each CT Profile and quantify what GO Terms decrease the most   // // First, open the adata object with all of the cells and  calculate the CT profiles for each cell. Save that
    calc_ct_profiles(adata_ch)
	
	// Create X dataframes, each dataframe has 1 CT Profile removed. X will be equal to the numberof unique CT Profiles we have
	remove_CTs_iteratively(calc_ct_profiles.out)
	// remove_CTs_iteratively.out.flatten().view()
	// Run EGAD on each of the X dataframes. A well as the ENTIRE dataframe
	dfs = remove_CTs_iteratively.out.flatten().concat(calc_ct_profiles.out)
    // dfs.view()
	EGAD(dfs, go_annotations_ch)
	// Perform Statistical Comparision to see what GO Terms decrease the most for each CT Profile
	// stats(EGAD.out, calc_ct_profiles.out)
	// Perhaps perform an assessment of what GO terms systematically decrease


}