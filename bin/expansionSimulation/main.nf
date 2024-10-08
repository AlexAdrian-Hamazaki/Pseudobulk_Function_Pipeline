
process calc_ct_profiles {
    publishDir "${params.publish}"
    conda params.python3_9

    input:
        path adata
    output:
        path "all_human_CT_profiles.csv"
    
    script:
    """
    makeCellTypeProfiles.py ${adata}
    """
}

process subsample_ct_profiles { // Subsample cell type profiles. Subsample different amounts of them
    publishDir "${params.publish}/ct_profiles"
    conda params.python3_9

    input:
        path all_human_CT_profiles_csv
        val ct_profile_step
    output:
        path "*_ct_pros.csv"
    script:
    """
    downsampleCTProfiles.py ${all_human_CT_profiles_csv}  ${ct_profile_step}
    """
}

process EGAD {
    executor "local"

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
min_ct_profiles = Channel.of(params.min_ct_profile)
step = Channel.of(params.ct_profile_step)
go_annotations_ch = Channel.fromPath(params.go_annotations)

workflow {
    // // First, open the adata object with all of the cells and  calculate the CT profiles for each cell. Save that
    calc_ct_profiles(adata_ch)
    // calc_ct_profiles.out.view()

    // Save a series of CSVs. In each CSV, sample cell type profiles. Start with 5 randomly smaple ct profiles. Add 5 each time, and run EGAD. Do this for control, but also do it for composition effects where we shuffle the variation that each CT profiles makes
    subsample_ct_profiles(calc_ct_profiles.out, step)

    // Run EGAD on each of those CSVs
    // subsample_ct_profiles.out.view()
    // subsample_ct_profiles.out.flatten().view()
    EGAD(subsample_ct_profiles.out.flatten(), go_annotations_ch)

    // Create a melted dataframe with all the info
    // melt_EGAD(EGAD.out)
}