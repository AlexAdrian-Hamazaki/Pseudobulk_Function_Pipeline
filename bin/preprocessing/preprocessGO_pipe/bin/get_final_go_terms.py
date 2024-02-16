#!/usr/bin/env python

import sys
import pandas as pd

def main():
    # Load annotated GO terms
    go_term_path = str(sys.argv[1])
    go_term_annotations_path = str(sys.argv[2])
    # Load parsed .obo file with all the GO terms
    ontology_name = str(sys.argv[3])
    
    # load annotations and ontology terms
    annotations = pd.read_csv(go_term_annotations_path)
    terms = pd.read_csv(go_term_path, sep = "\t")
    
	# get a LIST of just the GO terms that we have decided to keep
    lo_go_ids = getOnlyGOTermsToKeep(final_go_annotated = annotations)
    
    # get descriptions for all the GO terms we will be using
    terms_to_keep = terms[terms.id.isin(lo_go_ids)]
    
    terms_to_keep.to_csv(f"{ontology_name}_final_terms_with_metadata.csv", index=False)
 

def getOnlyGOTermsToKeep(final_go_annotated:pd.DataFrame) -> list:
    """Get a list of the GO Ids we want to keep

    Args:
        final_go_annotated (pd.DataFrame): go annotated df

    Returns:
        list: lo Go terms to keep
    """
    lo_go_ids = final_go_annotated.loc[:,"GO ID"].unique()
    return lo_go_ids

if __name__ == "__main__":
	main()