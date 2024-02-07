#!/usr/bin/env python

import sys
import pandas as pd

def main():
	go_term_path = str(sys.argv[1])
	go_term_annotations_path = str(sys.argv[2])
	ontology_name = str(sys.argv[3])

	annotations = pd.read_csv(go_term_annotations_path, header=None)
	terms = pd.read_csv(go_term_path, sep = "\t")
	print(terms.shape)
	# get descriptions for all the GO terms we will be using
	terms_to_keep = terms[terms.id.isin(annotations.loc[:,0])]
	print(terms_to_keep.shape)
	terms_to_keep.to_csv(f"{ontology_name}_final_qc_terms.csv", index=False)
 
if __name__ == "__main__":
	main()