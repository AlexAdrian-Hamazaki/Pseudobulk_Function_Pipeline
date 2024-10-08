#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys

def filterPC(PCGenes:pd.DataFrame):

    print("Filtering for high conf")
    # Filter for High Confidence Proteins
    PCGenes = PCGenes[PCGenes.loc[:,"Protein existence"] == "Evidence at protein level"]
    
    return PCGenes




def main():

    PCPath = sys.argv[1]
    print(PCPath)



    print("Opening uniprot")
    PCGenes = pd.read_csv(f'{PCPath}', sep = "\t", compression='gzip')
    print("Uniprot Opened")
    print(PCGenes.shape)
    
    print(PCGenes.head())
    
    filtered = filterPC(PCGenes=PCGenes)

    filtered.to_csv('filteredPC.tsv.gz', sep = "\t", compression = "gzip")
    

if __name__ == "__main__":
    print("running")
    main()
    