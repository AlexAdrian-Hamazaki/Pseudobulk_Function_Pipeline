#!/usr/bin/env python

import pandas as pd
import sys

# 
# Column names of the .gaf file can be found here: http://geneontology.org/docs/go-annotation-file-gaf-format-2.2/#aspect-column-9
# 
# Data downloaded from here: http://current.geneontology.org/products/pages/downloads.html



def main():
    """Process a go GAF file
    
    Saves two files. each has 4 columns
    DB_Object_Symbol: Gene Symbol
    GO ID : Go ID
    Aspect : MF, BP, CF
    DB Object Name: Gene Name
    
    
    bpProcessed.csv: A file with only the biological process GO relationships
    mfProcessed.csv: A file with only the  MFrelationships
    """
    GO_path = sys.argv[1]
    
    GO = pd.read_csv(GO_path, skiprows=41, sep = "\t", header = None)

    colnames = ['DB',
            'DB_Object_ID',
            'DB_Object_Symbol',
            'Qualifier',
            'GO ID',
            'DB:Reference (|DB:Reference)',
            'Evidence Code',
            'With (or) From',
            'Aspect',
            'DB Object Name',
            'DB Object Synonym (|Synonym)',
            'DB Object Type',
            'Taxon(|taxon)',
            'Date',
            'Assigned By',
            'Annotation Extension',
            'Gene Product Form ID']
    
    GO.columns = colnames
    
    GO = GO.loc[:,['DB_Object_Symbol', 'GO ID', 'Aspect', 'DB Object Name'	]]

    BP = GO[GO.loc[:, 'Aspect'] == "P"]
    
    MF = GO[GO.loc[:, 'Aspect'] == "F"]
    
    BP.to_csv('bp_annotations.csv.gz', index= False, compression='gzip')
    MF.to_csv('mf_annotations.csv.gz', index= False, compression='gzip')

    


if __name__ == "__main__":
    main()