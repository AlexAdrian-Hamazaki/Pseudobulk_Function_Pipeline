#!/usr/bin/env python3

import pandas as pd
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import csv
import numpy as np

def goHist(ser:pd.Series, GO_name:str):
    """Create a histogram showing the distribution of how many genes are in each GO term
    
    There is also a red vertical line a x = 19 to show the terms being removed

    Args:
        ser (pd.Series): named series. names are go ids, entries are ints showing how many genes arein that go term
        GO_name (str): the name of what GO ontology we are curating. for naming purposes only
    """
    # Clear the previous plot
    plt.clf()
    # Set the style to clean
    sns.set_style("white")
    # Create the histogram using distplot
    sns.histplot(ser, kde=False, color='steelblue', binwidth=5)
    
    # Draw a vertical red line at x = 20
    plt.axvline(x=19, color='red', linestyle = 'dashed')
    
    # Draw a vertical red line at x = 501
    plt.axvline(x=501, color='red', linestyle = 'dashed')

    # Set the labels for x-axis and y-axis
    plt.xlabel('GO Terms')
    plt.ylabel('Counts')

    # Set the title for the plot
    plt.title('Number of Genes Affiliated with GO Term')
    
        
    # Set the y-axis limit to 100
    
    plt.xlim(0, 650)

    # Save
    plt.savefig(f'{GO_name}_distribution.png')
    

def getSummaryStats(groupedGO:pd.Series, filteredBigGOs:pd.DataFrame, bigGOs:pd.Series, GO_name:str):
    """
    Get summary stats for the GO filtering
    
    Saves to a file called "summary_${gocsv.simpleName}.csv"

    Args:
        groupedGO (pd.Series): a grouped series of GO terms
        filteredBigGOs (pd.DataFrame): dataframe containing only the GO terms with >=20 genes
        bigGOs (pd.Series): A series of only GO terms with >=20 genes. names are go terms. values are number of genes
        GO_name (str): name of go namespace. For saving file purposes
    """
    with open(f'summary_{GO_name}.csv', mode = "w", newline='') as file:
        
        file.write(f'SUMMARY STATS FOR {GO_name}\n')

        
        # How many GO terms there were initially
        file.write(f'There were initially {len(groupedGO)} GO Terms \n')
        
        # How many GO terms were removed
        file.write(f'{len(groupedGO) - len(bigGOs)} GO Terms were removed because they did not contain >= 20 genes & <= 500 genes\n')

            
        # What are the 5 GO terms with the most Genes and how many genes are there
        file.write("The GO terms with the most Genes within them:\n")
        file.write(f'{takeTopFive(bigGOs)}')
        
        
def takeTopFive(ser:pd.Series):
    """returns the 5 largest values from a named series of ints

    Args:
        ser (pd.Series): series of ints. Names are go ids, values are how many genes are in that go id

    Returns:
        pd.series: series of top 5 go terms
    """
    # Sort the Series in descending order
    sorted_series = ser.sort_values(ascending=False)

    # Take the top 5 elements
    top_5_elements = sorted_series.head(5)
    
    return top_5_elements
    
    

def main():
    GOpath = sys.argv[1]
    GO_name = sys.argv[2]
    
    
    # Load GO
    GO = pd.read_csv(GOpath)
    
    # process GOs by removing duplicates
    GO = prelim_processing(GO)

    # Group the GO Terms by their GO IDs and get the size of each GO term
    groupedGO = GO.groupby("GO ID").size()
    # Get the GO terms that contain more than 20 GO terms or less than 500
    bigGOs_ser = groupedGO[(groupedGO >= 20) & (groupedGO <= 500)]
    
    # Graph the size of each GO Term just for QC
    goHist(ser = groupedGO, GO_name=f'unprocessed_{GO_name}')
    goHist(ser = bigGOs_ser, GO_name= f'processed_{GO_name}')
    
    # Filter the list of GO term annotations to get the GO terms and annotations for only the QCd GO Terms
    filteredBigGOs = GO[GO['GO ID'].isin(bigGOs_ser.index.to_list())]
    
    # Save the list of GO term to Gene Relationships
    filteredBigGOs.to_csv(f'{GO_name}_qc_annotations.csv', index=False)
    
    #Get the list of ALL the GO TERMS that match our criteria. Just the GO terms, not the annotations. Save this list
    np.savetxt(f'{GO_name}_qc_goTerms.csv', bigGOs_ser.index.to_list(), fmt='%s', delimiter='\n')

    
    
    # getSummaryStats(groupedGO=groupedGO,
    #                 filteredBigGOs=filteredBigGOs,
    #                 bigGOs=bigGOs_ser,
    #                 GO_name=GO_name)
    
def prelim_processing(GO:pd.DataFrame):
    """Remove bad duplicate rows

    Args:
        GO (pd.DataFrame): _description_
    """
    # remove row duplicates
    GO.drop_duplicates(inplace=True)
    
    # Remove instances where rows are duplicated because of bad GO annotations wioth the protein information
    # upon doing this, we are commiting to only using hgnc symbols
    lo_go_term_dfs = [] # list to hold go term annotations for each go term
    # group GO by GO term
    go_grouped = GO.groupby('GO ID')
    for go_term, annot_df in go_grouped:
        annot_df = annot_df.drop_duplicates(subset='hgnc_symbol', keep='first')
        lo_go_term_dfs.append(annot_df)
    annot_df = pd.concat(lo_go_term_dfs, axis = 0)
    return annot_df

if __name__ == "__main__":
    main()