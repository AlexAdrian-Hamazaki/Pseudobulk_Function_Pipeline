#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np


def main():
    # Load df of Go terms with annotations
    go_annot_path = sys.argv[1]
    go_name = sys.argv[2]
    go_annot = pd.read_csv(go_annot_path)
    
    # Count the genes annotated with each GO term
    bp_grouped = countGenesInGO(go_annot)
    
    # Calculate the overlap of genes between GO Terms
    bp_allPercDFs = createOverlapPercDF(bp_grouped)
    
    # Get only the dependent GO terms (GO terms with 70% gene overlap both ways)
    bp_dependent = getDependentGOs(bp_allPercDFs, 0.7)
    #save the dependent df with the % overlap values
    bp_dependent.to_csv(f"{go_name}_dependent_terms.csv", index=False)
    
    # Make sure that there are no GO terms that are dependent on 2 or more other GO terms. they are all only dependent on 1 go term
    assert checkDependentGOs_group3allterms(bp_dependent=bp_dependent)
    
    # Collapse the GO terms that are dependent by randomly selecting one of them
    lo_kept_go_terms, lo_removed_go_terms = getGOTermsToRemove(bp_dependent=bp_dependent)
    # save the kept go terms as a txt file
    save_txt(my_list=lo_kept_go_terms, go_name=go_name)

    # Get the GO terms to keep. This is the annotated GO terms
    final_go_annotated = removeDependentGOs(go_annotated = go_annot, lo_removed_go_terms = lo_removed_go_terms)
    # Save this as the final go annotated df
    final_go_annotated.to_csv(f"{go_name}_no_dependance.csv", index= False)


def save_txt(my_list, go_name):
    # save a new line delimited text file from a 
    # Write the list to a new-line delimited text file
    with open(f"{go_name}_kept_terms.txt", 'w') as file:
        for item in my_list:
            file.write(f"{item}\n")


def getOnlyGOTermsToKeep(final_go_annotated:pd.DataFrame) -> list:
    """Get a list of the GO Ids we want to keep

    Args:
        final_go_annotated (pd.DataFrame): go annotated df

    Returns:
        list: lo Go terms to keep
    """
    lo_go_ids = final_go_annotated.loc[:,"GO ID"].unique()
    return lo_go_ids
    
# now that I know what GO terms to keep and remove. I can remove the GOs that I want to remove
def removeDependentGOs(go_annotated:pd.DataFrame, lo_removed_go_terms:list)->pd.DataFrame:
    """Remove dependent GOs

    Args:
        go_annotated (pd.DataFrame): full list of go's with anotations
        lo_removed_go_terms (list): list of GO ids to remove

    Returns:
        pd.DataFrame: independent GOs
    """
    return go_annotated[~go_annotated.isin(lo_removed_go_terms)]
    
# Collapse the GO terms that are dependent by randomly selecting one of them
def getGOTermsToRemove(bp_dependent:pd.DataFrame, lo_kept_go_terms:list=[],  lo_removed_go_terms:list=[]) -> list:
    """For a df of GO terms that have been checked to ONLY have GO terms hat are dependent on 1 more GO term,
    
    randomly selects one of those dependent GO terms to keep 

    Args:
        bp_dependent (pd.DataFrame): _description_

    Returns:
        list: _description_
    """
    # shuffle
    bp_dependent = bp_dependent.sample(frac=1, random_state=100)
    
    # iterate through rows, keep target go term. and remove 'other' go term from dataframe
    
    top_row = bp_dependent.iloc[0,:]
    target = top_row['test_term']
    other = top_row['other_term']
    
    # add target to list of kept go terms
    lo_kept_go_terms.append(target)
    lo_removed_go_terms.append(other)
    
    # remove the 'other' call from the df so it is not also kept
    bp_dependent = bp_dependent[bp_dependent['test_term'] != other]
    bp_dependent = bp_dependent[bp_dependent['test_term'] != target]
    
    if bp_dependent.empty:
        return lo_kept_go_terms, lo_removed_go_terms
    else:
        return getGOTermsToRemove(bp_dependent, lo_kept_go_terms, lo_removed_go_terms)
    


def checkDependentGOs_group3allterms(bp_dependent:pd.DataFrame)->bool:
    """Check our dataframe of dependent GO terms for situations where there is actually 3 or more GO terms that are all dependent with eachother.
    
    If there are any, then raise an error. Currently I have not found this to be the case in my subset of GO terms, so I will not implement a way to deal with this problem because it does not exist

    Args:
        bp_dependent (pd.DataFrame): dataframe of dependent GO terms
    Returns:
        True if we find no sets of go terms of 3 or greater that are dependent

    Raises: ValueError: If we find that there may be a set of GO terms of length 3 or greater that are dependent on eachother
    """

    for test_term in bp_dependent.test_term:
        test_term_calls = bp_dependent[bp_dependent.loc[:,'test_term']==test_term]
        if test_term_calls.shape[0] == 1:
            continue
        else:
            raise ValueError(f'Found a set of GO terms of length 3> that are dependent {test_term_calls}')
    return True

    
def getDependentGOs(df_overlap: pd.DataFrame, overlapPerc: float):
    """ Return a dataframe containing only GO Terms to other GO Terms that overlap

    Args:
        df_merged (pd.DataFrame): _description_
        overlapPerc (float): How much percent you need to overlap
        
    Returns:
        filtered_df (pd.DataFrame): A pandas dataframe containing GO terms that have genes that overlap significantly
    """
    filtered_df = df_overlap[(df_overlap['perc1'] >= overlapPerc) & (df_overlap['perc2'] >= overlapPerc)]
    
    return filtered_df
    
def createOverlapPercDF(grouped_df: pd.DataFrame):
    """For one grouped df, calculate how many genes in each GO term overlap with the genes in each other GO term
    

    Args:
        df_grouped (pd.DataFrame): mf or bp that is grouped and sorted by ascending value for group size
        
    Returns:
        allPercDCs (pd.DataFrame): All of the GO to other GO term percent overlaps
        
    """
    # A list of percent overlap dataframes
    loPercDFs =[]
    
    for index, row in grouped_df.iterrows():
        
        # testID is the GO term we are testing overlap FOR
        testID = row["GO ID"]

        
        # testArray will be the genes affiliated with our testID GO Term
        testArray = row['loGenes'] #
        
        # otherRows is a dataframe of the other GO terms that are NOT our current test ID GO term
        otherRows = grouped_df.loc[grouped_df.index != index] # The other arrays
        
        # There are the intersecting GO Terms
        percDF = getIntersectingArraysFromOne(testID, testArray, otherRows)
        
        # append to list of perc DFs
        loPercDFs.append(percDF)
        
    # Row bind all of the percent dataframes
    allPercDFs = pd.concat(loPercDFs, axis = 0)
    
    return allPercDFs

def countGenesInGO(df: pd.DataFrame):
    """Count the number of genes in each GO ID

    Args:
        df (pd.DataFrame): bp or mf dataframe
    """
    
    grouped_df = df.groupby('GO ID').agg({'DB_Object_Symbol': list, 'GO ID': 'size'})
    
    grouped_df = grouped_df.rename(columns={'GO ID': 'Group_Size', 'DB_Object_Symbol': 'loGenes'})    
    
    grouped_df = grouped_df.sort_values('Group_Size')
    

    grouped_df = grouped_df.reset_index()

    grouped_df.reset_index(drop=True)
    return grouped_df


def getIntersectingArraysFromOne(testGOTerm: str, loGenesTest: np.array, otherRows: pd.DataFrame):
    
    """
    Get the arrays from a list that have an intersection greater than a given percentage with a test array.
    We are testing the intersection between testArray and the loGenes column of the otherRows dataframe

    Parameters:
        testGOTerm (str): The GO term that you are testing overlap for
        loGenesTest (np.ndarray[str]): List of genes that are afiliated with your testGOTerm
        otherRows: pd.DataFrame: Pandas dataframe where index is GO term, and loGenes is an array of genes contained by that go term

    Returns:
        df_percent_overlap (pd.DataFrame): Df where you have testTerm to other GO term percent overlaps
        
        
    """
    
    # A list of percent overlaps
    lo_percent_overlap = []
    

    for index, row in otherRows.iterrows():
        loGenesArray = row["loGenes"]
        assayGOTerm = row["GO ID"]
        
        ar1, ar2 = IntersectArrays(loGenesTest, loGenesArray) # Gets the percentage of itneractions between all the arrays
        lo_percent_overlap.append([testGOTerm, assayGOTerm, ar1, ar2]) # index will be the array GO Term
        
    
    df_percent_overlap = pd.DataFrame(data = lo_percent_overlap)
    df_percent_overlap.columns = ['test_term', 'other_term', 'perc1', 'perc2']
    
    return df_percent_overlap

def IntersectArrays(array1, array2):
    """
    Check if two numpy arrays have an intersection greater than a given percentage.

    Parameters:
        array1 (numpy.ndarray): The first numpy array.
        array2 (numpy.ndarray): The second numpy array.
        perc (int): The minimum intersection percentage required for the arrays to be considered intersecting.

    Returns:
        intersection_percentage_1, intersection_percentage_2: Two fractions. The first indicates how much % overlap there \
            when array 1 is in the denominator. The second indicates fraction whenarray 2 is in denominator
    """
    intersection = np.intersect1d(array1, array2)
    intersection_percentage_1 = len(intersection) / len(array1)
    intersection_percentage_2 = len(intersection) / len(array2)

    return intersection_percentage_1, intersection_percentage_2




if __name__ == "__main__":
    main()