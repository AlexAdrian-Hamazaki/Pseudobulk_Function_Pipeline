#!/usr/bin/env python3

import sys
import anndata as ad
import pandas as pd
import json
import random
import scanpy as sc
import numpy as np
import math


def main():
    CTProfile_path = sys.argv[1]
    CTProfile_name = sys.argv[2]
    proportions_json_path = sys.argv[3]
    variance_factor = float(sys.argv[4])
    num_simulations = int(sys.argv[5])
    totalSampleSize = 500
    
    # open cell type profile database
    df = pd.read_csv(CTProfile_path,index_col=0)
    
    # Retreive the baseline proportions from the proportions_json
    dictBaselineProportion = getBaselineProportion(CTProfile_name=CTProfile_name, proportions_json_path=proportions_json_path)
    
    # Where each element in the list is the cell type proportions I'm going for
    # Get a list of dictionaries where each dictionary defines the proportions of cell types we need to sample to make our simulated bulk sample Keys are the cell type and values are the proportion of what our 
    # simulated bulk sample should look like
    
    loDictsCellTypeProportions = getDictsOfCellTypeProportions(num_simulations,
                                                               dictBaselineProportion,
                                                               totalSampleSize,
                                                               variance_factor)


    # Actually get the number of cells we need to sample for each simulated bulk dataset
    loDictsNumberOfCells = [getNumberToSample(dictCellTypeProportion, totalSampleSize=totalSampleSize) for dictCellTypeProportion in loDictsCellTypeProportions]
    
    # Turn all the dictionaries into pandas series
    loSerNumberOfCells = [pd.Series(dictNumberOfCells) for dictNumberOfCells in loDictsNumberOfCells]
    
        
    # For each of the cell type compositions, get a simulated single cell dataset that is uncollapsed
    loUncollapsedSimulatedBulkSamples = [simulateBulk(cellTypeComposition, df = df) for cellTypeComposition in loSerNumberOfCells]
    
    # For each simulated single cell dataset, collapse into one single simulted bulk sample
    loSimulatedBulkSamples = [collapseSimulation(loCellTypeProfiles) for loCellTypeProfiles in loUncollapsedSimulatedBulkSamples]
    
    
    # Make a simulated bulk dataset where each row is a simulated bulk sample
    df_simulatedBulkDataset = pd.concat(loSimulatedBulkSamples, axis = 1)
    
    
    # Save the bulk simulated dataset
    df_simulatedBulkDataset.to_csv(f'{CTProfile_name}.csv.gz', compression='gzip')
    
    # Also write the compositons for each of the simulated bulk dataset
    for i, composition in enumerate(loSerNumberOfCells):
        composition.to_csv(f'{CTProfile_name}_n_sim_{i}_profiles.csv', index=True)
    
    # Also save the initial simulated bulks that are not collapsed
    for i, uncollapsedSimulatedBulkSample in enumerate(loUncollapsedSimulatedBulkSamples):
        uncollapsedSimulatedBulkSample.to_csv(f'{CTProfile_name}_n_sim_{i}_profiles_uncollapsed.csv', index=True)
        
def getProportionToSample(dictBaselineProportion:dict,  variance_factor:float, totalSampleSize:int) -> pd.Series:
    """For one tissue, and given the number of unique cell types, and given a total sample size for our simulated bulk,
    identify how many cells of each cell type we want to sample"
    
    The default number of cells you want to sample will be based on the dictbaselineProportion, the default percent composition that I'm setting for each cell type

    Args:
        dictBaselineProportion:dict: A dict where keys are cell types, and values are the default percent composition of that cell type.
        totalSampleSize (int): The total number of cells we want in each simulated bulk dataset
        variance_factor (float): The amount of the compositional effect we want to simulate. There is actually a direct interpretation of this number though. For a simulated bulk of 4 cell types and for 100 samples, \
            a variance of 0.01 indicates that the composition of each cell type can differ by 1% of the total number of cells. \
                For example, with zero composition effect, each cell type would have 25,25,25,25 cells. With variance=0.01, \
                    the number of cells in each cell type can vary by 100*0.01 = 1 cell. So each cell type could have
                    24-26 cells in it

    Returns:
        pd.Series: Named series, has the proportions
    """
    
    
    dict_ProportionToSample = {}  # dict to store the number of cells for each cell type. Keys are the cell type, and values are the number of cells of that cell type to sample

    # Get all the keys (cell types)
    loKeys = list(dictBaselineProportion.keys())
    
    # Loop through each cell type
    for key in loKeys:

        proportion_to_sample = getProportionCellsToSample(baseline_cell_proportion=dictBaselineProportion[key], variance_factor=variance_factor)
        
        
        # Add this number to our dictionary which contains the number of cels to sample for this simulated bulk sample
        dict_ProportionToSample [key] = proportion_to_sample
        
    # Remove proportions that are Negative (done before scaling to 1). Negatives will be reset to 0
    dict_ProportionToSample = removeNegatives(dict_ProportionToSample)
    
    # When the variance is relaly high, sometimes NO cells will be sampled from the method above. In this case, just choose to only sample 1 cell to the total smaple size
    if sum(dict_ProportionToSample.values()) == 0:
        random_key = np.random.choice(list(dict_ProportionToSample.keys()))
        dict_ProportionToSample[random_key] = totalSampleSize
    
       
    # Scale the proportions to summarize to 1
    scaled_dict_ProportionsToSample = scaleDictTo1(dict_ProportionToSample)
    
    # round based on the number of sig figs we have in totalSampleSize
    scaled_dict_ProportionsToSample = roundToSigFigs(scaled_dict_ProportionsToSample=scaled_dict_ProportionsToSample, totalSampleSize=totalSampleSize)
    
    return  scaled_dict_ProportionsToSample

def removeNegatives(dict_ProportionToSample):
    
    """Sets negative values to for a dictionary
    """
    
    for key in list(dict_ProportionToSample.keys()):
        if dict_ProportionToSample[key] <0:
            dict_ProportionToSample[key] = 0
                
    return dict_ProportionToSample

def getProportionCellsToSample(baseline_cell_proportion:float, variance_factor:float):

    proportion_to_sample = np.random.normal(baseline_cell_proportion, (float(variance_factor) * float(baseline_cell_proportion)), size = 1)
    return proportion_to_sample[0] # Return first element because we want a float not a np array

def roundToSigFigs(scaled_dict_ProportionsToSample:dict, totalSampleSize:int):
    
    for key in list(scaled_dict_ProportionsToSample.keys()):
        scaled_dict_ProportionsToSample[key] = round(scaled_dict_ProportionsToSample[key], len(str(totalSampleSize)))
        
    return scaled_dict_ProportionsToSample

def scaleDictTo1(dict_ProportionToSample:dict):
    """Scale dictionary values such that the sum of the dict is equal to 1

    Args:
        dict_ProportionToSample (dict): Contaisn the proportions of cel ltypes we need to sample
        
    Returns:
        scaled_dict(dict): the scaled dict where the proportions sum to 1
    """
    total_sum = sum(dict_ProportionToSample.values())

    scaled_dict = {key: value/total_sum for key,value in dict_ProportionToSample.items()}
    
    return scaled_dict


def getShuffledKeys(dict):
    """Randomly shuffle keys and return them as a list
    

    Args:
        dict (_type_): _description_
    """
    
    # Get the keys view and convert it to a list
    keys_list = list(dict.keys())
    
    # Shuffle the keys list randomly
    random.shuffle(keys_list)
    
    
    return keys_list


        
def getBaselineProportion(CTProfile_name:str, proportions_json_path:json) -> dict:
    """Get the list of proper cell type baseline proportions depending on if we are using brain or pbmc

    Args:
        CTProfile_name (str): name of the dataframe of cell type profiles
        proportions_json (json): json that contains baseline cell type profiles for brian or pbmc

    Returns:
        dict: A dict where keys are cell types, and values are their expected baseline proportions
    """
    
    # open json file
    with open(proportions_json_path, 'r') as json_file:
        dict_propostions_json = json.load(json_file)
    
    
    # Get the correct key
    for key in dict_propostions_json.keys():
        if key == CTProfile_name:
            correct_key = key
    # If no key was found,raise error
    if not correct_key:
        raise ValueError("Key not found")
    
    # Return the proportions of cell types based on our key (which is a tissue type)
    return dict_propostions_json[correct_key]

    
def collapseSimulation(df: pd.DataFrame) -> pd.Series:
    """Each dataframe is a bunch of cell type profiles.
    
    This function just sums up all othe expression and returns a pandas series which represents a simulated bulk sample

    Args:
        df (pd.DataFrame): _description_

    Returns:
        pd.Series: A bulk simulated sample
    """
    
    return df.sum(axis = 0)

def getNumberToSample(dictCellTypeProportion:dict, totalSampleSize:int) -> dict:
    """Given a dictionary of cell type proportions and a total sample size, get the total number of cells to subsmaple

    Args:
        dictCellTypeProportion (dict): _description_

    Returns:
        dict: values are total number of cells
    """
    dictNumberToSample = {}
    
    for key in dictCellTypeProportion:
        dictNumberToSample[key] = int(dictCellTypeProportion[key] * totalSampleSize)
        
    # If the sum of cells the dictionary is not equal to the totalSampleSize then re-scale
    if int(sum(dictNumberToSample.values())) != int(totalSampleSize):
        dictNumberToSample = fixRoundingProblem(dictNumberToSample, totalSampleSize)
        
    return dictNumberToSample

def fixRoundingProblem(dictNumberToSample, totalSampleSize):
    """The values in this dict should add upto the totalSampleSize. If they don't then fix the rounding problems by adding or subtracting cells

    Args:
        dictNumberToSample (_type_): _description_
    """
    # Calculate the remaining difference between the dictionary and the totals ample size
    remaining_difference = totalSampleSize - sum(dictNumberToSample.values())
    
    # If the remaining difference is 0. Then exit the recursion by returning the final dict
    if remaining_difference == 0:
        return dictNumberToSample
    
    # If the remaining difference is positive, then we need to add cells randomly to the dictionary
    elif remaining_difference > 0:
        # randomly get a key
        key = np.random.choice(list(dictNumberToSample.keys()))
    
        # Add one cell to the dictionary
        dictNumberToSample[key] +=1
        
        return fixRoundingProblem(dictNumberToSample, totalSampleSize)
        
    # If the remining difference is negative, then that means we need to randomly remove a cell
    elif remaining_difference < 0:
        # randomly get a key
        key = np.random.choice(list(dictNumberToSample.keys()))
    
        # Add one cell to the dictionary
        dictNumberToSample[key] -=1
        
        return fixRoundingProblem(dictNumberToSample, totalSampleSize)
        

    
    

def getDictsOfCellTypeProportions(num_simulations:int,
                                dictBaselineProportion:dict,
                                totalSampleSize:int,
                                variance_factor:float):
    
    # Initiate a list to hold dictionaries, where each dictionary cotnains the proportions of cells we need to sample to create a simulated bulk sample
    loDictsCellTypeProportions = []
        
        
    for num_simulation in range(1, num_simulations+1):
        dict_cellTypeComposition = getProportionToSample(dictBaselineProportion = dictBaselineProportion,
                                                   variance_factor = variance_factor,
                                                   totalSampleSize=totalSampleSize)

        loDictsCellTypeProportions.append(dict_cellTypeComposition)
        
    
    return loDictsCellTypeProportions

def simulateBulk(cellTypeComposition: pd.Series, df:pd.DataFrame) -> pd.DataFrame:
    """
    Subsamples an dataframe of cell type profiles based on a given cell type composition. 
    
    Returns an pdDataframe object that is basically a simulated bulk dataset
    
    cellTypeComposition = series where indexes are cell type, and value is the number of cells of that cell type to sample
    df =  cell type profile database
    """
    
    print(f"Cell type composition {cellTypeComposition}")
    print(f'df {df.head()}')
    
    loSubsampledCells = [] # Holds cell type profiles that have been multiplied by a scalar value (this is the same as sampling the cell n numer of times and then taking the sum.)
    
    for i, n_cells in enumerate(cellTypeComposition):
        
        # get the cell type we need to sample
        cell_type = cellTypeComposition.index[i]
        print(f'cell type {cell_type}')
        print(f'n_cells {n_cells}')
        #
        
        # For each cell type, sample that number of cells from "df"
        cell_type_profile_multiplied = repeatCellTypes(df = df, cell_type=cell_type, n_cells=n_cells)
        
        print(f'multiplied profile {cell_type_profile_multiplied}')
        
        
        # print diagnostic
        #print(f'{df_cells.shape[0]} cells of {cell_type} were subsampled')
        
        # Append those cells to our list of cells but only if the df_cells is not empty. If its empty then we did not smaple any of those cell types
        if len(cell_type_profile_multiplied) != 0:
            loSubsampledCells.append(cell_type_profile_multiplied)
    
    # Merge the list of dataframes into one large dataframe. This is one simulated sample that has yet to be compacted
    df_simulatedBulk_uncollapsed = pd.concat(loSubsampledCells, axis = 0)
        
    print(df_simulatedBulk_uncollapsed.head())
    
    return df_simulatedBulk_uncollapsed
        
def repeatCellTypes(df:pd.DataFrame, cell_type:str, n_cells:int) -> pd.Series:
    """For the df, which just has cell type profiles, and for one cell type. Return the cell type profile vector multiplied by n_cells

    If no cells are sampled. returns an empty dataframe
    Args:
        df (pd.DataFrame): df with cell type profiles
        cell_type (str): cell type name
        n_cells (int): how many times to multiple the cell type profile (same as sammpling the cell)

    Returns:
        pd.Series: Our cell type profile of our selected cell_type that has been multiplied by "n_cells".
    """

    # First filter the df for that cell type. This will be a series where values are gene expression levels
    cell_type_series = df[df.index  == cell_type]
    
    # CHeck to make sure we succesfully got our cell type profile
    if len(cell_type_series) == 0:
        raise ValueError(f"No Cell types of {cell_type} were found")
    
    # If n_cells is 0, then we will have to jsut return nothing
    if n_cells == 0:
        return pd.Series()

    return cell_type_series * n_cells
    

def checkEnoughCells(df:pd.DataFrame, n_cells: int) -> bool:
    """Retusn true if there is enough cells to sample without replacement
    False otherwise

    Args:
       df (pd.DataFrame): _description_
        n_cells (int): _description_

    Returns:
        bool: _description_
    """
    
    if df.shape[0] >= n_cells:
        return True
    else: 
        return False
    

def getSeed():

    # Generate a random seed for the random module
    random_seed = random.randint(1, 10000)

    print("Random Seed:", random_seed)
    
    return random_seed


    
    
def getBaselineProportion(CTProfile_name:str, proportions_json_path:json) -> dict:
    """Get the list of proper cell type baseline proportions depending on if we are using brain or pbmc

    Args:
        CTProfile_name (str): name of the dataframe of cell type profiles
        proportions_json (json): json that contains baseline cell type profiles for brian or pbmc

    Returns:
        dict: A dict where keys are cell types, and values are their expected baseline proportions
    """
    
    # open json file
    with open(proportions_json_path, 'r') as json_file:
        dict_propostions_json = json.load(json_file)
    
    # Initiate a correct key as False
    correct_key = False
    
    # Get the correct key
    for key in dict_propostions_json.keys():
        if key == f"{CTProfile_name}.csv":
            correct_key = key
    # If no key was found,raise error
    if not correct_key:
        raise ValueError("Key not found")
    
    # Return the proportions of cell types based on our key (which is a tissue type)
    return dict_propostions_json[correct_key]

    
    
    
    

    
    
    
    

if __name__ == "__main__":
    main()