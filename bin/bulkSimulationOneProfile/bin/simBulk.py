#!/usr/bin/env python3

import sys
import anndata as ad
import pandas as pd
import json
import random
import scanpy as sc


def main():
    CTProfile_path = sys.argv[1]
    CTProfile_name = sys.argv[2]
    proportions_json_path = sys.argv[3]
    variance_factor = sys.argv[4]
    
    
    num_simulations = 100
    totalSampleSize = 1000
    # CTProfile_path = "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/bulkSimulationOneProfile/work/2a/629409c1815c9a57d23fdc33b86185/brain_sc_with_metadata_cpm_pc_cell_type_profiles.csv"
    # CTProfile_name = "brain_sc_with_metadata_cpm_pc_cell_type_profiles"
    # variance_factor = 0.1
    # proportions_json_path = "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/bulkSimulationOneProfile/work/2a/629409c1815c9a57d23fdc33b86185/cell_type_proportions.json"
    
    # open cell type profile database
    df = pd.read_csv(CTProfile_path,index_col=0)

    
    # Retreive the baseline proportions from the proportions_json
    dictBaselineProportion = getBaselineProportion(CTProfile_name=CTProfile_name, proportions_json_path=proportions_json_path)
    
    
    # Initiate a list to hold a list to hold lists, where each list contains how many cells of each cell type to sample for a given simulation
    loloCellTypeCompositions = []
    
    # Where each element in the list is the cell type proportions I'm going for
    
    for num_simulation in range(1, num_simulations+1):
        cellTypeComposition = getNumberToSubsample(dictBaselineProportion = dictBaselineProportion,
                                                   totalSampleSize= totalSampleSize,
                                                   variance_factor = variance_factor)

        loloCellTypeCompositions.append(cellTypeComposition)
        
    # For each of the cell type compositions, get a simulated single cell dataset.
    loloCellTypeProfiles = [simulateBulk(cellTypeComposition, df = df) for cellTypeComposition in loloCellTypeCompositions]
    print(len(loloCellTypeProfiles))
    
        #~~~ UP TO HERE IS CHECKED AND MAKES SENSE

    # For each simulated single cell dataset, collapse into one single simulted bulk sample
    loSimulatedBulkSamples = [collapseSimulation(loCellTypeProfiles) for loCellTypeProfiles in loloCellTypeProfiles]
    print(len(loSimulatedBulkSamples))
    
    
    # Make a simulated bulk dataset where each row is a simulated bulk sample
    df_simulatedBulkDataset = pd.concat(loSimulatedBulkSamples, axis = 1)
    print(df_simulatedBulkDataset)
    
    # Save the bulk simulated dataset
    df_simulatedBulkDataset.to_csv(f'{CTProfile_name}.csv.gz', compression='gzip')
    
     # EVERYTHING IS CHECKED, THE ONLY DIFFERENCE SEEMS TO BE THE VARIANCE LEVEL
    # Also write the compositons for each of the simulated bulk dataset
    for i, composition in enumerate(loloCellTypeCompositions):
        composition.to_csv(f'{CTProfile_name}_n_sim_{i}_profiles.csv', index=True)
    
def collapseSimulation(df: pd.DataFrame) -> pd.Series:
    """Each dataframe is a bunch of cell type profiles.
    
    This function just sums up all othe expression and returns a pandas series which represents a simulated bulk sample

    Args:
        df (pd.DataFrame): _description_

    Returns:
        pd.Series: A bulk simulated sample
    """
    
    return df.sum(axis = 0)

def simulateBulk(cellTypeComposition: pd.Series, df:pd.DataFrame) -> pd.DataFrame:
    """
    Subsamples an dataframe based on a given cell type composition. 
    
    Returns an pdDataframe object that is basically a simulated bulk dataset
    """
    
    loCells = [] # Holds the subsampled cells from a variety of cell types
    
    for i, n_cells in enumerate(cellTypeComposition):
        
        # get the cell type we need to sample
        cell_type = cellTypeComposition.index[i]
        
        # For each cell type, sample that number of cells from "df"
        df_cells = repeatCellTypes(df = df, cell_type=cell_type, n_cells=n_cells)
        
        # print diagnostic
        print(f'{df_cells.shape[0]} cells of {cell_type} were subsampled')
        
        # Append those cells to our list of cells but only if the df_cells is not empty. If its empty then we did not smaple any of those cell types
        if df_cells.empty == False:
            loCells.append(df_cells)
    
    # Merge the list of dataframes into one large dataframe. This is one simulated sample that has yet to be compacted
    df_simulatedBulk = pd.concat(loCells, axis = 0)
    
    df_simulatedBulk.round(4)
        
    return df_simulatedBulk
        
def repeatCellTypes(df:pd.DataFrame, cell_type:str, n_cells:int) -> pd.DataFrame:
    """For the df, which just has cell type profiles. return a dataframe with n_cells number of repeats for that cell type profile

    If no cells are sampled. returns an empty dataframe
    Args:
        df (pd.DataFrame): _description_
        cell_type (str): _description_
        n_cells (int): _description_

    Returns:
        pd.DataFrame: _description_
    """
    
    # First filter the df for that cell type. This will be a series
    cell_type_series = df[df.index  == cell_type]
    
    # CHeck to make sure we succesfully got our cell type profile
    if len(cell_type_series) == 0:
        raise ValueError(f"No Cell types of {cell_type} were found")
    
    # If n_cells is 0, then we will have to jsut return nothing
    if n_cells == 0:
        return pd.DataFrame()
    
    
    # initiate a list to hold the cell type series X amount of times
    lo_CellTypeSeries = []

    # print(cell_type_series)
    # print(n_cells)
    # print(df.index)
    # Add the row n_cells-1 number of times 
    for _ in range(n_cells):
        lo_CellTypeSeries.append(cell_type_series)
        
    #  Create a pandas dataframe of those series
    df_CellType = pd.concat(lo_CellTypeSeries, axis = 0)
    
    return df_CellType
          

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

def getNumberToSubsample(dictBaselineProportion:dict, totalSampleSize:int, variance_factor:float) -> pd.Series:
    """For one tissue, and given the number of unique cell types, and given a total sample size for our simulated bulk,
    identify how many cells of each cell type we want to sample"
    
    The default number of cells you want to sample will be based on the dictbaselineProportion, the default percent composition that I'm setting for each cell type

    You first have to randomly shuffle the proportions dict though, otherwise at high variances the final cell type never gets sampled


    Args:
        dictBaselineProportion:dict: A dict where keys are cell types, and values are the default percent composition of that cell type.
        totalSampleSize (int): The total number of cells we want in each simulated bulk dataset
        variance_factor (float): The amount of the compositional effect we want to simulate. There is actually a direct interpretation of this number though. For a simulated bulk of 4 cell types and for 100 samples, \
            a variance of 0.01 indicates that the composition of each cell type can differ by 1% of the total number of cells. \
                For example, with zero composition effect, each cell type would have 25,25,25,25 cells. With variance=0.01, \
                    the number of cells in each cell type can vary by 100*0.01 = 1 cell. So each cell type could have
                    24-26 cells in it

    Returns:
        pd.Series: Named series, has the number of cells that each cell type needs to have.
    """
    
    
    dictNumberOfSamples = {}  # dict to store the number of cells for each cell type. Keys are the cell type, and values are the number of cells of that cell type to sample
    remaining_samples = totalSampleSize  # Initialize remaining samples to the total number. This is a counter 

    # GEt the number of unique cell types
    n_cell_types = len(dictBaselineProportion.keys())
    
    # Iniaite a counter to see how many cell types we've subsampled
    n = 0
    
    # Randomly shuffle all the dictionary keys
    loKeys = getShuffledKeys(dictBaselineProportion)
    
    # Loop through each cell type
    for key in loKeys:
        
        # uptick counter
        n+=1
        
        # Get the baseline number of cells.
        baselineCells = dictBaselineProportion[key] * totalSampleSize
        
        # Generate the minumim and maximum number of cells that we could get for this cell type
        min_samples = max(0, baselineCells - float(variance_factor) * float(totalSampleSize))
        max_samples = float(baselineCells) + float(variance_factor) * float(totalSampleSize) # If the number of remaining samples is less than variance. Then just use the remaining samples
        
        # We are on our last Cell Type. Set the number of cells to sample as the remaining cells. Unless that number is below 0. Then make the number of cells to sample 0
        if n == n_cell_types:
            num_samples = remaining_samples
            if num_samples < 0:
                num_samples = 0
        else:
            # Get the number of cells we will sample
            num_samples = int(random.uniform(min_samples, max_samples))
    
        
        # turn it into an int
        num_samples = int(num_samples)
        
        if num_samples < 0:
            raise ValueError("Number of cells to sample shuld not be negative")
        
        # Add the number of samples we will sample to the list
        dictNumberOfSamples[key] = num_samples
    
        # Remove these samples
        remaining_samples -= num_samples  # Subtract the allocated samples from the remainin
    
    return pd.Series(dictNumberOfSamples, name = 'numbers_to_sample')
    
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