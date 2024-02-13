#!/usr/bin/env python3


import sys
import scanpy as sc
import pandas as pd
import random
import anndata as ad
import numpy as np
import os
from scipy.sparse import csr_matrix


def main():
    #seed_value = 42  # Replace with your desired seed value
    #np.random.seed(seed_value)
    
    h5ad_path = sys.argv[1]
    h5ad_name = sys.argv[2]
    
    num_simulations = int(sys.argv[3])
    totalSampleSize = int(sys.argv[4])
    variance_factor = float(sys.argv[5])
    
    print(h5ad_name)
    print(num_simulations)
    print(totalSampleSize)
    print(variance_factor)

    # load adata
    adata = sc.read_h5ad(h5ad_path)
    # Identify the unique CTs that are in this sample
    loCTs = adata.obs.loc[:,"Cell type"].unique()
    
    print(f"List of all the unique Cell Types")
    
    # init list
    loCellTypeCompositions = [] # List to hold the different cell type compositions for each simulation. These should be pd.Series
    
    # get NUMBERS for the amount of cells in each simulation that we want to subsample
    for n_simulation in range(1, num_simulations+1):
        cellTypeComposition = getProportionToSample(dictBaselineProportion=dictBaselineProportions,
                                                    variance_factor = variance_factor,
                                                    totalSampleSize = totalSampleSize)

        loCellTypeCompositions.append(cellTypeComposition)
    

    # Get a list that contains subsampled H5ADs. Each element is a simulated bulk dataset
    loSubsampledH5ADs = [subsampleH5AD(cellTypeComposition = cellTypeComposition,
                                     adata = adata) for cellTypeComposition in loCellTypeCompositions]
    
    
    
    # Make the right directory
    os.mkdir(h5ad_name)
    
    # Write each simulated bulk dataset.
    # For writing purposes, compress adata.X as well
    for i, subsampledH5AD in enumerate(loSubsampledH5ADs):
        subsampledH5AD.X = csr_matrix(subsampledH5AD.X)
        subsampledH5AD.write_h5ad(f'n_sim_{i}.h5ad', compression='gzip', as_dense=("X"))
    
    # Also write the compositons for each of the simulated bulk dataset
    for i, composition in enumerate(loCellTypeCompositions):
        composition.to_csv(f'n_sim_{i}.csv', index=True)
        
    
        


def subsampleH5AD(cellTypeComposition: pd.Series, adata:sc.AnnData) -> sc.AnnData:
    """
    Subsamples an h5ad object based on a given cell type composition. 
    
    Returns an anndata object that is basically a simulated bulk sample
    """
    
    loH5ADs = [] # Holds the subsampled h5ads
    
    for i, n_cells in enumerate(cellTypeComposition):
        
        # get the cell type we need to sample
        cell_type = cellTypeComposition.index[i]
        

        # For each cell type, get n_cells number of those cells from the adata object
        cells = subsampleCellType(adata = adata, cell_type=cell_type, n_cells=n_cells)
        
        # print diagnostic
        print(f'{len(cells)} cells of {cell_type} were subsampled')
        
        #Append those cells to 
        loH5ADs.append(cells)
    
    # Merge the list of h5ads into 1 h5ad object 
    adata = ad.concat(loH5ADs, axis = 0)
        
    return adata
        
    
def subsampleCellType(adata : sc.AnnData, cell_type:str, n_cells:int) -> sc.AnnData:
    """
    
    For each cell type, get n_cells number of those cells from the adata object


    Args:
        adata (sc.AnnData): adata object
        cell_type (str): a string cell type
        n_cells (int): the number of cells to subsample

    Returns:
        sc.AnnData: ann anndata object that has only our sampled cells
    """

    # Slice data to only get one cell type
    adata = adata[adata.obs.loc[:,"Cell type"] == cell_type]
    print(adata)
    print(cell_type)
    # check if there is enough cells to sample with replacement
    is_enough = checkEnoughCells(adata = adata, n_cells=n_cells)
    

    
    # if there is more cells in the adata object then n_cells, subsample to n_cells
    if is_enough:
        # Get a random seed 
        rand_seed = getSeed()
        sc.pp.subsample(adata, n_obs = n_cells, random_state = rand_seed)
    # if there is not enough cells to subsample, do a resampling method
    else: 
        adata = resamplingMethod(adata, n_obs = n_cells)
    
    
    return adata

def resamplingMethod(adata:ad.AnnData, n_obs:int)-> ad.AnnData:
    loAdataCells = []
    print(adata)
    
    # For each n_obs, we want to sample one cell. That cell is that n_obs we've sampled
    for _ in range(n_obs):
        # Randomly get a seed and set it
        rand_seed = getSeed()
        np.random.seed(rand_seed)
        
        # Indicate an indice to sample from the adata object
        indice = np.random.choice(range(adata.shape[0]), size=1, replace=True)
        print(f'Random Indice {indice}')
        
        # sample that cell from the adata object
        sampled_cell = adata[indice,:]
        print(f'Random Cell {sampled_cell.obs_names}')
        print(f'Random Cell {sampled_cell.X}')

        
        # Add the sampled cell to the adata object
        loAdataCells.append(sampled_cell)
        print('')
    print("~~~~~~~~~~~~~")
    print(adata)
    
    print(loAdataCells)
    # Merge the cells together into an adata object
    adata_subsampled = ad.concat(loAdataCells, axis = 0)
    adata_subsampled.obs_names_make_unique()

    return adata_subsampled

def getSeed():

    # Generate a random seed for the random module
    random_seed = random.randint(1, 10000)

    print("Random Seed:", random_seed)
    
    return random_seed
    
def checkEnoughCells(adata:sc.AnnData, n_cells: int) -> bool:
    """Retusn true if there is enough cells to sample without replacement
    False otherwise

    Args:
        adata (sc.AnnData): _description_
        n_cells (int): _description_

    Returns:
        bool: _description_
    """
    
    if adata.shape[0] >= n_cells:
        return True
    else: 
        return False
    
# def getNumberToSubsample(loCTs:list, totalSampleSize:int, variance_factor:float) -> pd.Series:
#     """For one tissue, and given the number of unique cell types, and given a total sample size for our simulated bulk,
#     identify how many cells of each cell type we want to sample"

#     Args:
#         loCTs (list): List of the names of the unique cel ltypes
#         totalSampleSize (int): The total number of cells we want in each simulated bulk dataset
#         variance_factor (float): The amount of the compositional effect we want to simulate. There is actually a direct interpretation of this number though. For a simulated bulk of 4 cell types and for 100 samples, \
#             a variance of 0.01 indicates that the composition of each cell type can differ by 1% of the total number of cells. \
#                 For example, with zero composition effect, each cell type would have 25,25,25,25 cells. With variance=0.01, \
#                     the number of cells in each cell type can vary by 100*0.01 = 1 cell. So each cell type could have
#                     24-26 cells in it

#     Returns:
#         pd.Series: Named series, has the number of cells that each cell type needs to have.
#     """
    
    
#     loNumberOfSamples = []  # List to store the number of samples per group
#     remaining_samples = totalSampleSize  # Initialize remaining samples to the total number

#     # Loop through each group except the last one
#     for i, group in enumerate(loCTs[:-1]):
#         # Generate a random percentage within a certain range influenced by variance_factor
#         min_samples = max(0, remaining_samples / (len(loCTs) - i) - variance_factor * totalSampleSize)
#         max_samples = min(remaining_samples / (len(loCTs) - i) + variance_factor * totalSampleSize, remaining_samples)
        
#         # Get the number of cells we will sample
#         samples = int(random.uniform(min_samples, max_samples))
        
#         # print(f'min {min_samples}')
#         # print(f'max {max_samples}')
#         # print(f'per {samples}')
        
#         # turn it into an int
#         num_samples = int(samples)
        
#         # Add the number of samples we will sample to the list
#         loNumberOfSamples.append(num_samples) 
    
#         # Remove these samples
#         remaining_samples -= num_samples  # Subtract the allocated samples from the remaining

#     # Allocate the remaining samples to the last group so that the last group is there
#     loNumberOfSamples.append(remaining_samples)
    
#     random.shuffle(loNumberOfSamples)  # Shuffle the order of allocated samples among groups
    
    
#     return pd.Series(loNumberOfSamples, index = loCTs, name = 'numbers_to_sample')


def getProportionToSample(dictBaselineProportion:dict,  variance_factor:float, totalSampleSize:int) ->  dict:
    """
    Identify how many cells of each cell type we want to sample for 1 simulated bulk sample
    
    The default number of cells you want to sample will be based on the dictbaselineProportion, the default percent composition that I'm setting for each cell type
    dictbaseilne proprotion
    
    Also performs some scaling to ensure the number of cells we sample will be ==  totalSampleSize

    Args:
        dictBaselineProportion:dict: dictionary where keys are cell types. Values are a list where the first element is the mean proportion for that cell type, the second element is its stdev. This will be used to inform how many cells I will sample
        totalSampleSize (int): The total number of cells we want in each simulated bulk dataset
        variance_factor (float): The amount of the compositional effect we want to simulate.
    Returns:
        scaled_dict_ProportionsToSample (dict): A Simulated bulk sample in the form of a dictionary, Keys are cell types and values are how many cells of that cell type we want to sample
    """

    dict_ProportionToSample = {}  # dict to store the number of cells for each cell type. Keys are the cell type, and values are the number of cells of that cell type to sample

    # Get all the keys (cell types)
    loCts = list(dictBaselineProportion.keys())
    
    # Identify how many cells of each cell type we want to sample, (keys are cell types)
    for ct in loCts:
        # for each cell type, get a value to sample
        proportion_to_sample = getProportionCellsToSample(cell_proportion_info=dictBaselineProportion[ct], variance_factor=variance_factor, totalSampleSize=totalSampleSize)
        
        # Add this number to our dictionary which contains the number of cels to sample for this simulated bulk sample
        dict_ProportionToSample[ct] = proportion_to_sample
        
    # Remove proportions that are Negative, Negatives will be reset to 0
    dict_ProportionToSample = removeNegatives(dict_ProportionToSample)
    
    # When the variance is really high, sometimes NO cells will be sampled from ANY cell types. In this case, just choose to only sample 1 cell to the total smaple size
    if sum(dict_ProportionToSample.values()) == 0:
        random_key = np.random.choice(list(dict_ProportionToSample.keys()))
        dict_ProportionToSample[random_key] = totalSampleSize
    
    return  dict_ProportionToSample 


def getProportionCellsToSample(cell_proportion_info:list, variance_factor:float, totalSampleSize:int) -> float:
    
    """
    Get an int value that represents how many cells we want to sample.

    This is done by sampling an int from a normal distribution of a set mean and standard deviation
    Args:
        cell_proportion_info (list[mean,stdev]): A List where the 0th element is the mean of the sampling distribution, and the 1rst element is the stdev
        variance_factor (float): The amount of the compositional effect we want to simulate.

    Returns:
        float: The proportion of cells of this cell type that we want to sample
    """
    # Get the "baseline" or "mean" value of this cell type's proportion
    baseline_cell_proportion = cell_proportion_info[0]
    
    # Get the stdev of this cell's sampling
    # stdev_cell_proportion = cell_proportion_info[1]
    stdev_cell_proportion = 1

    proportion_to_sample = np.random.normal(baseline_cell_proportion, float(variance_factor)*2*float(stdev_cell_proportion), size = 1)
    #proportion_to_sample = round(proportion_to_sample[0], len(str(totalSampleSize))) # round to the sig figs or how many cells we are sampling
    return proportion_to_sample[0] # Return first element because we want a float not a np array


def removeNegatives(dict_ProportionToSample):
    
    """Sets negative values to for a dictionary
    """
    
    for key in list(dict_ProportionToSample.keys()):
        if dict_ProportionToSample[key] <0:
            dict_ProportionToSample[key] = 0
                
    return dict_ProportionToSample

if __name__ == "__main__":
    main()