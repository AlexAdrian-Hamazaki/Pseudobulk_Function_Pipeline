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
    
    # Stop the pipeline if I have a cell type in my dictBaselineProportion that I do not have a cell type profile for
    assert checkSameCTs(dictBaselineProportion, df)
    
    
    # Get a list of dictinoaries, where each dictionary represents a simulated bulk sample. 
    # Each dictionary tells us the proportion of each cell type we want to sample for that simulated bulk sample
    loDictsCellTypeProportions = getDictsOfCellTypeProportions(num_simulations,
                                                               dictBaselineProportion,
                                                               totalSampleSize,
                                                               variance_factor)


    # Actually get the number of cells we need to sample for each simulated bulk dataset
    loDictsNumberOfCells = [getNumberToSample(dictCellTypeProportion, totalSampleSize=totalSampleSize) for dictCellTypeProportion in loDictsCellTypeProportions]
    
    # Turn all the dictionaries into pandas series
    loSerNumberOfCells = [pd.Series(dictNumberOfCells) for dictNumberOfCells in loDictsNumberOfCells]
    
    # Writethe compositons for each of the simulated bulk dataset for diagnostic purposes
    for i, composition in enumerate(loSerNumberOfCells):
        composition.to_csv(f'{CTProfile_name}_n_sim_{i}_profiles.csv', index=True)
    
    # Simulate bulk sample and save results
    master_simulate_bulk_wrapper(loSerNumberOfCells=loSerNumberOfCells,
                                 CT_profile_df=df,
                                 CTProfile_name=CTProfile_name,
                                 control = False)
    
    # Shuffle around the CT Profile information (random control)
    control_CT_profile_df = shuffle_CT_Profile_df(df)
    # Simulate bulk samples with the shuffled CT profiles and save results
    master_simulate_bulk_wrapper(loSerNumberOfCells=loSerNumberOfCells,
                                 CT_profile_df=control_CT_profile_df,
                                 CTProfile_name=CTProfile_name,
                                 control = True)
    
def shuffle_CT_Profile_df(CT_profile_df:pd.DataFrame) -> pd.DataFrame:
    """Randomly shuffle each column in the CT profile DF such that the CT profiles will no longer be comprised of actual true CT profile expression.
    
    THe expression of a gene in each CT profile will ranodmly become the expression of one of our CT profiles
    

    Args:
        CT_profile_df (pd.DataFrame): columns are genes, rows are CT profiles, values are expression

    Returns:
        pd.DataFrame: control_CT_profile_df: Expression of each gene has been shuffled between CT profiles
    """
    control_CT_profile_df = CT_profile_df.apply(np.random.permutation)
    return control_CT_profile_df

    
def master_simulate_bulk_wrapper(loSerNumberOfCells:list,
                                 CT_profile_df:pd.DataFrame,
                                 CTProfile_name:str,
                                 control: bool = False):
    if control == True:
        control_profiles = "_cntrl_"
    else:
        control_profiles = "_exp_"
    
    # For simulated bulk dataset (but its currently just numbers), scale the CT profiles by the number of cells of that cell type we want to sample
    loUncollapsedSimulatedBulkSamples = [simulateBulk(cellTypeComposition, df = CT_profile_df) for cellTypeComposition in loSerNumberOfCells]
    
    # Save the initial simulated bulks that are not collapsed
    for i, uncollapsedSimulatedBulkSample in enumerate(loUncollapsedSimulatedBulkSamples):
        uncollapsedSimulatedBulkSample.to_csv(f'{CTProfile_name}{control_profiles}n_sim_{i}_profiles_uncollapsed.csv', index=True)
    
    # For each simulated single cell dataset, collapse into one single simulted bulk SAMPLE
    loSimulatedBulkSamples = [collapseSimulation(loCellTypeProfiles) for loCellTypeProfiles in loUncollapsedSimulatedBulkSamples]
    
    # Concatenate all the simulated bulk SAMPLES into one simulated bulk DATASET
    df_simulatedBulkDataset = pd.concat(loSimulatedBulkSamples, axis = 1)
    
    # Save the bulk simulated dataset
    df_simulatedBulkDataset.to_csv(f'{CTProfile_name}{control_profiles}.csv.gz', compression='gzip')
    
    
    
def checkSameCTs(dictBaselineProportion:dict, df:pd.DataFrame) -> bool:
    """Returns true if all the keys in dictBaselineProportion can be found as indexes in the df of cell type profiles

    Args:
        dictBaselineProportion (dict): keys are cell type profiles, values are lists of [mean,stdev]
        df (pd.DataFrame): dfof cell type profiles. Indexes are cell type names

    Returns:
        bool: True if all keys in dictBaselineProportion are found as indexes in df, False otherwise
    """
    
    # Get keys in dictbaselineProportion
    lo_CTs = list(dictBaselineProportion.keys())
    
    for CT_name in lo_CTs:
        if not CT_name in (df.index):
            print(f"{CT_name} is found in my dictBaselineProportion, but not as a CT Profile index")
            return False
    return True
        
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
        proportion_to_sample = getProportionCellsToSample(cell_proportion_info=dictBaselineProportion[ct], variance_factor=variance_factor)
        
        # Add this number to our dictionary which contains the number of cels to sample for this simulated bulk sample
        dict_ProportionToSample[ct] = proportion_to_sample
        
    # Remove proportions that are Negative, Negatives will be reset to 0
    dict_ProportionToSample = removeNegatives(dict_ProportionToSample)
    
    # When the variance is really high, sometimes NO cells will be sampled from ANY cell types. In this case, just choose to only sample 1 cell to the total smaple size
    if sum(dict_ProportionToSample.values()) == 0:
        random_key = np.random.choice(list(dict_ProportionToSample.keys()))
        dict_ProportionToSample[random_key] = totalSampleSize
    
       
    # Scale the proportions to summarize to 1
    #scaled_dict_ProportionsToSample = scaleDictTo1(dict_ProportionToSample)
    
    # round based on the number of sig figs we have in totalSampleSize
    #scaled_dict_ProportionsToSample = roundToSigFigs(scaled_dict_ProportionsToSample=scaled_dict_ProportionsToSample, totalSampleSize=totalSampleSize)
    
    return  dict_ProportionToSample 

def removeNegatives(dict_ProportionToSample):
    
    """Sets negative values to for a dictionary
    """
    
    for key in list(dict_ProportionToSample.keys()):
        if dict_ProportionToSample[key] <0:
            dict_ProportionToSample[key] = 0
                
    return dict_ProportionToSample

def getProportionCellsToSample(cell_proportion_info:list, variance_factor:float) -> float:
    
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
    stdev_cell_proportion = cell_proportion_info[1]
    
    proportion_to_sample = np.random.normal(baseline_cell_proportion, (float(variance_factor) * float(stdev_cell_proportion)), size = 1)
    return proportion_to_sample[0] # Return first element because we want a float not a np array

def roundToSigFigs(scaled_dict_ProportionsToSample:dict, totalSampleSize:int):
    
    for key in list(scaled_dict_ProportionsToSample.keys()):
        scaled_dict_ProportionsToSample[key] = round(scaled_dict_ProportionsToSample[key], len(str(totalSampleSize)))
        
    return scaled_dict_ProportionsToSample

def scaleDictTo(dict_ProportionToSample:dict, totalSampleSize:int):
    """Scale dictionary values such that the sum of the dict is equal to totalSampleSize

    Args:
        dict_ProportionToSample (dict): Contaisn the proportions of cel ltypes we need to sample
        
    Returns:
        scaled_dict(dict): the scaled dict where the proportions sum to 1
    """
    total_sum = sum(dict_ProportionToSample.values())
    
    # Scale the dict, round as well
    scaled_dict = {key: round(totalSampleSize*(value/total_sum)) for key,value in dict_ProportionToSample.items()}
        
    # If the sum of cells the dictionary is not equal to the totalSampleSize then re-scale
    if int(sum(scaled_dict.values())) != int(totalSampleSize):
        scaled_dict = fixRoundingProblem(scaled_dict,  totalSampleSize)
    
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
        df (pd.DataFrame): a dataframe where rows are cell type profiles that have been scaled, and columns are genes, values are expression levels

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
        
    # Scale the values to be equal to the total sample size
    scaled_dictNumberToSample = scaleDictTo(dict_ProportionToSample = dictCellTypeProportion, totalSampleSize=totalSampleSize)


        
    return scaled_dictNumberToSample

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
    """
    
    I want to create a bunch of simulated bulk samples. As the first step to do this, I will create a large list.
    
    Each element in the list represents a simulated bulk sample. It is a dictionary. 
    The dictionary's keys are the cell type, and the values are how many of that cell type I want to sample for my simulated bulk sample
    
    

    Args:
        num_simulations (int): The number of bulk samples I want to simulate. loDictsCellTypeProportions will be this length.
        dictBaselineProportion (dict): dictionary where keys are cell types. Values are a list where the first element is the mean proportion for that cell type, the second element is its stdev. This will be used to inform how many cells I will sample
        totalSampleSize (int): How many cell type profiles I want to "sample" for each simulated bulk sample
        variance_factor (float): factor indicating how much composition effect I want to simulate.

    Returns:
        loDictsCellTypeProportions:List(dict). A List of len num_simulations, where each element is a dictionary represents how many cell types of each cell type I want to sample for each simulated bulk sample
    """
    
    # Initiate a list to hold dictionaries, where each dictionary cotnains the proportions of cells we need to sample to create a simulated bulk sample
    loDictsCellTypeProportions = []
        
    # For each simulation, get a dictionary indicating how many cells of each cell type I want to sample
    for num_simulation in range(1, num_simulations+1):
        dict_cellTypeComposition = getProportionToSample(dictBaselineProportion = dictBaselineProportion,
                                                   variance_factor = variance_factor,
                                                   totalSampleSize=totalSampleSize)

        loDictsCellTypeProportions.append(dict_cellTypeComposition)
        
    
    return loDictsCellTypeProportions

def simulateBulk(cellTypeComposition: pd.Series, df:pd.DataFrame) -> pd.DataFrame:
    """
    Simulate One Bulk Dataset. 
    
    The bulk dataset is made by scaling cell type profiles by a number provided by cellTypeComposition
    
    This dataset is NOT collapsed, so we can still see the contribution of each cell type profile
    

    cellTypeComposition = series where indexes are cell type, and value is the number of cells of that cell type to sample
    df =  cell type profile database
    
    Returns: df_simulatedBulk_uncollapsed (pd.DataFrame), a dataframe where columns are genes, rows are ct profiles, values are expression. This is an uncollapsed simulated bulk dataset
    """
 
    # Holds cell type profiles that have been multiplied by a scalar value (this is the same as sampling the cell n numer of times and then taking the sum.)
    loMultipliedProfiles = []
    
    # For each cell type in the cellTypeComposition series, multiple that cell type profiles by the value in cellTypeComposition
    for i, n_cells in enumerate(cellTypeComposition):
        #i is index
        #n_cells is the number of cells we want to sample/scale by for this cell type
        
        # get the name of the cell type that we need to sample
        cell_type = cellTypeComposition.index[i]
        
        # For each cell type, scale the cell type profile by n_cells
        cell_type_profile_multiplied = repeatCellTypes(df = df, cell_type=cell_type, n_cells=n_cells)
        
        # Append those cells to our list of cells but only if the df_cells is not empty. If its empty then we did not smaple any of those cell types
        if len(cell_type_profile_multiplied) != 0:
            
            loMultipliedProfiles.append(cell_type_profile_multiplied)
    
    # Merge the list of dataframes into one large dataframe. This is one simulated bulk sample that has yet to be compacted
    df_simulatedBulk_uncollapsed = pd.concat(loMultipliedProfiles, axis = 0)
    
    # round all the values
    df_simulatedBulk_uncollapsed = df_simulatedBulk_uncollapsed.round(4)
        
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
    cell_type_profile = df[df.index  == cell_type]
    
    # CHeck to make sure we succesfully got our cell type profile
    if len(cell_type_profile) == 0:
        raise ValueError(f"No Cell types of {cell_type} were found")
    
    # If n_cells is 0, then we will have to jsut return nothing
    if n_cells == 0:
        return pd.Series()

    # If all is good, jsut scale the cell type profile by n_cells
    return cell_type_profile * n_cells
    

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