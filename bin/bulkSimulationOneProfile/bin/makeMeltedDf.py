#!/usr/bin/env python3

import pandas as pd
import sys



def main():
    lo_EGADs = sys.argv[1:]
    #lo_EGADs = os.listdir("/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/bulkSimulationOneProfile/data/dev/EGAD/brain_sc_with_metadata_cpm_pc_cell_type_profiles_cntrl_.csv") + os.listdir("/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/bulkSimulationOneProfile/data/dev/EGAD/brain_sc_with_metadata_cpm_pc_cell_type_profiles_exp_.csv")

    # First, split the EGADs into organism parts such that keys are OPs and values are paths
    split_EGAD_dict = split_list_into_OPs(lo_EGADs)


    # Process one organism part at a time. Each element will be a melted dataframe for an organism part
    for OP_key, lo_paths in split_EGAD_dict.items():
        melted_df = process_organism_part(OP_key, lo_paths)
        melted_df.to_csv(f"{OP_key}_melted_df.csv", index=False)




def process_organism_part(OP_key:str, lo_paths:list):
    # First split the lists into experimental and control
    lo_exp_paths = list(filter(lambda x: 'exp' in x, lo_paths))
    lo_contrl_paths = list(filter(lambda x: 'cntrl' in x, lo_paths))
    
    # process each of the paths
    melted_exp_df = get_melted_df(lo_exp_paths, 'exp')
    melted_contrl_df = get_melted_df(lo_contrl_paths, 'contrl')
    
    # concate the melted dfs row wise
    return pd.concat([melted_exp_df, melted_contrl_df], axis = 0).sort_values(by = 'variable', ascending=True)


    
def get_melted_df(lo_paths, type:str):
	# FOr each of these paths, load pandas dataframe
	lo_dfs = [pd.read_csv(path) for path in lo_paths]

	# Parse the path for the variance level
	lo_variance_levels = [get_variance_int(path) for path in lo_paths]

	# MAke the AUC column equal to the variance level. Also gets a list of pd series
	lo_auc_series = process_auc_colnames(lo_dfs, lo_variance_levels)
 
	# Concat all of the dfs together
	concat_df = pd.concat(lo_auc_series, axis = 1)

	# Melt the df
	melted_df = concat_df.reset_index().melt(id_vars = ['index'])

	# Add a colum indicating if the experiment is contrl or exp
	melted_df['type'] = type
	return melted_df

def process_auc_colnames(lo_dfs:list, lo_variance_levels:list):
	lo_ser_pro = [] # list hold AUC vectors wwhere names are the variance level
	for i, df in enumerate(lo_dfs):
		var_int = lo_variance_levels[i]
		# Renaming the 'auc' column to 'var'
		
		df_pro = df.rename(columns={'auc': str(var_int)})
		auc_ser = df_pro.loc[:,str(var_int)]

		lo_ser_pro.append(auc_ser)

	# if len(lo_ser_pro) < 2:
	# 	print(lo_dfs)
	# 	assert False
	
	return lo_ser_pro

def get_variance_int(path:str) -> str:
    return str(path.split('_')[-2])

def split_list_into_OPs(lo_EGADs:list) -> dict:
    # Get a list of the unique OP NAmes
    unique_OPs = get_unique_OPs(lo_EGADs)
    
    # Make a dict where the keys are OPs and values is a list of all the paths affiliated with that OP
    EGAD_dict = make_EGAD_dict(unique_OPs, lo_EGADs)
    
    return EGAD_dict

    
def get_unique_OPs(lo_EGADs):
    # Get the organism nane
    lo_names = []
    
    for path in lo_EGADs:
        name = path.split("_")[1]
        lo_names.append(name)
        
    lo_names = list(set(lo_names))
    
    return lo_names

def make_EGAD_dict(unique_OPs:list, lo_EGADs:list):
    EGAD_dict = {}
    
    for OP in unique_OPs:
        # Using the filter function along with a lambda function to check for the substring
        filtered_list = list(filter(lambda x: OP in x, lo_EGADs))
        
        # Add to dict
        EGAD_dict[OP] = filtered_list
        
    return EGAD_dict


if __name__ == "__main__":
	main()