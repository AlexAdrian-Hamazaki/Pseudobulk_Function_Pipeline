#!/usr/bin/env python3


import pandas as pd
import sys
import numpy as np



def main():
	ct_profile_path = sys.argv[1]
	step = sys.argv[2]
 
	# load CT Profiles
	ct_pros = pd.read_csv(ct_profile_path, index_col=0	)
	

	# Sample {step} number of ct profiles and save those DFs iteratively
	lo_subsampled_dfs = sample_ct_profiles(ct_pros=ct_pros, step = step)
 

	# save each df
	[df.to_csv(f"{df.shape[0]}_ct_pros.csv") for df in lo_subsampled_dfs]
  
	# When I add this info I get more performance on EGAD. How can I model composition effects here with only CT profiles
	# I can add
 
	

def sample_ct_profiles(ct_pros:pd.DataFrame, step:int):
	"""Sample cell type profiles. start with step, and then add step
 
 Returns a list of dataframes of subsampled ct profiles

	Args:
		ct_pros (pd.DataFrame): _description_
		step (int): _description_
	"""
	step=int(step)
 
	# init list to hold dfs we sample
	lo_subsampled_dfs = []
 
	# Init value to keep us in while look
	sample_again = True
	
	# init how many ct profiles we want to sample
	ct_pro_counter = step
 
	while sample_again:
		if ct_pro_counter >= int(ct_pros.shape[0]): # If we are being asked to subsample more CT profiles than exist in the ct profile dataframe...
			lo_subsampled_dfs.append(ct_pros)
			sample_again = False
     
		subsampled_df = ct_pros.sample(n=ct_pro_counter, axis = 0)
		lo_subsampled_dfs.append(subsampled_df)
		ct_pro_counter += step

	return lo_subsampled_dfs

if __name__ == "__main__":
	main()