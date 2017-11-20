# -*- coding: utf-8 -*-
"""
Analysis of the output of the simulation
"""
os.chdir('/home/cmee19/Documents/Project/Sim_TH/Metacom_sim/')
import os
import pandas as pd
import scipy as sc
import pickle
import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.pyplot as plt
from simulation_metacom_class import *

#%%
#abundance_vect = abund_data.iloc[:,0].values #debug
#abundance_vect_nozeros = abundance_vect[abundance_vect != 0] #debug

def oct_calc(abundance_vect):
    """Calculates the octaves of an abundance array.
    abundance_vect: pandas.series of abundances
    """    
    #converts abundances into which bin they belong to. vector operation
    octaves_sps_belong = sc.floor(sc.log2(abundance_vect)).astype(int)
    #how many octaves for this community? how long the octaves vect needs to be add
    #for the 0. the max + 1 to count the 0 as well
    nbins = sc.amax(octaves_sps_belong) + 1
    oct_vect = [0] * nbins
    for i in octaves_sps_belong:
        oct_vect[i] = oct_vect[i] + 1
    
    return oct_vect

def samples_oct_list(abund_data):
    """Calculates the octaves of a set of communities in a pandas dataframe
    abund_data: Pandas dataframe where the communies are in columns and the species
    in rows
    """    
    oct_list = []
    for i in abund_data.columns:
        abundance_vect_nozeros = abund_data[abund_data.loc[:,i] != 0].loc[:,i]        
        oct_list.append(oct_calc(abundance_vect_nozeros))
    
    return oct_list


#oct_data_list = samples_oct_list(abund_data)
#oct_sim_list = samples_oct_list(loaded_sim.sp_abund_df.iloc[:,1:11])
##oct_abund_list1 = oct_data_list
#oct_abund_list2 = oct_sim_list
#i = 0

def diff_index_calc(oct_abund_list1, oct_abund_list2):
    """Calculates and index of disimilarity between two lists of octaves vectors.
    Oct_abund_list: List of list of octaves, one list of octaves for each tree hole
    It's meant to compare the octaves of real data and simulated data.
    The global difference between the octaves of all tree holes is summarize summing
    the difference for each tree hole.
    Output: tuple(relative index / abusolute index) 
    relative index: is how big is the difference between each of the bins relative 
    to the mean between the two     
    absolute index: total number of individulas that are in different bins    
    """
    rel_index_list = []
    abs_index_list = []
    smty_index_list = []
    for i in range(10):
        abund_data_array = sc.asarray(oct_abund_list1[i], dtype='double')
        abund_sim_array = sc.asarray(oct_abund_list2[i], dtype = 'double')
        
        # make the length of the arrays similar to each other
        if len(abund_data_array) < len(abund_sim_array):
            small_len = abund_data_array
            long_len = abund_sim_array
        else:
            small_len = abund_sim_array
            long_len = abund_data_array
        diff = len(long_len) - len(small_len)        
        small_len = sc.append(small_len, [0]*diff)
           
        relative_index_vect = abs(long_len - small_len)/long_len   
        rel_index_list.append(sum(relative_index_vect)/len(relative_index_vect))
        
        absolute_index_vect = abs(long_len - small_len)
        abs_index_list.append(sum(absolute_index_vect)/len(absolute_index_vect))
        
        similarity_index_vect = []
        for i in range(len(long_len)):
            similarity_index_vect.append(sc.minimum(long_len[i], small_len[i])/sc.amax([long_len[i], small_len[i]]))
            
        smty_index_list.append(sum(similarity_index_vect)/len(similarity_index_vect))         
        
    rel_index_final = sum(rel_index_list)/10
    abs_index_final = sum(abs_index_list)/10
    smty_index_final = sum(smty_index_list)/10
    
    return (rel_index_final, abs_index_final, smty_index_final)

#diff_index_calc(oct_data_list, oct_sim_list)

#loaded_sim = pickle.load(open('TH_metacom_sim1.pkl', 'rb'))
#b = pickle.load(open('Third_sim/TH_metacom_sim2.pkl', 'rb'))

def parms_comparator(path_to_simfiles, abund_data, n_sim):
    """Builds a dataframe with the parameters and the differece indices between the 
    octaves.
    path_to_simfiles: path to the folder where the files of the simulation are
    abund_data: dataframe with the abundance data of only one time point, and there 
    should be 10 samples/communities.
    n_sim: number of simulations
    """
    
    
    diff_df = pd.DataFrame(columns = ['m', 'cln', 'rel_diff_index', 'abs_diff_index', 'smty_index_index'])
    for i in range(1,26):
        try:
            loaded_sim = pickle.load(open(os.path.join(path_to_simfiles,'TH_metacom_sim%s.pkl'%i), 'rb'))
        except:
            continue
            
        oct_sim_list = samples_oct_list(loaded_sim.sp_abund_df.iloc[:,1:11])    

        #I need the size of each tree hole to sample the same amount of ind from the data
        th_size = loaded_sim.sp_abund_df.sum(0)        
        # in order to sample abund_data I need to build like a metacommunity, each
        # sps has a weight
        abund_data_total_th = abund_data.sum(axis = 0)
        metacom = abund_data.div(abund_data_total_th, axis = 1)
        # sampling the each of the holes
        i = 0
        oct_sampled_superlist = []
        for i in range(10):        
            sampled_data = pd.DataFrame()        
            for i in range(10):
                sampled_data = pd.concat([sampled_data, pd.Series(sc.random.choice(metacom.index, p = metacom.iloc[:,i], size = th_size[i+1])).value_counts()], axis = 1)
            sampled_data.columns = abund_data.columns    
            sampled_data.fillna(0, inplace = True)
            
            oct_sampled_list = samples_oct_list(sampled_data)
            [len(i) for i in oct_sampled_list]
            oct_sampled_superlist.append(pd.Series(samples_oct_list(sampled_data)))
        oct_sampled_superlist[2][3][0]
            
        diff = diff_index_calc(oct_sim_list, oct_data_list)        
        #row = pd.DataFrame([loaded_sim.m, loaded_sim.cln, diff])    
        diff_df = diff_df.append({'m': loaded_sim.m, 'cln':loaded_sim.cln, 'rel_diff_index':diff[0], 'abs_diff_index':diff[1], 'smty_index_final':diff[2]},ignore_index = True)

    return diff_df



#%%

abund = pd.read_csv("../../Data/Processed_data/Samples1K_OTUs100/TH_abund.csv", header = 0, index_col = 0)
abiotic = pd.read_csv('../../Data/Processed_data/Samples1K_OTUs100/TreeHoles_Abiotic.csv', header = 0)
# I cant index the dates because the name of the column has a dot, and python
# recognizes as a method
#abund_data = abund[abiotic[abiotic.R\.date == '2015-09-29'].loc[:,'Location']]

# Change the name of the R.date column
name_list = abiotic.columns.tolist()
name_list[3] = "R_date"
abiotic.columns = name_list
# Select the samples from only one time point that has samples for the 10 tree holes
samples = abiotic[abiotic.R_date == '2015-09-29'].loc[:,['Samples_shortname', 'Location']]
samples.sort_values('Location', axis = 0, inplace = True)
samples = samples.loc[:,'Samples_shortname']
abund_data = abund[samples]
abund_data

loaded_sim = pickle.load(open('TH_metacom_sim1.pkl', 'rb'))

#%% Running the function

# parameter comparator function only works for a group of functions

#%% plotting

#oct_data_list = samples_oct_list(abund_data)
#oct_vect = oct_data_list[2]
#   
#upper = [1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192]
#plt.bar(sc.arange(len(oct_vect)), oct_vect)
#plt.xticks(sc.arange(len(oct_vect)), upper[0:len(oct_vect)+1])
#plt.show()


