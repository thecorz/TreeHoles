
# coding: utf-8

# ## Analysis of the output of the simulation

# The notebook should be run with the simulation_metacom_class in the same folder

# In[5]:

import os
import pandas as pd
import scipy as sc
import pickle
import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.pyplot as plt
from simulation_metacom_class import *


# Function to calculate the octaves

# In[6]:

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


# Function to calculate the octaves of an array

# In[7]:

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


# Function to calculate the indices that measure the difference between the octaves of the real data and the simulated data

# In[49]:

def diff_index_calc(oct_abund_list1, oct_abund_list2):
    """Calculates and index of disimilarity between two lists of octaves vectors
    It's meant to compare the octaves of real data and simulated data.
    oct_abund_list1/2 list of octaves, where octaves are lists as well
    Output: tuple(relative index / abusolute index) 
    relative index: is how big is the difference between each of the bins relative 
    to the mean between the two     
    absolute index: total number of individulas that are in different bins    
    """
    rel_index_list = []
    abs_index_list = []
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
           
        relative_index_vect = abs(long_len - small_len)/((long_len + small_len)/2)   
        rel_index_list.append(sum(sc.asarray(relative_index_vect)))
        
        absolute_index_vect = abs(long_len - small_len)
        abs_index_list(sum(sc.asarray(absolute_index_vect)))
        
    rel_index_final = sum(sc.asarray(rel_index_list))
    abs_index_final = sum(sc.asarray(abs_index_list))
    
    return (rel_index_final, abs_index_final)
diff_index_calc(oct_data_list, oct_data_list)


# Function to compare the octaves of the data and the octaves of the simulation. It requires that the real data to compare with has to have 10 samples/communities like the number of communities in the simulation, to be able to compare 1:1

# In[37]:

def parms_comparator(path_to_simfiles, abund_data, n_sim):
    """Builds a dataframe with the parameters and the differece indices between the 
    octaves.
    path_to_simfiles: path to the folder where the files of the simulation are
    abund_data: dataframe with the abundance data of only one time point, and there 
    should be 10 samples/communities.
    n_sim: number of simulations
    """
    oct_data_list = samples_oct_list(abund_data)
    
    diff_df = pd.DataFrame(columns = ['m', 'cln', 'rel_diff_index', 'abs_diff_index'])
    for i in range(1,26):
        try:
            loaded_sim = pickle.load(open(os.path.join(path_to_simfiles,'TH_metacom_sim%s.pkl'%i), 'rb'))
        except:
            continue
            
        oct_sim_list = samples_oct_list(loaded_sim.sp_abund_df)    

        diff = diff_index_calc(oct_sim_list, oct_data_list)        
        #row = pd.DataFrame([loaded_sim.m, loaded_sim.cln, diff])    
        diff_df = diff_df.append({'m': loaded_sim.m, 'cln':loaded_sim.cln, 'rel_diff_index':diff[0], 'abs_diff_index':diff[1]},ignore_index = True)

    return diff_df


# Let's run the functions. First we have to select samples from one time point of the real data.

# In[32]:

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
abund_data.columns


# Let's run it!

# In[38]:

third_sim_comp = parms_comparator(path_to_simfiles="Third_sim",abund_data=abund_data, n_sim=25)


# In[33]:

third_sim_comp


# No all the simulations are present. Let's order them according to the indices. Let's start with the relative index (third column) measures de difference but as a proportion of the mean of the two octaves.

# In[39]:

third_sim_comp.sort_values('rel_diff_index', axis = 0)


# According to this index, the best is the number 12, with m=0.059142 and cln = 530.
# Let's do the same with the absoulute index, that measures the diference between octaves in number of individuals

# In[40]:

third_sim_comp.sort_values('abs_diff_index', axis = 0)


# The results are consistent with the previous result, both indices are very similar, so it was a bit expected, but I'm glad that that is the case, for simplicity. Let's plot the winner simulation and the real data

# In[44]:

oct_data_list = samples_oct_list(abund_data)
#oct_vect = oct_data_list[11]
#plt.bar(sc.arange(len(oct_vect)), oct_vect)
#plt.xticks(sc.arange(len(oct_vect)), upper[0:len(oct_vect)+1])
#plt.show()

