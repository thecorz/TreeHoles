#!/bin/python

"""
Tree Holes Neutral Simulation
10 treeholes
"""
from __future__ import division
import sys
import csv
import time
import pandas as pd
import scipy as sc
import os
#import lhsmdu
import pickle
import math

pd.__version__

# import logger

class Neutral_sim:
    """ 
        Class for the whole simulation. It will store the info of parameters and everything.
    """
    
    def __init__(self, abund_file, sample_proportion, sim_time, seed, m, cln):
        """ Returns an simulation object with some of the parameters"""
        self.abund_file = abund_file        
        self.sim_time = sim_time * 60 * 60        
        self.seed = seed #random.seed(PBS_index) #HPC
        self.sample_proportion = sample_proportion # between 0 and 1
        self.m = m
        self.cln = cln
        self.sp_abund = []
    
    def individual_lineages_samples(self, abund, abiotic):
        """Assigns the number of ind in the samples from the number of reads in 
        of the samples in the real data. Calculates the total number of ind
        in each community (the population of ind in each tree hole) for the 
        simulation, as a multiple of the number of ind in the sample. The size in num
        of ind of the tree holes is similar then, and doesn't take in account the 
        difference in biomass or volume of the holes"""
        
        TH_abund_tpoint = abund.filter(like = '7_7_15', axis = 1)
        self.N_TH = TH_abund_tpoint.sum(0)
        self.N_TH.sort_index(inplace = True)
        swap = self.N_TH.loc['RNA_ID7_7_15_THS_10']
        self.N_TH.drop('RNA_ID7_7_15_THS_10', inplace = True)        
        self.N_TH = self.N_TH.append(pd.Series(swap))
        self.N_TH.index = ["TH0", "TH1", "TH2", "TH3", "TH4", "TH5", "TH6", "TH7", "TH8", "TH9"]
        self.N_Tot = sum(self.N_TH)
        
        
        self.J_TH = []
        for i in range(10):
            self.J_TH.append(int(max(round(self.N_TH[i] * 1/self.sample_proportion),2)))
        
        self.J_Tot = sum(self.J_TH)        
        
        return 0        

    def individual_lineages_weighted_samples(self, abund, abiotic):
        """Assigns the number of ind in the samples from the number of reads in 
        of the samples in the real data. Calculates the total number of ind
        in each community (the population of ind in each tree hole) for the 
        simulation, as a multiple of the number of ind in the sample. It weights
        the size of the communities in the tree holes with the volume of detritus 
        and the bjomass value"""
        TH_abund_tpoint = abund.filter(like = '7_7_15', axis = 1)
        #number of ind in the sample(lineages) = number of reads in the samples of real data        
        self.N_TH = TH_abund_tpoint.sum(0)
        # order of the tree holes        
        self.N_TH.sort_index(inplace = True)
        swap = self.N_TH.loc['RNA_ID7_7_15_THS_10']
        self.N_TH.drop('RNA_ID7_7_15_THS_10', inplace = True)        
        self.N_TH = self.N_TH.append(pd.Series(swap))
        self.N_TH.index = ["TH0", "TH1", "TH2", "TH3", "TH4", "TH5", "TH6", "TH7", "TH8", "TH9"]
        
        sim.N_Tot = sum(sim.N_TH)
        
        # Calculate the size of the community from the samples with the sample_proportion
        # Here we use the sample proportion in the opposite way. That's why we is inverted.                
        # We use a loop to calculate the each individual because we want to have         
#        equal_weights = []
#        for i in range(10):
#            equal_weights.append(int(max(round(sim.N_TH[i] * 1/sim.sample_proportion),2)))
#        equal_weights_tot = sum(equal_weights)
        equal_weights_tot = sim.N_Tot * 1/sim.sample_proportion

        # Create a matrix with the Vol_detritus and the biomass of all the treeholes
        vol_biomass = abiotic[abiotic.R_date == '2015-07-07'].loc[:,['Samples_shortname', 'Location', 'Vol_detritus', 'Biomass']]
        vol_biomass = vol_biomass.sort_values(by = 'Location')
        vol_biomass.index = range(10)
        # we want to keep the total size all the communities together(J_TH or equal_weights_tot equal
        # but change the relative size of the each tree hole accordingly to their vol of detritus 
        # Here it doesn't make sense multiply by the density because we don't want
        # to estimate the number of individuals in each tree hole
        relative_vol_detritus = vol_biomass.loc[:,'Vol_detritus']/vol_biomass.loc[:,'Vol_detritus'].sum()        
        
        self.J_TH = []
        for i in range(10):
            self.J_TH.append(int(max(round(equal_weights_tot * relative_vol_detritus[i]), 2)))
        
        self.J_Tot = sum(self.J_TH)
    
    def individual_lineages_TH(self, abund, abiotic):
        """Calculates the total number of individuals for the community in the simulation,
        based on the biomass parameter in the abiotic file and Vol_detritus of the hole"""
        # Create a matrix with the Vol_detritus and the biomass of all the treeholes
        vol_biomass = abiotic[abiotic.R_date == '2015-07-07'].loc[:,['Samples_shortname', 'Location', 'Vol_detritus', 'Biomass']]
        vol_biomass = vol_biomass.sort_values(by = 'Location')
        vol_biomass.index = range(10)
        
#            # move the TH10 to the end of the frame
#            moverow = vol_biomass.ix[[1],:]
#            vol_biomass.iloc[10] = moverow.squeeze()
#            vol_biomass = vol_biomass.drop([1,11,12])
#            vol_biomass.index = range(0,10)
#            
        # I looked up the density of detritus in the literature and it's on average
        # 0.2 gr/cm3. the volume of detritus is in cm^3 and the biomass in cfu/gram
        
        self.J_TH = []
        for i in range(10):
            self.J_TH.append(int(max(round(vol_biomass.iloc[i,2]*vol_biomass.iloc[i,3]),2)))
        
        self.J_Tot = sum(self.J_TH)

        TH_abund_tpoint = abund.filter(like = '7_7_15', axis = 1)
        self.N_TH = self.sample_proportion * pd.DataFrame(TH_abund_tpoint.sum(0))
        # order by Tree hole number        
        self.N_TH.sort_index(inplace = True)
        swap = self.N_TH.loc['RNA_ID7_7_15_THS_10']
        self.N_TH.drop('RNA_ID7_7_15_THS_10', inplace = True)        
        self.N_TH = self.N_TH.append(swap).iloc[:,0]
        self.N_TH.index = ["TH0", "TH1", "TH2", "TH3", "TH4", "TH5", "TH6", "TH7", "TH8", "TH9"]
        self.N_Tot = sum(self.N_TH)
        
        return 0       
    
    def lineages_matrix_builder(self):
        """ This is sthe main table, the first column has the total abundance 
        of each lineage. The second column shows were the individual of that lineage
        is in each time step. From the 3rd to the 12th column are the abundances 
        of each lineage but at the bottom of the genealogy"""
        # First column, the total abundance of the lineages
        Ind_Tot = sc.array(sc.repeat(1,self.N_Tot), dtype= '<u4')
        
        # Second column
        b = [i for i in range(0, 10) for item in range(self.N_TH[i])]
        # Put it in an array
        TH_timestep = sc.array(b, dtype = '<u4')
        
        # from the 3rd to the 12th column, the abundance of each lineage in each hole
        TH_origin = sc.zeros([self.N_Tot,10], dtype = '<u4')
        for i in range(TH_origin.shape[0]):
            TH_origin[i,TH_timestep[i]] = 1
        
        # Put the three arrays together
        self.lineages = sc.column_stack((Ind_Tot, TH_timestep, TH_origin))

        return 0
        
    def dispersal_matrix_builder(self):
        """ Builds the dispersal matrix from the distance matrix and the cln"""
        # calculates the denominator out of the loop    
        denominator = (1+(dist_matrix/self.cln))*(self.J_Tot - 1)
        # the DistMatrix is a pandas dataframe this converts it in a numpy array, no so fast but it's not in the loop     
        self.dispersal_matrix = denominator.values    
        for r in range(10):
            for c in range(10):
                if self.cln != 0:
                    self.dispersal_matrix[r,c] = self.J_TH[c]/self.dispersal_matrix[r,c]
                else:
                    self.dispersal_matrix[r,c] = 0
        # For the diagonal of the dispersion matrix, i.e. the individuals that remain in the tree hole where they currently are, this is the probability of remain in the same tree hole.            
        for i in range(10):
            self.dispersal_matrix[i,i] = 1 - (sum(sc.delete(self.dispersal_matrix[i,:],i)))
        
        return 0
    
    def metacommunity_builder(self):
        """Builds the metacommunity from the abundance of the samples"""
        total_OTU_abundTH = abund.sum(axis = 1)
        total_ind_abundTH = total_OTU_abundTH.sum()
        self.metacom = total_OTU_abundTH / total_ind_abundTH
        return 0       
        
    def dispersion_event(self, start, rndm, ind):
        """ We are going to decide where the individual is going to disperse accoring to the random number we originated. The probabilies of dispersion to each of the tree holes, including the probability of no dispersion add up to 1. These probability are going to be our bins. If we consider this probabilities as bins, the individual is going to migrate to the TH that has a range that includes rndm. In order to do this process as quickly as possible, basically we have to classify the rndm into 10 bins. I think that one of the most efficient ways of doing this is dividing the possibilities in two, making a classification tree. In order to do so we will use nested if else statements."""
         
        # classification tree with nested if and else statements. I think it can be done in a faster way working out the acumulated probability matrix, where each element of the matrix is the sum of the previous ones by row. #profiling #improve
        first_half = sum(self.dispersal_matrix[start, 0:5])     
        if rndm <= first_half:
            if rndm <= sum(self.dispersal_matrix[start, 0:3]):
                if rndm <= self.dispersal_matrix[start, 0]:
                    self.lineages[ind, 1] = 0
                else:
                    if rndm <= sum(self.dispersal_matrix[start, 0:2]):
                        self.lineages[ind, 1] = 1
                    else:
                        self.lineages[ind, 1] = 2
            else:
                if rndm <= sum(self.dispersal_matrix[start, 0:4]):
                    self.lineages[ind, 1] = 3
                else:
                    self.lineages[ind, 1] = 4
        else:
            if rndm <= first_half + sum(self.dispersal_matrix[start, 5:8]):
                if rndm <= first_half + self.dispersal_matrix[start, 5]:
                    self.lineages[ind, 1] = 5
                else:
                    if rndm <= first_half + sum(self.dispersal_matrix[start, 5:7]):
                        self.lineages[ind, 1] = 6
                    else:
                        self.lineages[ind, 1] = 7
            else:
                if rndm <= first_half + sum(self.dispersal_matrix[start, 5:9]):
                    self.lineages[ind, 1] = 8
                else:
                    self.lineages[ind, 1] = 9
        # put the final destination of the migration in an object
        return self.lineages[ind, 1]

    def migration_event(self, ind):
        """ Migration from the metacommunity event. Similar to speciation"""
        
        # we randomly choose an element of our metacommunity, weighted according to his relative abundance. 
        # relative_OTU_abund it's index is the sps and the only column is the relative abundance. 
        # random.choice takes weights with the parameter p            
        ind_metacom = sc.random.choice(self.metacom.index, p = self.metacom)
        # In lineages we use the second column to store the position of the individual in each simulation step,
        # so once we stablish that a imigration occurs we don't need it anymore and we use that position to 
        # store the identity of the taxon.
        self.lineages[ind, 1] = ind_metacom
        # here I could transform the array slice into a list I don't know what is more efficient. 
        # I thought this could be the memmory leak (if it was pandas but it's an array)            
        self.sp_abund.append(self.lineages[ind])
        self.lineages = sc.delete(self.lineages, ind, 0)  #profiling
        
        return 0
            
    def prob_coalescence(self, start, end):
        """ Probability of coalescence"""
        n_lineages_end = sum(self.lineages[:, 1] == end)
        if start == end:
            coa_prob = (n_lineages_end-1)/(self.J_TH[end]-1) 
        else:
            coa_prob = (n_lineages_end)/(self.J_TH[end])
        
        # coa_prob = 0.01  #debug
        return coa_prob
        
    def coalescence_event(self, coa_prob, ind, end):
        """Coalescence event """        
                       
        # choose a random ind from the end treehole to coalesce with, the indv has to be different from the first ind 
        range_ind = sc.flatnonzero(self.lineages[:,1] == end)  # range: the positions of the indv to choose from
        range_ind = sc.delete(range_ind,sc.flatnonzero(range_ind == ind))  # get rid of the first ind in range
        ind2 = sc.random.choice(range_ind)
        
        self.lineages[ind, 0] = self.lineages[ind, 0] + self.lineages[ind2, 0]
        self.lineages[ind, 2:12] = self.lineages[ind, 2:12] + self.lineages[ind2, 2:12]
        self.lineages = sc.delete(self.lineages, ind2, 0)

        return 0        

    def coalescence_simulation(self, start_time): #profiling
        """Main loop of the simulation, where migration, coalescence and specitation take place """
        
        while self.N_Tot > 1 and ((start_time + self.sim_time) > time.time()):
            
            # choose an individual from the matrix (a lineage)
            ind = sc.random.choice(sc.flatnonzero(self.lineages[:,0]))#profiling
            # choose a random number from 0 to 1
            rndm = sc.random.uniform()
    
            # Dispersion        
            start = self.lineages[ind, 1]
            end = self.dispersion_event(start = start, rndm = rndm, ind = ind)
                       
            # prob of coalescenc
            coa_prob = self.prob_coalescence(start, end)

            # migration from the metacommunity (speciation)
            if sc.random.uniform() <= self.m:
                self.migration_event(ind)
            #coalescence
            elif sc.random.uniform() <= coa_prob and sum(self.lineages[:, 1] == end) > 1:
                self.coalescence_event(coa_prob, ind, end)

            self.N_Tot = len(self.lineages[:, 1])
        
        if (start_time + self.sim_time) > time.time():
            # append the last individual to sp_abund
            self.lineages = self.lineages.flatten()    
            self.lineages[1] = sc.random.choice(self.metacom.index, p = self.metacom)   
            self.sp_abund.append(self.lineages)
        self.sp_abund_df = pd.DataFrame.from_records(self.sp_abund, columns = ["Tot_abund", "ID", "TH0", "TH1", "TH2", "TH3", "TH4", "TH5", "TH6", "TH7", "TH8", "TH9"])
        self.sp_abund_df = self.sp_abund_df.groupby(by = 'ID').sum()

        return 0
    
    def save_results(self):
        self.sp_abund_df.to_csv("Sim.csv", header = True, index = True)
        variables = [self.cln, self.m, self.J_TH, self.N_TH, self.sample_proportion]
        with open('variables%s.csv'%PBS_index, 'wb') as myfile:
            wr = csv.writer(myfile, quoting=csv.QUOTE_MINIMAL)
            wr.writerow(variables)

        params = [self.cln, self.m, self.sample_proportion]
        with open('params.csv', 'wb') as myfile:
            wr = csv.writer(myfile, quoting=csv.QUOTE_MINIMAL)
            wr.writerow(params)
#        os.system("cat params.csv Sim.csv > Sim_$PBS_ARRAY_INDEX.csv")
        os.system("cat params.csv Sim.csv > BestSim.csv")
        os.system("rm sim.csv")
        os.system("rm params.csv")
#def latin_hypercube_sampling(m_min, m_max, cln_min, cln_max, n_param_samples):
#    """
#    Samples efficiently and completely from the parameter space specified. Returns a list of pairs of parameters which
#    simulations should run with.
#    :param m_min: the minimum value of the first dimension
#    :param m_max: the maximum value of the first dimension
#    :param cln_min: the minimum value of the second dimension
#    :param cln_max: the maximum value of the second dimension
#    :param n_param_samples: the number of samples to draw in each dimension
#    :return: a list of parameter pairs generated using the lhsmdu module.
#    """
#    samples = lhsmdu.sample(2, n_param_samples, randomSeed=1001)
#    # debug assertions
#    assert (m_min < m_max)
#    assert (cln_min < cln_max)
#    
#    samples[0] = (m_max - m_min) * samples[0] + m_min
#    samples[1] = (cln_max - cln_min) * samples[1] + cln_min
#    output_list = []
#    # Now create our pairs of parameters
#    for i in range(0, n_param_samples, 1):
#        output_list.append([samples[0, i], samples[1, i]])
#    
#    return (output_list)

def choose_parameters_log(m_min, m_max, cln_min, cln_max, n_param_samples):
    """Distributes the parameters evenly between the ranges, but the paramaters
    are not rounded."""
    m = sc.linspace(math.log10(m_min), math.log10(m_max), n_param_samples)
    cln = sc.linspace(math.log10(cln_min), math.log10(cln_max), n_param_samples)
    
    output_list = []
    for i in range(0, n_param_samples, 1):
        for j in range(0, n_param_samples, 1):        
            output_list.append([math.pow(10,m[j]), math.pow(10,cln[i])])
    
    return (output_list)
    
choose_parameters_log(0.1, 0.001, 1, 0.1, 5)
  
def choose_parameters(m_min, m_max, cln_min, cln_max, n_param_samples):
    """Distributes the parameters evenly between the ranges"""
    m = sc.linspace(m_min, m_max, n_param_samples)
    cln = sc.linspace(cln_min, cln_max, n_param_samples)
    
    output_list = []
    for i in range(0, n_param_samples, 1):
        for j in range(0, n_param_samples, 1):        
            output_list.append([m[j], cln[i]])
    
    return (output_list)
    
#choose_parameters(0.1, 0.001, 1, 0.1, 5)
  

if  __name__ =='__main__':
    # start measuring time
    start_time = time.time()
    # set working directory
    os.chdir('/home/cmee19/Documents/Project/Sim_TH/Metacom_sim/')
       
    #TH_names = ["TH0", "TH1", "TH2", "TH3", "TH4", "TH5", "TH6", "TH7", "TH8", "TH9"]
    
#    pathWork = os.path.join(os.environ.get('WORK'), 'TreeHoles/Simulation')   
#    abund_pathfile = os.path.join(pathWork, 'TreeHoles_Abundances2.csv')    
#    abund = pd.read_csv(abund_pathfile)
#    dist_matrix = pd.read_csv(os.path.join(pathWork,'Circular_dist_matrix.csv'), header = 0, index_col = 0)
#    abiotic = pd.read_csv(os.path.join(pathWork,'TreeHoles_Abiotic.csv'), header = 0)   
    abund_pathfile = '../../Data/Processed_data/Samples1K_OTUs100/TH_abund.csv'
    abund = pd.read_csv(abund_pathfile)
    dist_matrix = pd.read_csv('../../Data/Processed_data/Circular_dist_matrix.csv', header = 0, index_col = 0)
    abiotic = pd.read_csv('../../Data/Processed_data/Samples1K_OTUs100/TreeHoles_Abiotic.csv', header = 0)
# Take the tag for the paralize processing as an argument when you run the script   
#    if len(sys.argv) < 2:
#        PBS_index = 1
#    else:
#        PBS_index = int(sys.argv[1])

    PBS_index = 1
    #param_list = [choose_parameters_log(0.0001, 0.1, 1, 1000, 25), choose_parameters_log(0.0001, 0.1, 1, 1000, 25), choose_parameters_log(0.0001, 0.1, 1, 1000, 25)]
    #param_list = [item for sublist in param_list for item in sublist]    
    param_list = [(0.0001, 1)]
        
    if PBS_index > 0 and PBS_index <= 25:
        sample_proportion = 0.1
    if PBS_index > 25 and PBS_index <= 50:
        sample_proportion = 0.001
    if PBS_index > 50 and PBS_index <= 75:
        sample_proportion = 0.0001
        
        
    sim = Neutral_sim(abund_pathfile, sim_time=56, sample_proportion = sample_proportion, seed = PBS_index, m = param_list[PBS_index-1][0], cln = param_list[PBS_index-1][1])    
    sim.individual_lineages_samples(abund, abiotic)
    sim.lineages_matrix_builder()
    sim.dispersal_matrix_builder()
    sim.metacommunity_builder()
    sim.coalescence_simulation(start_time)
    sim.save_results()
    pickle.dump(sim, open('TH_metacom_sim%s.pkl' %PBS_index, 'wb'))
    #loaded_sim = pickle.load(open('TH_metacom_sim%s.pkl' %PBS_index, 'rb'))
