####################################################################
##### Implementation of Sean's code for the tree holes simulation
####################################################################


#### Load Data ####

rm(list=ls())    # clear up workspace from any previous functions, variables, etc.

setwd("~/Documents/Project/Code/Mechanistic_model")

DistMatrix <- read.csv("../../Data/Processed_data/TH_Spatialdist_matrix.csv")
abiotic <- read.csv("../../Data/Processed_data/Samples1K_OTUs100/Shortnames/Abiotic.csv")

# Parameters values
Density = 100; Sample_scale = 0.001; generation = 2; v = 0.01; cln = 100; tau=-1; verbose=1; rand_seed = 1 


# TreeHolesSim<-function(Density,Sample_scale,Collision,generation,v,cln,tau=-1,verbose=1,rand_seed){

    
#   v - probability of speciation and is the green letter nu (which looks like v) in the main text
  #   cln - colonisation parameter c in the main text
  #   verbose - if =1 then full data will be printed to screen and saved to a file (see below)
  #   seed - alows user to optionally set the random seed
  #####################################  


  ptm<-proc.time()[3]
  # start timer to say how long run took

  ##################################################
  #### Set up variables and parameters ####

  # Set the sizes of each Tree Hole measured in numbers of individuals
  #if it's too small though it will be assigned a value of 2 for the purpose of small scale testing
  # The number of ind is calculated: Volume * Density
  # Function max so it's going to use a minimum of 2 indv. this is in case the Density is set to a very low number.
  
  # Take just one time point: 1-9-2015 because it has the most of the samples
  Vol_Biomass <- unique(abiotic[abiotic$R.date == "15-09-01", c(2,15,23,36)])
  Vol_Biomass <- Vol_Biomass[order(Vol_Biomass$Location),]
  Vol_Biomass <- Vol_Biomass[c(1,3:10,2),-1]
  
  
  # TH1 <- max(2,round(Vol_Biomass[1,2]*Vol_Biomass[1,3]))
  # TH2 <- max(2,round(Vol_Biomass[2,2]*Vol_Biomass[2,3]))
  # TH3 <- max(2,round(Vol_Biomass[3,2]*Vol_Biomass[3,3]))
  # TH4 <- max(2,round(Vol_Biomass[4,2]*Vol_Biomass[4,3]))
  # TH5 <- max(2,round(Vol_Biomass[5,2]*Vol_Biomass[5,3]))
  # TH6 <- max(2,round(Vol_Biomass[6,2]*Vol_Biomass[6,3]))
  # TH7 <- max(2,round(Vol_Biomass[7,2]*Vol_Biomass[7,3]))
  # TH8 <- max(2,round(Vol_Biomass[8,2]*Vol_Biomass[8,3]))
  # TH9 <- max(2,round(Vol_Biomass[9,2]*Vol_Biomass[9,3]))
  # TH10 <- max(2,round(Vol_Biomass[10,2]*Vol_Biomass[10,3]))
  
  densitytest <- 100
  
  TH1 <- max(2,round(Vol_Biomass[1,2]*densitytest))
  TH2 <- max(2,round(Vol_Biomass[2,2]*densitytest))
  TH3 <- max(2,round(Vol_Biomass[3,2]*densitytest))
  TH4 <- max(2,round(Vol_Biomass[4,2]*densitytest))
  TH5 <- max(2,round(Vol_Biomass[5,2]*densitytest))
  TH6 <- max(2,round(Vol_Biomass[6,2]*densitytest))
  TH7 <- max(2,round(Vol_Biomass[7,2]*densitytest))
  TH8 <- max(2,round(Vol_Biomass[8,2]*densitytest))
  TH9 <- max(2,round(Vol_Biomass[9,2]*densitytest))
  TH10 <- max(2,round(Vol_Biomass[10,2]*densitytest))
  
  
  # calculate the global community size as the sum of all communities
  J_Tot<-sum(TH1,TH2,TH3,TH4,TH5,TH6,TH7,TH8,TH9,TH10)

  # this tells us how many lineages there are being sampled per land mass
  N_TH1<-max(2,round(TH1*Sample_scale))
  N_TH2<-max(2,round(TH2*Sample_scale))
  N_TH3<-max(2,round(TH3*Sample_scale))
  N_TH4<-max(2,round(TH4*Sample_scale))
  N_TH5<-max(2,round(TH5*Sample_scale))
  N_TH6<-max(2,round(TH6*Sample_scale))
  N_TH7<-max(2,round(TH7*Sample_scale))
  N_TH8<-max(2,round(TH8*Sample_scale))
  N_TH9<-max(2,round(TH9*Sample_scale))
  N_TH10<-max(2,round(TH10*Sample_scale))
  
  # 2a: first set up parameters for coalescence simulation
  
  set.seed(rand_seed)
  
  # the global lineage count as the sum of all lineages
  N_Tot<-sum(N_TH1, N_TH2, N_TH3, N_TH4, N_TH5, N_TH6, N_TH7, N_TH8, N_TH9, N_TH10)
  
  # create data frame to store how many lineages there are in each continent
  lineage_matrix<-data.frame(N_TH1, N_TH2, N_TH3, N_TH4, N_TH5, N_TH6, N_TH7, N_TH8, N_TH9, N_TH10)
  
  # this variable will hold the species richness of the system - starting at 0 it is incrememted whenever speciation occurs on a lineage that is known to survive and be sampled. Speciated lineages are then removed from the simulations.
  sp_abund <- as.data.frame(matrix(ncol = 10))
  colnames(sp_abund) <- c("TH1", "TH2", "TH3", "TH4", "TH5", "TH6", "TH7", "TH8", "TH9", "TH10")
  
#####################################################################
#### Main DF #####
  
  # build a dataframe for information on the lineages for coalescnece - each row is a lineage, the second column gives the continent that lineage is in. The first column is alwyas 1 in this implementation, but could be edited later to keep track of species abundacnes as well as species richnesses 
  Lineages <- rep(1,N_Tot)
  TH <- c((rep(1,N_TH1)),(rep(2,N_TH2)),(rep(3,N_TH3)),(rep(4,N_TH4)),(rep(5,N_TH5)),(rep(6,N_TH6)),(rep(7,N_TH7)),(rep(8,N_TH8)),(rep(9,N_TH9)),(rep(10,N_TH10)))
  DF<-data.frame(Lineages,TH)
  
  # The abundances of the different lineages but in the present day, at the bottom of the genealogy. The number of lineages for each simulation time step is kept in the vector "Lineages" in DF 
  present <- matrix(0, ncol = 10, nrow = nrow(DF))
  colnames(present) <- c("TH1", "TH2", "TH3", "TH4", "TH5", "TH6", "TH7", "TH8", "TH9", "TH10")
  
  for (i in 1:nrow(DF)){
    present[i,DF[i,2]] <- 1
  }
  
  DF <- cbind(DF, present)
  
  # J_TH will be used to index the community size that each land mass has 
  J_TH <- data.frame(TH1,TH2,TH3,TH4,TH5,TH6,TH7,TH8,TH9,TH10)

  # total_time<-0
  # total time of simulation in simulation time steps
  
##############################################################
##### Dispersal Matrix #####
  # 2b - dispersal: determine where the current time sits within the geological history so that the correct distances can be deteremined.
  
  denominator<-(1+(DistMatrix/cln))*(J_Tot-1) # It's a constant so it calculates it out the loop once better than calculate it in the loop many times(more time)
  # begin to form the colonisation/dispersal probabilities by calculating the dividend.
  
  # re-normalise the dispersal probabilities matrix 
  dispersal_matrix<-denominator
  for(r in 1:10){
    for(c in 1:10){
      if(cln!=0){
        dispersal_matrix[r,c]<-J_TH[c]/dispersal_matrix[r,c]
        # each respective landmass then has its community size divided by the dividend thus yielding a probability. If cln==0 though, no colonisation occurs and it is instead assigned a value of zero.
      } else {
        dispersal_matrix[r,c]<-0
      }
    }
  }
  
  for(i in 1:10){
    dispersal_matrix[i,i]<-1-(sum(dispersal_matrix[i,-i]))
    # determines the probability of an individual remaining within it's starting land mass
  }

#################################################################  
#### WHILE LOOP ####
    
  #while(N_Tot > 1){
    # main while loop of the phase 2 part of the simulation.

    # Reorder DF according to the TH the lineage is in
    DF <- DF[order(DF$TH),]
    
    #### Should the TH that get to only one indv carry on in the simulation?
    
    # # Inlclude only TH that have more than 1 linneage in the group for choosing an individual
    # 
    # DFselect <- DF
    # 
    # if(length(lineage_matrix))
    
    
    #### Choose individual ####
    # 2c -figure out which individual dies and is replaced (in reverse time because this is coalescence)
    
    N_Tot <- length(DF[,1])
    # recount the current number of lineages remaining in system
    
    individual <- sample(seq(1,N_Tot,1),1)
    # use the lineage count to randomly pick an individual lineage on which dispersal may act
    
    ind1 <- rownames(DF)[individual]
    
    ##############################################################
    #### Dispersion ####
    rndm<-runif(1,0,1)
    # random number used to figure out where that individal chosen will disperse to (if anywhere)
    # dispersal_matrix
    start<-DF[individual,2]
    #starting location of the individual
    
    if(rndm<=dispersal_matrix[start,1]){DF[individual,2] <- 1}
    
    for (i in 2:10){
      if(rndm <= sum(dispersal_matrix[start,1:i]) && rndm > sum(dispersal_matrix[start,1:i-1])){DF[individual,2] <- i}
    }
    
    End<-DF[individual,2]
    #end location of individual
    
    # 2d - lineage movement (if applicable)

    if(start!=End){
        # If the end location is different from the starting location, then the respective communities need to have their lineage counts adjusted
      lineage_matrix[start]<-lineage_matrix[start]-1
      lineage_matrix[End]<-lineage_matrix[End]+1
    }
    #   lineage_matrix
    
    N_Tot<-sum(lineage_matrix)
    #   N_Tot Total number of lineages in all Tree Holes 
    
    #### Prob speciation and coalescence  ####
    # 2e - speciation and coalescence probabilities to be calculated

    if(start==End){
      prob_coa <- (lineage_matrix[End]-1)/(J_TH[End]-1) # Apparently the probability of coalescence is the number of lineages in that continent divided by the number of individuals in that continent
      # if the start == End then there is one less individual that our random chosen individual could coalesce with
    } else {
      prob_coa<-lineage_matrix[End]/J_TH[End]
    }

    prob_spec<-v
    #provided dispersal doesn't influence spec prob
    
    # dispersal influenced speciation with tau
    
    if(start!=End){
      # possibly do vicariance speciation
      if (tau > 0) { # tau is a parameter of the model: the probability of speciation due to vicariance. If it's set to a negative number 1then no vicariance is added to the the model.  
        # speciation depends on distance and tau
        prob_spec<-v+((1-v)*(DistMatrix/(DistMatrix+tau)))
      } else if (tau == 0) {
         # speciation guaranteed because tau = 0 and dispersal is cross-continent.
         prob_spec<- 1
      }
      # note that if tau < 0 then nothing changes here so we have no vicariance speciation
    }
    
    #### Perform Speciation and coalescence ####
    # 2f - do speciation and coalescence as needed
    prob_coa <- 0.01 #debug
   if(runif(1,0,1)<=prob_spec){
      # using a random number, test if speciation or coalescence occured.
      sp_abund <- rbind(sp_abund, DF[individual,3:12])
      DF<-DF[-individual,]
      lineage_matrix[End]<-lineage_matrix[End]-1
   }
   if(runif(1,0,1)<=prob_coa){
      #individual2 <- sample(seq(1,lineage_matrix[[End]]), 1)
      range <- which(DF$TH == End)
      
      # Choose a second individual to coalesce with
      # Which position has the first ind in the TH[End]. the var "indvidual" holds the position of the first individual in DF, i.e. in the whole lot of TH's.
      individual1_TH <- individual - (min(which(DF$TH == End))-1)
      
      # Take out the individual1 from the range of individual to choose the second individual from
      
      range <- range[-(individual1_TH)]
      
      # Choose the second individual among these
      
      individual2 <- sample(range, 1)
      
      
      # individual2 <- sample(which(DF$TH == End), 1)
      # # If the first chosen and the second ind are the same, sample randombly 
      # if (individual == individual2){
      #   while (individual == individual2){
      #     if (length(which(DF$TH == End)) == 1) {break}
      #     individual2 <- sample(which(DF$TH == End), 1)    
      #   }
      # }
      if (length(range) > 1){     
      DF[individual,1] <- DF[individual, 1] + DF[individual2, 1]
      DF[individual,3:12] <- DF[individual, 3:12] + DF[individual2, 3:12]
      DF <- DF[-individual2,]
      lineage_matrix[End]<-lineage_matrix[End]-1
      }
    }
    #total_time<-total_time+((1/N_Tot)*(J_Tot/(2*generation)))

    if (verbose == 1) {
      print("Linneages:")
      print(lineage_matrix)
      print("N_Tot")
      print(N_Tot)
      # print("Species abundance")
      # print(sp_abund)
    }
    
    
  }# Termina el while loop
  # The first column are NA so I just get rid off them. Im sure there is a more elegant way to do this.
  sp_abund <- sp_abund[-1,]
  
  if (verbose == 1) {
      
      total_time_secs <-proc.time()[3]-ptm
      run_time<-total_time_secs/60
      # figure out how long simulations took
      
      print("simulations completed")
      print("random seed =")
      print(rand_seed)
      # print("simualtion run time (minutes) = ")
      # print(as.numeric(run_time))
      print("species richness result =")
      print(sp_abund)
      print("outputting data to .RDA file named")
      filenametemp <- paste(rand_seed,"GEO_sample",Sample_scale,"_v",v,"_c",cln,".Rdata",sep="")
   
      print(filenametemp)
      
      save(rand_seed, sp_abund, file=filenametemp)

  }
  
#   # return the species richnes result
#   return(sp_rich)
# #}
# 
# # Test simulation run
# GeoSim(Density=1/100000,Sample_scale=0.1,Collision=10*10^6,generation=50,v=0.1,cln=0,rand_seed = 1)
# 


#####################
# FOR USE WITH CLUSTER
#
# iter<-as.numeric(Sys.getenv("PBS_ARRAY_INDEX")) # PBS array 1-4699
#   #   iter<-299 # for local testing. At current settings will result in 4 Ma being the time point simulated 
# point<-((iter-iter%%100)/100)  # 100 repetitions will be perform for each point in time
# time_point<-(point-1)*4   # -1 enables you tart prior to pangaea's break up for nicer looking figures.
#
# GeoSim(Density=1/100000,Sample_scale=0.1,Collision=time_point*10^6,generation=50,v=0.1,cln=0,rand_seed = iter)

####################



