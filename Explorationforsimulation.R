################################################################
#                     Explorationforsimulation.R
################################################################

# We want to explore which are the best suitable time points of the data, to take the numbers of individuals and run the simulation with them.

# ******************************
# January 27th, 2017
# Imperial College London
# Silwood Park
# Javier Cuadrado Corz
# ******************************

######### SETUP  #######
# Packages
library(ggplot2)
library(dplyr)
library(reshape2)
rm(list = ls())
setwd("~/Documents/Project/Code")

# Paths
pathIn <- "../Data/Processed_data/Samples1K_OTUs100/"
pathOutExploratory <- "../Exploratory_analysis/Factors/"

# Loading data
abiotic <- read.csv(paste(pathIn,"Abiotic.csv", sep = ""))
TH_abund <- read.csv(paste(pathIn,"TH1K_OTUs100.csv", sep = ""), header = TRUE, row.names = 1)
TH_abi <- read.csv(paste(pathIn,"TreeHoles_Abiotic.csv", sep = ""))
R_abund <- read.csv(paste(pathIn,"Rain_Abundances.csv", sep = ""), header = TRUE, row.names = 1)
R_abi <- read.csv(paste(pathIn,"Rain_Abiotic.csv", sep = ""))

#### EXPLORING DATA FOR THE SIMULATION ####
#### Biomass over time
#tree holes
pdf(paste(pathOutExploratory, "TH_Biomass_time.pdf"))
ggplot(data = TH_abi)+
  geom_line(aes(as.Date(TH_abi$R_date), y = TH_abi$Biomass, color = as.factor(TH_abi$Location)))
dev.off()
#abiotic
pdf(paste(pathOutExploratory, "Abund_Biomass_time.pdf"))
ggplot(data = abiotic)+
  geom_line(aes(as.Date(abiotic$R_date[abiotic$R_date == "2015-09-29" or abiotic$R_date == "2015-11-17"]), y = abiotic$Biomass, color = as.factor(abiotic$Location))))
dev.off()
# get rid off some the samples with the enormous biomass values
ggplot(data = abiotic[abiotic$R_date != c("2015-09-29", "2015-09-01"),])+
  geom_line(aes(x = as.Date(R_date), y = Biomass, color = as.factor(Location)))

####### Dataframe with all of the criteria
# We want a time point with all the tree holes samples, all of the biomass values, no very different biomass values and uniform samples regarding drough or flooded
# let's have a look to the highst values of biomass
biomass_df <- abiotic[,c('Samples_shortname', 'R_date', 'Biomass')]
biomass_df <- biomass_df[order(biomass_df$Biomass),]
sort(unique(abiotic$R_date))

# build the criteria table
# number of tree hole samples
n_samples <- group_by(TH_abi, R_date) %>% summarise(n())
# sum of the biomass values of the tree holes for each time point
sum_biomass <- group_by(biomass_df, R_date) %>% summarise(sum(Biomass, na.rm = TRUE))
# number of NAs in biomass for each time point
n_biomass <- TH_abi[,c('R_date', 'Biomass')] %>% group_by(R_date) %>% summarise(sum(is.na(Biomass))) 
# number of Flooded and Drought holes
n_FD <- TH_abi %>% group_by(R_date) %>% count(F_D)
n_FD2 <-dcast(n_FD, formula = R_date ~ F_D, value.var = 'n') # to add them in the general table

biomass_df <- data.frame("n_samples" = n_samples, n_biomass, "Biomass" = sum_biomass, "nFD" = n_FD2)
biomass_df = biomass_df[,c(1,2,4,6,8,9)]
colnames(biomass_df) <- c('R_date',"n_samples","n_biomass", "Biomass", "nD", "nF")

# We want samples with all the tree hole samples and all the biomass values
biomass_df[biomass_df$n_samples == 10 & biomass_df$n_biomass == 0,]
