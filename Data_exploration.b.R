################################################################
#                     DATA_EXPLORATION.R
################################################################

# Exploration of the simulated data. Similar script than for the real data


# ******************************
# January 27th, 2017
# Imperial College London
# Silwood Park
# Javier Cuadrado Corz
# ******************************

## Packages

library(corrgram)
library(ggplot2)
library(sads)
library(dplyr)

rm(list = ls())
setwd("~/Documents/Project/Code")

# Paths

pathIn <- "../Neutral_Sim/Metacom_sim/Results/"
pathOutExploratory <- "../Neutral_Sim/Metacom_sim/Results/"

# Loading data

abiotic <- read.csv(paste(pathIn,"First_run.csv", sep = ""))


TH_abund <- read.csv(paste(pathIn,"TreeHoles_Abundances.csv", sep = ""), header = TRUE, row.names = 1)

TH_abi <- read.csv(paste(pathIn,"TreeHoles_Abiotic.csv", sep = ""))

R_abund <- read.csv(paste(pathIn,"Rain_Abundances.csv", sep = ""), header = TRUE, row.names = 1)

R_abi <- read.csv(paste(pathIn,"Rain_Abiotic.csv", sep = ""))


#########################################
#####     DESCRIPTION OF THE DATA   #####

## Number of NA's
NAs_TH <- apply(TH_abi, 2, function(x) sum(is.na(x)))
write.csv(NAs_TH, paste(pathOutExploratory,"NAs_TH.csv", sep = ""))
NAs_Rain <- apply(R_abi, 2, function(x) sum(is.na(x)))
write.csv(NAs_Rain, paste(pathOutExploratory,"NAs_Rain.csv", sep = ""))



## Number of 0's
zerosTH <- apply(TH_abund==0, 2, sum)
write.csv(zerosTH, paste(pathOutExploratory,"zeros_TH.csv", sep = ""))
zerosRain <- apply(R_abund==0, 2, sum)
write.csv(zerosRain, paste(pathOutExploratory,"zeros_Rain.csv", sep = ""))

svg(paste(pathOutExploratory,"zeros.pdf", sep = ""))
par(mfrow = c(2,2))
plot(zerosTH)
hist(zerosTH)
plot(zerosRain)
hist(zerosRain)
dev.off()


#####  Genus's and Reads per sample ####                         
#### Tree Holes
hist(TH_abi$TotalReads)
hist(TH_abi$nGenus)
# Order the samples in relation to the number of reads
Ordered_TH_abi <- TH_abi[order(TH_abi$TotalReads),]

#### Rain

hist(R_abi$TotalReads)
hist(R_abi$nGenus)

##### problems with R.date ######
# There are samples where the date in the name of the sample doesn't match with the date in the column R.date

# datesfromname <- regmatches(TH_abi$ID.OTU, gregexpr("\\d\\d\\.\\d\\d\\.\\d\\d", TH_abi$ID.OTU))
# 
# compareDates <- cbind(datesfromname, as.character.Date(TH_abi$R.date))
# write.csv(compareDates, paste(pathOutExploratory,"CompareDates.csv", sep = ""))

# The good dates are in the R.date column

####### Number of samples in each time point for Processed data ##########

Samples_time_TH <- table(TH_abi$R.date[TH_abi$Sample.type == "THS"])

# for rain is a littel bit more complicated because there are 0 counts, that the function table wont show unless is transformed in a factor, and indicate explicitly the dates (levels =)
Samples_time_R <- table(factor(R_abi$R.date[R_abi$Sample.type == "RAIN"], levels = unique(TH_abi$R.date)))
# Because the dates are a factor the order is different. We order them with sort
Samples_time_R <- Samples_time_R[sort(names(Samples_time_R))]

# ready to put them in a df
Samples_time <- data.frame(TreeHoles = Samples_time_TH, Rain = Samples_time_R, row.names = 1)
# Samples_time$Rain.Var1 <- as.Date(Samples_time$Rain.Var1, format = "%d.%m.%Y")
# Samples_time <- Samples_time[order(Samples_time$Rain.Var1),]
# Samples_time <- Samples_time[c(2,1,3)]
write.csv(Samples_time, paste(pathOutExploratory,"ProcessedSamples_time.csv", sep = ""))

svg(paste(pathOutExploratory,"ProcessedSamples_time.svg", sep = ""))
ggplot(data = Samples_time)+ 
  geom_step(aes(x = as.Date(rownames(Samples_time)), y = Samples_time$TreeHoles.Freq, color = "Tree Holes"))+
  geom_step(aes(x = as.Date(rownames(Samples_time)), y = Samples_time$Rain, color = "Rain"))+
  scale_y_continuous(breaks = seq(0, 10, by = 1))+
  scale_x_date(breaks = as.Date(rownames(Samples_time)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

##### Samples per time point and hole ####

svg(paste(pathOutExploratory,"ProcessedSamples_time_hole.svg", sep = ""))
ggplot(data = abiotic)+
  geom_point(aes(x = as.Date(abiotic$R.date), y = abiotic$Location, size = abiotic$TotalReads, color = abiotic$nGenus))+
  scale_x_date(breaks = as.Date(abiotic$R.date))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#### Huge differences in biomass accross samples ####
# 30 degrees of magnitude of difference between samples is a lot. Damian told me that Biomass refered to colony forming units per gram of sample. Wondering if instead per gram is total in the hole, taking in account the volume of the hole.

biomass <- lm(abiotic$Biomass ~ abiotic$Volume)
summary(biomass) 
library(lmer)


plot(abiotic$Biomass ~ abiotic$Volume)

###########################################################
############## FACTORS OVER TIME #########################

svg(paste(pathOutExploratory,"TH_Temperature_time.svg", sep = ""))
ggplot(data = TH_abi[TH_abi$Location == 3,])+
  geom_line(aes(x = as.Date(TH_abi$R.date[TH_abi$Location == 3]), y = TH_abi$Air.Temp[TH_abi$Location == 3], color = "Air.Temp.Mean")) +
  geom_line(aes(x = as.Date(TH_abi$R.date[TH_abi$Location == 3]), y = TH_abi$Air.Temp.max[TH_abi$Location == 3], color = "Air.Temp.Max"))+
  geom_line(aes(x = as.Date(TH_abi$R.date[TH_abi$Location == 3]), y = TH_abi$Air.Temp.min[TH_abi$Location == 3], color = "Air.Temp.Min"))+
  geom_line(aes(x = as.Date(TH_abi$R.date[TH_abi$Location == 3]), y = TH_abi$Grass.Temp[TH_abi$Location == 3], color = "Grass.Temp.Mean"))+
  geom_line(aes(x = as.Date(TH_abi$R.date[TH_abi$Location == 3]), y = TH_abi$Grass.Temp.max[TH_abi$Location == 3], color = "Grass.Temp.max"))+
  geom_line(aes(x = as.Date(TH_abi$R.date[TH_abi$Location == 3]), y = TH_abi$Grass.Temp.min[TH_abi$Location == 3], color = "Grass.Temp.min"))+
  geom_line(aes(x = as.Date(TH_abi$R.date[TH_abi$Location == 3]), y = TH_abi$Soil.Temp2[TH_abi$Location == 3], color = "Soil.Temp2"))+
  geom_line(aes(x = as.Date(TH_abi$R.date[TH_abi$Location == 3]), y = TH_abi$Soil.Temp4[TH_abi$Location == 3], color = "Soil.Temp4"))+
  geom_line(aes(x = as.Date(TH_abi$R.date[TH_abi$Location == 3]), y = TH_abi$Grass.Temp.max[TH_abi$Location == 3], color = "Grass.Temp.max"))+
  scale_x_date()
dev.off()

# It's not possible to plot two data sets with different scales
# I can use the basic plot


TH_abi_ordered <- TH_abi[order(TH_abi$R.date),]


svg(paste(pathOutExploratory,"TH_Temp_Rain_time.svg", sep = ""))
par(mar = c(5,5,2,5))
with(TH_abi_ordered, plot(x = as.Date(TH_abi_ordered$R.date[TH_abi_ordered$Location == 3]), TH_abi_ordered$Rainfall[TH_abi_ordered$Location == 3], type="l", col="blue", ylab = "Rainfall"))

par(new = T)
with(TH_abi_ordered, plot(x = as.Date(TH_abi_ordered$R.date[TH_abi_ordered$Location == 3]), TH_abi_ordered$Air.Temp[TH_abi_ordered$Location == 3], type="l", col="2", axes=F, xlab=NA, ylab=NA))
axis(side = 4)
mtext(side = 4, line = 3, 'Mean Air Temperature', col = 2)
dev.off()

###########################################################
############## EXPLORATORY ANALYIS ########################

#####  All against all ploting #####
### TREE HOLES
# abiotic data temperatures and rain
svg(paste(pathOutExploratory,"Correlation_THabioticVar.svg", sep = ""))
pairs(TH_abi[,17:26])
dev.off()
# other var
svg(paste(pathOutExploratory,"Correlation_THVar.svg", sep = ""))
pairs(TH_abi[,7:16])
dev.off()

##### correlation between variables  #####
Corr_var <- cor(TH_abi[c(7,10,11,13:25,27,28)],use = "complete.obs")
write.csv(Corr_var, paste(pathOutExploratory,"CorrVar_TH", sep = ""))

##### Correlogram
svg(paste(pathOutExploratory,"Correlogram_TH.svg", sep = ""))
corrgram(TH_abi[c(7,10,11,13:25,27,28)], lower.panel = panel.pts, upper.panel = panel.conf, diag.panel = panel.minmax,  cex.labels = 1) 
dev.off()

###############################################################
#####         SPECIES RICHNESS                           ######

##### Sps Abundance Distributions ####

# # I'll take one the smallish samples (the smallest is 0 reads) and one of the largest sample (with more reads)
# 
# # F.D.RNA.P1.D3.C.05.ID10.2.15.THS.6.16S.4 has 1204 reads 115 OTU's
# TH_abi[TH_abi$ID.OTU == "F.D.RNA.P1.D3.C.05.ID10.2.15.THS.6.16S.4",]
# 
# svg(paste(pathOutExploratory,"Sps_dist1.svg", sep = ""))
# hist(TH_abund$F.D.RNA.P1.D3.C.05.ID10.2.15.THS.6.16S.4[TH_abund$F.D.RNA.P1.D3.C.05.ID10.2.15.THS.6.16S.4!= 0])
# dev.off()

# F.D.RNA.P3.B4.H.01.ID09.06.15.THS.8.16S.4 has 86751 reads 1693 OTU's
TH_abi[TH_abi$ID.OTU == "F.D.RNA.P3.B4.H.01.ID09.06.15.THS.8.16S.4",]

svg(paste(pathOutExploratory,"Sps_distribution_9_6_15_TH8.svg", sep = ""))
hist(TH_abund$F.D.RNA.P3.B4.H.01.ID09.06.15.THS.8.16S.4[TH_abund$F.D.RNA.P3.B4.H.01.ID09.06.15.THS.8.16S.4!=0])
dev.off()
max(TH_abund$F.D.RNA.P3.B4.H.01.ID09.06.15.THS.8.16S.4)
# The most frequent sps had 3309 reads

##### Octaves ##### 
# Using the sads package

oct_vect_sadsPack <- octav(TH_abund$F.D.RNA.P3.B4.H.01.ID09.06.15.THS.8.16S.4)

svg(paste(pathOutExploratory,"TH_oct_dist_sadsPack1.svg", sep = ""))
barplot(oct_vect_sadsPack$Freq, names.arg = oct_vect_sadsPack$upper)
dev.off()

# It gives a different result when I apply the function from the HPC week 

# Using the octaves function from the HPC week
octaves <- function(a){
  
  # a is the abundance vector. We are gonna transform the abundance vector doing an operation on it. This operation is gonna transform each value of abundance into the octave it belongs to. The abundance and the octaves are related by this formula. Floor takes the largest integer.
  octave_a_belongs <- floor(log2(a))+1
  
  # We need to create a empty vector where we are gona store the counts, or the frecuency, of the different abundances.creates a vector of 0. how many terms do we need the vector to have? we need to do some calculations. We need to work out what is the highst octave we are gonna get according to the abundances we've got.
  n_bins <- floor(log2(max(a)))+1
  
  # Then we create the vector of 0 with the right length
  octaves_bins = rep(0,n_bins)
  
  # The final vector we want to get is that one where the position in the vector indicates the bin. so the 3rd position of the vector is the the 3rd octave. The info we have is to which octaves the abundances belong to(that's octave_a_belongs). What we want now is to, kind of, work out the frecuencies = how many abundaces are in the 3rd octave = that number should go in the 3rd position of the vector we want ultimally to get (octaves). So we iterate over a_belongs_octave, each number give us the position in the octave vector we want to count, thats is "octaves[i]". then we want to add one every time because we are counting the number of abundances for every octave class or bin
  for(i in octave_a_belongs){
    octaves_bins[i] = octaves_bins[i] + 1
  }
  return(octaves_bins)
}
oct_vec1 <- octaves(TH_abund$F.D.RNA.P3.B4.H.01.ID09.06.15.THS.8.16S.4)
svg(paste(pathOutExploratory,"TH_Oct_dist1.svg", sep = ""))
barplot(oct_vec1, main = "Species abundance distribution", xlab = "Octaves", names.arg = c("1-2", "2-4", "4-8", "8-16", "16-32", "32-64", "64-128", "128-256", "512-1024", "1024-2048", "2048-4096", "4096-8192"), width = 10)
dev.off()
