# Create population file for netlings1.1

# Author: Natasha Puzovic
# Date created: 07.12.2022
# Date modified: 06.01.2023
# Mail: npuzovic[at]protonmail.com

# This script generates a population file in the format of a table.
# It takes as input a network topology and creates a clonal population with
# hardcoded values for noise and sbasal genotypes (0 and 20, respectively).
# The phenotypes and fitness are set to 0.

gen = 1
numNodes = 40
popSize = 1000

for(netNum in 1:numNets)
{
  setwd(paste0("net_", netNum))
  netEstPop <- as.matrix(read.table("netEst.pop", sep = " ", header = FALSE))
  
  # choose fittest net
  fittestIndividual <- netEstPop[which(netEstPop[, ncol(netEstPop)] == max(netEstPop[, ncol(netEstPop)]))[1], ]
  fittestIndividual[1] <- gen
  MSBPop <- as.data.frame(lapply(fittestIndividual, rep, popSize), col.names = NULL)
  MSBPop[, 2] <- 1:popSize
  
  write.table(MSBPop, "pop.MSB.firstgen", sep = " ", row.names = F, col.names = F)
  print(paste0("Written population file to pop.MSB.firstgen in net ", netNum, "."))
  
  setwd("..")
}

