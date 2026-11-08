#!/usr/bin/env Rscript

# Set parameters for netlings1.1

# Author: Natasha Puzovic
# Date created: 07.12.2022
# Mail: npuzovic@protonmail.com
# Date last edited: 14.02.2023

# This script generates parameter files for multi-scenario simulation runs.
# It takes as input the output files from a sim that established the mutation-selection 
# equilibrium of noise and sbasal genotypes. These equilibrium values become the starting 
# point of the new, multi-scenario sim runs.

# for further sims
options(scipen = 999)
numReplicates = 10
sOptDisplacement = 20
numNets = 2000

options(warn=2) 


# go into each net folder
for(netNum in 1:numNets)
{
  setwd(paste0("net_", netNum))  
  print(paste0("working on net ", netNum)) 
  
  # read params from mutSelBalance establishment
  prmsV0 <- read.table(file = "prmsV_forMutSelBalance", sep = " ", row.names = 1)

  # choose 1 gene for dir. and fluct. sel
  differentlySelectedGene = sample(1:NUM_NODES, 1)
  while((prmsV0["S_OPT", differentlySelectedGene] + sOptDisplacement > S_MAX) ||
        (prmsV0["S_OPT", differentlySelectedGene] - sOptDisplacement < 0))
  {
    print("New S_OPT is below or above expression level boundaries. Picking another gene at random to be selected.")
    differentlySelectedGene = sample(1:NUM_NODES, 1)
  }
  
  ##### directional sel pos
  prmsV <- prmsV0
  prmsV["S_OPT_DELTA", differentlySelectedGene] <- sOptDisplacement
  write.table(prmsV, file = paste0("prmsV_39stab_1dirSelPos"), sep = " ", row.names = TRUE, col.names = FALSE, quote = FALSE)

  ##### directional sel neg
  prmsV <- prmsV0
  prmsV["S_OPT_DELTA", differentlySelectedGene] <- -sOptDisplacement
  write.table(prmsV, file = paste0("prmsV_39stab_1dirSelNeg"), sep = " ", row.names = TRUE, col.names = FALSE, quote = FALSE)
  
  setwd("..")
}

