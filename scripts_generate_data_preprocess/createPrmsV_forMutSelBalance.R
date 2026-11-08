#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# Create population file for netlings1.1

# Author: Natasha Puzovic
# Date created: 07.12.2022
# Date modified: 06.01.2023
# Mail: npuzovic@protonmail.com

# This script generates a parameter file in the format of a table.
# It takes as input a realized population file and creates a vector parameter file.

if (length(args) != 2)
{
  stop("Please specify results folder; network topology type, e.g. ER, BA, WS; density; number of networks.", call. = FALSE)
} else if (length(args) == 2)
{
  baseFolder = args[1]
  numNets = args[2]
}

options(scipen = 999)
gen = 1
numNodes = 40
popSize = 1000

for(netNum in 1:numNets)
{
  setwd(paste0("net_", netNum))
  netEstPop <- as.matrix(read.table("netEst.pop", sep = " ", header = FALSE))
  
  # choose fittest net
  fittestNet <- netEstPop[which(netEstPop[, ncol(netEstPop)] == max(netEstPop[, ncol(netEstPop)]))[1], ]
  
  # calculate indices of genotypes, phenotypes, and fitness
  numIDCols = 2
  indIDCols = c(1, 2)
  #indGenotNoiseCols = seq(from = numIDCols  + 1, to = numIDCols  + numNodes, by = 1)
  #indGenotSBasalCols = seq(from = numIDCols  + numNodes + 1, to = numIDCols  + 2*numNodes, by = 1)
  #indGenotNetCols = seq(from = numIDCols  + 2*numNodes + 1, to = numIDCols  + 2*numNodes + numNodes^2, by = 1)
  indPhenotypeCols = seq(from = numIDCols  + 2*numNodes + numNodes^2 + 1, to = numIDCols  + 3*numNodes + numNodes^2, by = 1)
  #indFitnessCol = numIDCols + length(indGenotNoiseCols) + length(indGenotSBasalCols) + 
                  #length(indGenotNetCols) + length(indPhenotypeCols) + 1
  sSteady = fittestNet[indPhenotypeCols]
  
  # make prmsVector
  SELPRESS = rep(1, numNodes)
  S_OPT = sSteady
  S_OPT_DELTA = rep(0, numNodes)
  
  prmsV = as.data.frame(rbind(SELPRESS, S_OPT, S_OPT_DELTA))
  colnames(prmsV) = NULL
  write.table(prmsV, file = paste0("prmsV_forMutSelBalance"),
              sep = " ", row.names = TRUE, col.names = FALSE, quote = FALSE)
  print(paste0("Written vector parameters to prmsV_forMutSelBalance for net ", netNum, "."))
  setwd("..")
}