#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# Create population file for netlings1.1

# Author: Natasha Puzovic
# Date created: 07.12.2022
# Date modified: 06.01.2023
# Mail: npuzovic@protonmail.com

# This script generates a population file in the format of a table.
# It takes as input a network topology and creates a clonal population with
# hardcoded values for noise and sbasal genotypes (0 and 20, respectively).
# The phenotypes and fitness are set to 0.

if (length(args) != 2)
{
  stop("Please specify results folder; number of networks.", call. = FALSE)
} else if (length(args) == 2)
{
  baseFolder = args[1]
  numNets = args[2]
}

setwd(paste0("sim_data/", baseFolder))

Rseed = 451
set.seed(Rseed)
options(scipen = 999)
gen = 1
numNodes = 40
popSize = 1000

startingNoise = 0

for(netNum in 1:numNets)
{
  setwd(paste0("net_", netNum))
  
  rmtrx <- as.matrix(read.table("rmtrxFoundingNet", sep = " ", header = FALSE))
  if(dim(rmtrx)[2] != numNodes)
  {
    rmtrx <- rmtrx[, -dim(rmtrx)[2]]
  }
  
  startingRmtrx = c(t(rmtrx))
  genotNoises <- rep(startingNoise, numNodes)
  genotSbasals <- floor(runif(n = numNodes, min = 0, max = 100.99))
  genotRmtrx <- startingRmtrx
  phenot <- rep(0, numNodes)
  fitness <- 0
  
  oneIndividual <- c(gen, 1, genotNoises, genotSbasals, genotRmtrx, phenot, fitness)
  firstGenPop <- as.data.frame(lapply(oneIndividual, rep, popSize), col.names = NULL)
  firstGenPop[, 2] <- 1:popSize

  write.table(firstGenPop, "pop.netEst.firstgen", sep = " ", row.names = F, col.names = F)
  print(paste0("Written population file to pop.netEst.firstgen for net ", netNum, "."))
  setwd("..")
}


