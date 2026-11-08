#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# Set parameters for netlings1.1

# Author: Natasha Puzovic
# Date created: 07.12.2022
# Mail: npuzovic@protonmail.com
# Date last edited: 26.01.2023

# This script generates parameter files for multi-scenario simulation runs.
# It takes as input the output files from a sim that established the mutation-selection 
# equilibrium of noise and sbasal genotypes. These equilibrium values become the starting 
# point of the new, multi-scenario sim runs.

if (length(args) != 2)
{
  stop("Please specify results folder; number of networks.", call. = FALSE)
} else if (length(args) == 2)
{
  baseFolder = args[1]
  numNets = args[2]
}

setwd(paste0("sim_data/", baseFolder))

Rseed = 454
set.seed(Rseed)
options(scipen = 999)
gen = 1
numNodes = 40
popSize = 1000

numReplicates = 10

#options(warn=2) 

# make prmsScalar
MUTABLE_NOISE = 1
MUTABLE_S_BASAL = 1
MUTABLE_NET = 0
SEL = 1
CONST_SEL = 1
NUM_NODES = 40
NUM_TIMESTEPS = 50
S_MAX = 100
POP_SIZE = 1000
NUM_GENERATIONS = 10000
NOISE_EXT = 0
MUTATION_RATE_NOISE_PER_GENE = 0.01
MUTATION_DISTR_NOISE_LOWERBOUND = 0
MUTATION_DISTR_NOISE_UPPERBOUND = 200
MUTATION_RATE_SBASAL_PER_GENE = 0.01
MUTATION_DISTR_SBASAL_LOWERBOUND = 0
MUTATION_DISTR_SBASAL_UPPERBOUND = 100
MUTATION_RATE_NET_PER_LINK = 0.001
MUTATION_DISTR_NET_LOWERBOUND = -4
MUTATION_DISTR_NET_UPPERBOUND = 4
REC_RATE = 0.05
NUM_OUT_GENS = 500

# dirSelPos
# go into each net folder
for(netNum in 1:numNets)
{
  setwd(paste0("net_", netNum))
  
  # replicates
  for(replNum in 1:numReplicates)
  {
    PRNG_SEED = floor(runif(1, 1, 2e+8))
    prmsS = as.data.frame(rbind(MUTABLE_NOISE, MUTABLE_S_BASAL, MUTABLE_NET, SEL, CONST_SEL, 
                                NUM_NODES, NUM_TIMESTEPS, S_MAX, POP_SIZE, NUM_GENERATIONS, 
                                NOISE_EXT, MUTATION_RATE_NOISE_PER_GENE, MUTATION_DISTR_NOISE_LOWERBOUND,
                                MUTATION_DISTR_NOISE_UPPERBOUND, MUTATION_RATE_SBASAL_PER_GENE,
                                MUTATION_DISTR_SBASAL_LOWERBOUND, MUTATION_DISTR_SBASAL_UPPERBOUND, 
                                MUTATION_RATE_NET_PER_LINK, MUTATION_DISTR_NET_LOWERBOUND, MUTATION_DISTR_NET_UPPERBOUND,
                                REC_RATE, NUM_OUT_GENS, PRNG_SEED))
    colnames(prmsS) = NULL
    write.table(prmsS, file = paste0("prmsS_dirSelPos_rep", replNum), sep = " ", row.names = TRUE, col.names = FALSE, quote = FALSE)
  }
  print(paste0("Written scalar parameters for net ", netNum, "."))
  setwd("..")
}

# dirSelNeg
# go into each net folder
for(netNum in 1:numNets)
{
  setwd(paste0("net_", netNum))
  
  # replicates
  for(replNum in 1:numReplicates)
  {
    PRNG_SEED = floor(runif(1, 1, 2e+8))
    prmsS = as.data.frame(rbind(MUTABLE_NOISE, MUTABLE_S_BASAL, MUTABLE_NET, SEL, CONST_SEL, 
                                NUM_NODES, NUM_TIMESTEPS, S_MAX, POP_SIZE, NUM_GENERATIONS, 
                                NOISE_EXT, MUTATION_RATE_NOISE_PER_GENE, MUTATION_DISTR_NOISE_LOWERBOUND,
                                MUTATION_DISTR_NOISE_UPPERBOUND, MUTATION_RATE_SBASAL_PER_GENE,
                                MUTATION_DISTR_SBASAL_LOWERBOUND, MUTATION_DISTR_SBASAL_UPPERBOUND, 
                                MUTATION_RATE_NET_PER_LINK, MUTATION_DISTR_NET_LOWERBOUND, MUTATION_DISTR_NET_UPPERBOUND,
                                REC_RATE, NUM_OUT_GENS, PRNG_SEED))
    colnames(prmsS) = NULL
    write.table(prmsS, file = paste0("prmsS_dirSelNeg_rep", replNum), sep = " ", row.names = TRUE, col.names = FALSE, quote = FALSE)
  }
  print(paste0("Written scalar parameters for net ", netNum, "."))
  setwd("..")
}

# fluctSel
CONST_SEL = 0
# go into each net folder
for(netNum in 1:numNets)
{
  setwd(paste0("net_", netNum))
  
  # replicates
  for(replNum in 1:numReplicates)
  {
    PRNG_SEED = floor(runif(1, 1, 2e+8))
    prmsS = as.data.frame(rbind(MUTABLE_NOISE, MUTABLE_S_BASAL, MUTABLE_NET, SEL, CONST_SEL, 
                                NUM_NODES, NUM_TIMESTEPS, S_MAX, POP_SIZE, NUM_GENERATIONS, 
                                NOISE_EXT, MUTATION_RATE_NOISE_PER_GENE, MUTATION_DISTR_NOISE_LOWERBOUND,
                                MUTATION_DISTR_NOISE_UPPERBOUND, MUTATION_RATE_SBASAL_PER_GENE,
                                MUTATION_DISTR_SBASAL_LOWERBOUND, MUTATION_DISTR_SBASAL_UPPERBOUND, 
                                MUTATION_RATE_NET_PER_LINK, MUTATION_DISTR_NET_LOWERBOUND, MUTATION_DISTR_NET_UPPERBOUND,
                                REC_RATE, NUM_OUT_GENS, PRNG_SEED))
    colnames(prmsS) = NULL
    write.table(prmsS, file = paste0("prmsS_fluctSel_rep", replNum), sep = " ", row.names = TRUE, col.names = FALSE, quote = FALSE)
  }
  print(paste0("Written scalar parameters for net ", netNum, "."))
  setwd("..")
}