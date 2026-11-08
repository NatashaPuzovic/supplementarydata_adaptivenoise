#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# Create parameter files for netlings1.1

# Author: Natasha Puzovic
# Date created: 12.02.2022
# Date modified: 12.02.2023
# Mail: npuzovic@protonmail.com

if (length(args) != 2)
{
  stop("Please specify results folder; number of networks.", call. = FALSE)
} else if (length(args) == 2)
{
  baseFolder = args[1]
  numNets = args[2]
}

setwd(paste0("sim_data/", baseFolder))

# This script generates two parameter files in the format of a table.
Rseed = 450
set.seed(Rseed)
options(scipen = 999)
gen = 1
numNodes = 40
popSize = 1000

# make prmsV (only one is needed for the entire batch)
SELPRESS = rep(1, numNodes)
S_OPT = rep(50, numNodes)
S_OPT_DELTA = rep(0, numNodes)

prmsV = as.data.frame(rbind(SELPRESS, S_OPT, S_OPT_DELTA))
colnames(prmsV) = NULL
write.table(prmsV, file = paste0("scripts_and_params/prmsV_forNetEst"),
            sep = " ", row.names = TRUE, col.names = FALSE, quote = FALSE)
print("Written vector parameters to prmsV_forNetEst.")

# make prmss (only one is needed per network topo)
MUTABLE_NOISE = 0
MUTABLE_S_BASAL = 0
MUTABLE_NET = 1
SEL = 1
CONST_SEL = 1
NUM_NODES = 40
NUM_TIMESTEPS = 50
S_MAX = 100
POP_SIZE = 1000
NUM_GENERATIONS = 5000
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

for(netNum in 1:numNets)
{
  setwd(paste0("net_", netNum))
  PRNG_SEED = floor(runif(1, 1, 2e+8))
  prmsS = as.data.frame(rbind(MUTABLE_NOISE, MUTABLE_S_BASAL, MUTABLE_NET, SEL, CONST_SEL, 
                              NUM_NODES, NUM_TIMESTEPS, S_MAX, POP_SIZE, NUM_GENERATIONS, 
                              NOISE_EXT, MUTATION_RATE_NOISE_PER_GENE, MUTATION_DISTR_NOISE_LOWERBOUND,
                              MUTATION_DISTR_NOISE_UPPERBOUND, MUTATION_RATE_SBASAL_PER_GENE,
                              MUTATION_DISTR_SBASAL_LOWERBOUND, MUTATION_DISTR_SBASAL_UPPERBOUND, 
                              MUTATION_RATE_NET_PER_LINK, MUTATION_DISTR_NET_LOWERBOUND, MUTATION_DISTR_NET_UPPERBOUND,
                              REC_RATE, NUM_OUT_GENS, PRNG_SEED))
  colnames(prmsS) = NULL
  write.table(prmsS, file = paste0("prmsS_forNetEst"),
              sep = " ", row.names = TRUE, col.names = FALSE, quote = FALSE)
  print(paste0("Written scalar parameters to prmsS_forNetEst for net ", netNum, "."))
  setwd("..")
}







