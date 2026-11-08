# Set parameters for netlings1.1

# Author: Natasha Puzovic
# Date created: 07.12.2022
# Mail: npuzovic@protonmail.com
# Date last edited: 26.01.2023

# This script generates parameter files for multi-scenario simulation runs.
# It takes as input the output files from a sim that established the mutation-selection 
# equilibrium of noise and sbasal genotypes. These equilibrium values become the starting 
# point of the new, multi-scenario sim runs.

# for further sims
set.seed(1e4)
options(scipen = 999)
numReplicates = 10
sOptDisplacement = 20
numNets=500

options(warn=2) 

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
MUTATION_RATE_NET_PER_LINK = 0.0001
MUTATION_DISTR_NET_LOWERBOUND = -4
MUTATION_DISTR_NET_UPPERBOUND = 4
REC_RATE = 0.05
NUM_OUT_GENS = 500
PRNG_SEED = 1 #floor(runif(1, 1, 2e+8))

prmsS0 = as.data.frame(rbind(MUTABLE_NOISE, MUTABLE_S_BASAL, MUTABLE_NET, SEL, CONST_SEL, 
                            NUM_NODES, NUM_TIMESTEPS, S_MAX, POP_SIZE, NUM_GENERATIONS, 
                            NOISE_EXT, MUTATION_RATE_NOISE_PER_GENE, MUTATION_DISTR_NOISE_LOWERBOUND,
                            MUTATION_DISTR_NOISE_UPPERBOUND, MUTATION_RATE_SBASAL_PER_GENE,
                            MUTATION_DISTR_SBASAL_LOWERBOUND, MUTATION_DISTR_SBASAL_UPPERBOUND, 
                            MUTATION_RATE_NET_PER_LINK, MUTATION_DISTR_NET_LOWERBOUND, MUTATION_DISTR_NET_UPPERBOUND,
                            REC_RATE, NUM_OUT_GENS, PRNG_SEED))

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
  
  ##### directional sel. (39, 1)
  prmsV <- prmsV0
  prmsV["S_OPT_DELTA", differentlySelectedGene] <- sOptDisplacement
  write.table(prmsV, file = paste0("prmsV_39stab_1dir"), sep = " ", row.names = TRUE, col.names = FALSE, quote = FALSE)
  prmsS <- prmsS0
  for(replNum in 1:numReplicates)
  {
    prmsS[which(rownames(prmsS) == "PRNG_SEED"), ] = floor(runif(1, 1, 2e+8))
    write.table(prmsS, file = paste0("prmsS_sel_rep", replNum), sep = " ", row.names = TRUE, col.names = FALSE, quote = FALSE)
  }
  
  ##### neutrality
  prmsV <- prmsV0
  prmsV["SELPRESS", ] <- 0
  write.table(prmsV, file = paste0("prmsV_neu"), sep = " ", row.names = TRUE, col.names = FALSE, quote = FALSE)
  prmsS <- prmsS0
  prmsS[which(rownames(prmsS) == "SEL"), ] <- 0
  prmsS[which(rownames(prmsS) == "CONST_SEL"), ] <- 0
  
  for(replNum in 1:numReplicates)
  {
    prmsS[which(rownames(prmsS) == "PRNG_SEED"), ] = floor(runif(1, 1, 2e+8))
    write.table(prmsS, file = paste0("prmsS_neu_rep", replNum), sep = " ", row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  setwd("..")
}

