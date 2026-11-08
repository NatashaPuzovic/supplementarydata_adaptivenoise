# Author: Natasha Puzovic (puzovic@evolbio.mpg.de)
# Date created: 07 Jan 2023
# Last modified: 07 Feb 2023

# give folder name
#path = paste0("/home/npuzovic/expressionnoise/expressionnoise_master/results/", resultsFolder)
#setwd(path)

numReps = 10
adaptationResponseCriteria = 0.1 # the interval in which adaptation is successful

library(dplyr)
library(reshape2)

# read params
prmsS <- read.table(file = "prmsS_dirSelPos_rep1", sep = " ", header = FALSE, row.names = 1)
NUM_NODES = prmsS[which(rownames(prmsS) == "NUM_NODES"), ]
NUM_GENERATIONS = prmsS[which(rownames(prmsS) == "NUM_GENERATIONS"), ]
NUM_OUT_GENS = prmsS[which(rownames(prmsS) == "NUM_OUT_GENS"), ]
S_MAX = prmsS[which(rownames(prmsS) == "S_MAX"), ]
node_names = paste("G", 1:NUM_NODES, sep = "")
prmsV <- read.table("prmsV_39stab_1dirSelPos", sep = " ", header = FALSE, row.names = 1)
colnames(prmsV) <- node_names
# which genes had a change of selection regime?
indexSelGene <- as.integer(which(prmsV["S_OPT_DELTA", ] != 0))
sOptSelGene <- prmsV["S_OPT", indexSelGene] + prmsV["S_OPT_DELTA", indexSelGene]

# names
meanGenotEtaColnames <- paste("meanGEta_", node_names, sep = "")
meanGenotSbasalColnames <- paste("meanGSbasal_", node_names, sep = "")
meanPhenotColnames <- paste("meanP_", node_names, sep = "")
varPhenotColnames <- paste("varP_", node_names, sep = "")

allReplStats <- data.frame(matrix(nrow = 0, ncol = 17))
colnames(allReplStats) <- c("scen", "rep", "generation", "gene", 
                            "meanGenotEta", "meanGenotEtaFirstGen",
                            "meanGenotSbasal", "meanGenotSbasalFirstGen",
                            "meanP", "varP", 
                            "CVP", "noiseP", "FanoP",
                            "relDeltaVar", "relDeltaCV", "relDeltaNoise", "relDeltaFano")

# read directional selection results
for(replNum in 1:numReps) 
{
  # read sim results
  hist <- read.table(file = paste0("dirSelPos_rep", replNum, ".hist"), sep = " ", header = FALSE)
  hist <- hist[, -(dim(hist)[2])]
  colnames(hist) <- c("generation", meanPhenotColnames, varPhenotColnames, "meanF")
  PVF_long_meanP <- melt(hist[, c("generation", meanPhenotColnames)], 
                      id.vars = c("generation"),
                      variable.name = "gene",
                      value.name = "meanP")
  PVF_long_meanP$gene <- gsub("meanP_", "", PVF_long_meanP$gene)
  PVF_long_varP <- melt(hist[, c("generation", varPhenotColnames)],
                           id.vars = c("generation"), 
                           variable.name = "gene",
                           value.name = "varP")
  noises <- read.table(file = paste0("dirSelPos_rep", replNum, ".noises"), sep = " ", header = FALSE)
  noises <- noises[, -(dim(noises)[2])]
  colnames(noises) <- c("generation", meanGenotEtaColnames)
  noisesSumm <- melt(noises[, c("generation", meanGenotEtaColnames)],
                           id.vars = c("generation"), 
                           variable.name = "gene",
                           value.name = "meanGenotEta")
  
  sbasals <- read.table(file = paste0("dirSelPos_rep", replNum, ".sbasals"), sep = " ", header = FALSE)
  sbasals <- sbasals[, -(dim(sbasals)[2])]
  colnames(sbasals) <- c("generation", meanGenotSbasalColnames)
  sbasalsSumm <- melt(sbasals[, c("generation", meanGenotSbasalColnames)],
                     id.vars = c("generation"), 
                     variable.name = "gene",
                     value.name = "meanGenotSbasal")
  generations <- unique(hist$generation)
  
  # make df
  oneReplStats <- data.frame(matrix(nrow = NUM_NODES*(length(generations)), ncol = 17))
  colnames(oneReplStats) <- c("scen", "rep", "generation", "gene", 
                              "meanGenotEta", "meanGenotEtaFirstGen",
                              "meanGenotSbasal", "meanGenotSbasalFirstGen",
                              "meanP", "varP", 
                              "CVP", "noiseP", "FanoP",
                              "relDeltaVar", "relDeltaCV", "relDeltaNoise", "relDeltaFano")
  oneReplStats$scen <- "dirSelPos"
  oneReplStats$rep <- replNum
  oneReplStats$generation <- PVF_long_meanP$generation
  oneReplStats$gene <- PVF_long_meanP$gene
  oneReplStats$meanGenotEta <- noisesSumm$meanGenotEta
  oneReplStats$meanGenotSbasal <- sbasalsSumm$meanGenotSbasal
  
  for(geneNum in 1:NUM_NODES)
  {
    geneName = node_names[geneNum]
    oneReplStats[oneReplStats$gene == geneName, "meanGenotEtaFirstGen"] <-
      oneReplStats[oneReplStats$generation == 1, "meanGenotEta"][geneNum]
    oneReplStats[oneReplStats$gene == geneName, "meanGenotSbasalFirstGen"] <-
      oneReplStats[oneReplStats$generation == 1, "meanGenotSbasal"][geneNum]
  }
  oneReplStats$meanP <- PVF_long_meanP$meanP
  oneReplStats$varP <- PVF_long_varP$varP
  
  # calculate additional measures
  oneReplStats$CVP <- sqrt(oneReplStats$varP)/oneReplStats$meanP
  oneReplStats$noiseP <- oneReplStats$varP/((oneReplStats$meanP)^2)
  oneReplStats$FanoP <- oneReplStats$varP/oneReplStats$meanP
  
  oneReplStats$relDeltaVar <- oneReplStats$varP - rep(oneReplStats[oneReplStats$generation == 1, "varP"], each = length(generations))
  oneReplStats$relDeltaVar <- 
    oneReplStats$relDeltaVar/(oneReplStats$varP + rep(oneReplStats[oneReplStats$generation == 1, "varP"], each = length(generations)))
  
  oneReplStats$relDeltaCV <- oneReplStats$CVP - rep(oneReplStats[oneReplStats$generation == 1, "CVP"], each = length(generations))
  oneReplStats$relDeltaCV <- 
    oneReplStats$relDeltaCV/(oneReplStats$CVP + rep(oneReplStats[oneReplStats$generation == 1, "CVP"], each = length(generations)))
  
  oneReplStats$relDeltaNoise <- oneReplStats$noiseP - rep(oneReplStats[oneReplStats$generation == 1, "noiseP"], each = length(generations))
  oneReplStats$relDeltaNoise <- 
    oneReplStats$relDeltaNoise/(oneReplStats$noiseP + rep(oneReplStats[oneReplStats$generation == 1, "noiseP"], each = length(generations)))
  
  oneReplStats$relDeltaFano <- oneReplStats$FanoP - rep(oneReplStats[oneReplStats$generation == 1, "FanoP"], each = length(generations))
  oneReplStats$relDeltaFano <- 
    oneReplStats$relDeltaFano/(oneReplStats$FanoP + rep(oneReplStats[oneReplStats$generation == 1, "FanoP"], each = length(generations)))
  
  # add to big df
  allReplStats <- rbind(allReplStats, oneReplStats)
}

#### dplyr
allReplStats <- allReplStats %>%
  group_by(gene, generation) %>%
  summarise(meanGenotEta = mean(meanGenotEta),
            #meanGenotEtaFirstGen = mean(meanGenotEtaFirstGen), # mean func just for consistency
            meanGenotSbasal = mean(meanGenotSbasal),
            #meanGenotSbasalFirstGen = mean(meanGenotSbasalFirstGen), # mean func just for consistency
            meanP = mean(meanP),
            varP = mean(varP),
            CVP = mean(CVP),
            noiseP = mean(noiseP),
            FanoP = mean(FanoP),
            relDeltaVar = mean(relDeltaVar),
            relDeltaCV = mean(relDeltaCV),
            relDeltaNoise = mean(relDeltaNoise),
            relDeltaFano = mean(relDeltaFano)) 

# add average var, noise genotype, sbasal genotype from MSBalance
MSBHist <- read.table("mutSelBalance.hist", header = FALSE)
colnames(MSBHist) <- c("generation", meanPhenotColnames, varPhenotColnames, "meanF")
# read noises file
MSBNoises <- read.table("mutSelBalance.noises", header = FALSE)
colnames(MSBNoises) <- c("generation", meanGenotEtaColnames)
# read sbasals file
MSBSbasals <- read.table("mutSelBalance.sbasals", header = FALSE)
colnames(MSBSbasals) <- c("generation", meanGenotSbasalColnames)
prmsS_forMutSelBalance <- read.table(file = "prmsS_forMutSelBalance", sep = " ", header = FALSE, row.names = 1)
NUM_GENERATIONS_forMutSelBalance = prmsS_forMutSelBalance[which(rownames(prmsS_forMutSelBalance) == "NUM_GENERATIONS"), ]

aveVarP <- colSums(MSBHist[MSBHist$generation > NUM_GENERATIONS_forMutSelBalance/2, varPhenotColnames])/
  length(which(MSBHist$generation > NUM_GENERATIONS_forMutSelBalance/2))
aveMeanGEta <- colSums(MSBNoises[MSBNoises$generation > NUM_GENERATIONS_forMutSelBalance/2, meanGenotEtaColnames])/
  length(which(MSBNoises$generation > NUM_GENERATIONS_forMutSelBalance/2))
aveMeanGSBasal <- colSums(MSBSbasals[MSBSbasals$generation > NUM_GENERATIONS_forMutSelBalance/2, meanGenotSbasalColnames])/
  length(which(MSBSbasals$generation > NUM_GENERATIONS_forMutSelBalance/2))
MSBAve <- data.frame(gene = node_names, MSBAveVarP = aveVarP, MSBAveMeanGEta = aveMeanGEta, MSBAveMeanGSBasal = aveMeanGSBasal)

aveMeanPheno <- colSums(MSBHist[MSBHist$generation > NUM_GENERATIONS_forMutSelBalance/2, meanPhenotColnames])/
  length(which(MSBHist$generation > NUM_GENERATIONS_forMutSelBalance/2))
MSBAveForOut <- data.frame(gene = node_names, 
                           S_OPT = unlist(prmsV[rownames(prmsV) == "S_OPT", ]), 
                           MSBAveMeanPheno = aveMeanPheno,
                           MSBAveVarP = aveVarP,
                           MSBAveMeanGEta = aveMeanGEta,
                           MSBAveMeanGSBasal = aveMeanGSBasal)
# write summary of this network results
write.table(MSBAveForOut,
            file = "MSBAve.txt",
            sep = "\t", row.names = FALSE)

### Directional selection
# find time to adaptation
selGeneStats <- allReplStats[allReplStats$gene == node_names[indexSelGene], ]
adaptedGenerationsIndices <- which(selGeneStats$meanP > (sOptSelGene - S_MAX*adaptationResponseCriteria) &
                                     selGeneStats$meanP < (sOptSelGene + S_MAX*adaptationResponseCriteria))
# if the optimal range is reached at any point
if(length(adaptedGenerationsIndices) != 0)
{
  # if the mean expression level stays in the optimal range after the first instance of adaptation
  if(identical(adaptedGenerationsIndices, adaptedGenerationsIndices[1]:adaptedGenerationsIndices[length(adaptedGenerationsIndices)]))
  {
    timeToAdapt <- selGeneStats$generation[adaptedGenerationsIndices[1]]
  } else
  {
    #timeToAdapt <- "NC_adaptation"
    timeToAdapt <- selGeneStats$generation[adaptedGenerationsIndices[1]] # it's not strictly convergent, but convergent nonetheless
  }
} else # if adaptation never happens
{
  timeToAdapt <- "Nonadaptive"
}

# get selpress
selPress <- allReplStats
selPress <- inner_join(selPress, MSBAve, by = "gene")

selPress$changeVar <- selPress$varP - selPress$MSBAveVarP
selPress$changeGenotEta <- selPress$meanGenotEta - selPress$MSBAveMeanGEta
selPress$changeGenotSbasal <- selPress$meanGenotSbasal - selPress$MSBAveMeanGSBasal

selPress$normChangeVar <- selPress$changeVar/selPress$MSBAveVarP
selPress$normChangeGenotEta <- selPress$changeGenotEta/selPress$MSBAveMeanGEta
selPress$normChangeGenotSbasal <- selPress$changeGenotSbasal/selPress$MSBAveMeanGSBasal

# keep only first and last gen, and cast to long format 
firstLastgen <- allReplStats[allReplStats$generation %in% c(1, NUM_GENERATIONS), ]
firstLastgen <- melt(firstLastgen, id.vars = c("gene", "generation"))
firstLastgen <- dcast(firstLastgen, id.vars = c("gene"), gene ~ variable + generation)

dirSelPos <- firstLastgen[match(node_names, firstLastgen$gene), ] # reorder gene names

### directional selection
# if adaptation occurs
if(timeToAdapt != "Nonadaptive")
{
  # adaptation phase
  selPress_dir_A <- selPress[selPress$generation <= timeToAdapt, ] # include the gen where the new sopt was reached
  numTotalOutGenerations <- dim(hist)[1]
  numOutGenerationsDuringAdaptation <- numTotalOutGenerations - length(adaptedGenerationsIndices) + 1
  selPress_dir_A <- selPress_dir_A %>%
    group_by(gene) %>%
    summarise(aveNormChangeVar_A = sum(normChangeVar)/numOutGenerationsDuringAdaptation,
              selpress_eta_A = sum(normChangeGenotEta)/numOutGenerationsDuringAdaptation,
              selpress_sbasal_A = sum(normChangeGenotSbasal)/numOutGenerationsDuringAdaptation)
  
  selPress_dir_A$selpress_eta_A_abs = abs(selPress_dir_A$selpress_eta_A)
  selPress_dir_A$selpress_sbasal_A_abs = abs(selPress_dir_A$selpress_sbasal_A)
  selPress_dir_A$gene <- factor(selPress_dir_A$gene, levels = node_names)

  # postadaptation phase
  numOutGenerationsAfterAdaptation <- numTotalOutGenerations - numOutGenerationsDuringAdaptation
  selPress_dir_PA <- selPress[selPress$generation > timeToAdapt, ]
  selPress_dir_PA <- selPress_dir_PA %>%
    group_by(gene) %>%
    summarise(aveNormChangeVar_PA = sum(normChangeVar)/numOutGenerationsAfterAdaptation,
              selpress_eta_PA = sum(normChangeGenotEta)/numOutGenerationsAfterAdaptation,
              selpress_sbasal_PA = sum(normChangeGenotSbasal)/numOutGenerationsAfterAdaptation)
  
  selPress_dir_PA$selpress_eta_PA_abs = abs(selPress_dir_PA$selpress_eta_PA)
  selPress_dir_PA$selpress_sbasal_PA_abs = abs(selPress_dir_PA$selpress_sbasal_PA)
  selPress_dir_PA$gene <- factor(selPress_dir_PA$gene, levels = node_names)
}

# if adaptation never occurs
if(timeToAdapt == "Nonadaptive")
{
  dirSelPos$aveNormChangeVar_A = NA
  dirSelPos$selpress_eta_A = NA
  dirSelPos$selpress_sbasal_A = NA
  dirSelPos$selpress_eta_A_abs = NA
  dirSelPos$selpress_sbasal_A_abs = NA
  
  dirSelPos$aveNormChangeVar_PA = NA
  dirSelPos$selpress_eta_PA = NA
  dirSelPos$selpress_sbasal_PA = NA
  dirSelPos$selpress_eta_PA_abs = NA
  dirSelPos$selpress_sbasal_PA_abs = NA
  
  dirSelPos$timeToAdapt_dirSel <- NA
  dirSelPos$timeToAdapt_dirSel[dirSelPos$gene == node_names[indexSelGene]] <- timeToAdapt
} else 
{
  dirSelPos <- inner_join(dirSelPos, selPress_dir_A, by = "gene")
  dirSelPos <- inner_join(dirSelPos, selPress_dir_PA, by = "gene")
}

dirSelPos$timeToAdapt_dirSel <- "Not_applied"
dirSelPos$timeToAdapt_dirSel[dirSelPos$gene == node_names[indexSelGene]] <- timeToAdapt
dirSelPos$sOptDelta <- unlist(prmsV[which(rownames(prmsV) == "S_OPT_DELTA"), ])
dirSelPos <- cbind(rep("dirSelPos", NUM_NODES), dirSelPos)
colnames(dirSelPos)[1] <- "scen"

#if(file.exists("allReplEvolResults_dirSelPos.txt"))
#{
#  stop("Summarized data file already present.")
#}
# write summary of this network results
write.table(dirSelPos,
            file = "allReplEvolResults_dirSelPos.txt",
            sep = "\t", row.names = FALSE)
