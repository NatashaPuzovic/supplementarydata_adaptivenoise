# Author: Natasha Puzovic (puzovic@evolbio.mpg.de)
# Date created: 15 Jan 2023
# Last modified: 04 Feb 2023

# give folder name
#path = paste0("/home/npuzovic/expressionnoise/expressionnoise_master/results/", resultsFolder)
#setwd(path)

#if(file.exists("allReplEvolResults.txt")){
#  stop("Summarized data file already present.")
#}

numReps = 10
#adaptationResponseCriteria = 0.1 # the interval in which adaptation is successful

library(dplyr)
library(reshape2)

# read params
prmsS <- read.table(file = "prmsS_fluctSel_rep1", sep = " ", header = FALSE, row.names = 1)
NUM_NODES = prmsS[which(rownames(prmsS) == "NUM_NODES"), ]
NUM_GENERATIONS = prmsS[which(rownames(prmsS) == "NUM_GENERATIONS"), ]
NUM_OUT_GENS = prmsS[which(rownames(prmsS) == "NUM_OUT_GENS"), ]
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

# read fluctuating selection results
for(replNum in 1:numReps) 
{
  # read sim results
  hist <- read.table(file = paste0("fluctSel_rep", replNum, ".hist"), sep = " ", header = FALSE)
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
  noises <- read.table(file = paste0("fluctSel_rep", replNum, ".noises"), sep = " ", header = FALSE)
  noises <- noises[, -(dim(noises)[2])]
  colnames(noises) <- c("generation", meanGenotEtaColnames)
  noisesSumm <- melt(noises[, c("generation", meanGenotEtaColnames)],
                     id.vars = c("generation"), 
                     variable.name = "gene",
                     value.name = "meanGenotEta")
  
  sbasals <- read.table(file = paste0("fluctSel_rep", replNum, ".sbasals"), sep = " ", header = FALSE)
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
  oneReplStats$scen <- "fluctSel"
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

# get selpress
selPress <- allReplStats
selPress <- inner_join(selPress, MSBAve, by = "gene")

selPress$changeVar <- selPress$varP - selPress$MSBAveVarP
selPress$changeGenotEta <- selPress$meanGenotEta - selPress$MSBAveMeanGEta
selPress$changeGenotSbasal <- selPress$meanGenotSbasal - selPress$MSBAveMeanGSBasal

selPress$normChangeVar <- selPress$changeVar/selPress$MSBAveVarP
selPress$normChangeGenotEta <- selPress$changeGenotEta/selPress$MSBAveMeanGEta
selPress$normChangeGenotSbasal <- selPress$changeGenotSbasal/selPress$MSBAveMeanGSBasal

### fluctuating selection
selPress_fluct_A <- selPress
numTotalOutGenerations <- dim(hist)[1]
selPress_fluct_A <- selPress_fluct_A %>%
  group_by(gene) %>%
  summarise(aveNormChangeVar_A = sum(normChangeVar)/numTotalOutGenerations,
            selpress_eta_A = sum(normChangeGenotEta)/numTotalOutGenerations,
            selpress_sbasal_A = sum(normChangeGenotSbasal)/numTotalOutGenerations)

selPress_fluct_A$selpress_eta_A_abs = abs(selPress_fluct_A$selpress_eta_A)
selPress_fluct_A$selpress_sbasal_A_abs = abs(selPress_fluct_A$selpress_sbasal_A)
selPress_fluct_A$gene <- factor(selPress_fluct_A$gene, levels = node_names)

# keep only first and last gen, and cast to long format 
firstLastgen <- allReplStats[allReplStats$generation %in% c(1, NUM_GENERATIONS), ]
firstLastgen <- melt(firstLastgen, id.vars = c("gene", "generation"))
firstLastgen <- dcast(firstLastgen, id.vars = "gene", gene ~ variable + generation)

fluctSel <- firstLastgen[match(node_names, firstLastgen$gene), ] # reorder gene names

fluctSel <- inner_join(fluctSel, selPress_fluct_A, by = "gene")
fluctSel$sOptDelta <- unlist(prmsV[which(rownames(prmsV) == "S_OPT_DELTA"), ])
fluctSel <- cbind(rep("fluctSel", NUM_NODES), fluctSel)
colnames(fluctSel)[1] <- "scen"

# write summary of this network results
write.table(fluctSel,
            file = "allReplEvolResults_fluctSel.txt",
            sep = "\t", row.names = FALSE)
