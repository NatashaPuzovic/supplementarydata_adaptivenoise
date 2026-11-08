# give folder name
#path = paste0("/home/npuzovic/expressionnoise/expressionnoise_master/results/", resultsFolder)
#setwd(path)

if(file.exists("allReplEvolResults.txt")){
  stop("Summarized data file already present.")
}

# test if all the arguments are present: if not, return an error
if (length(args) != 1) {
  stop("Please specify results folder; network topology type, e.g. ER, BA, WS; density; number of networks.", call. = FALSE)
} else if (length(args) == 1) {
  numReps = args[1]
}


# read params
params <- read.table(file = "paramsEvolSel_rep1.txt", 
                     sep = " ", header = FALSE, col.names = c("param", "val"))
num_nodes = params$val[params$param == "NUM_NODES"]
num_generations = params$val[params$param == "NUM_GENERATIONS"]
starting_noise = params$val[params$param == "STARTING_NOISE_INT_ALL"]
num_out_gens = params$val[params$param == "NUM_OUT_GENS"]

# names
node_names = paste("G", 1:num_nodes, sep = "")
meanGenotColnames <- paste("meanG_", node_names, sep = "")
meanPhenotColnames <- paste("meanP_", node_names, sep = "")
varPhenotColnames <- paste("varP_", node_names, sep = "")

allReplStats <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(allReplStats) <- c("scen", "rep", "generation", "gene", 
                            "meanG", "meanP", "varP", 
                            "CVP", "noiseP", "FanoP",
                            "relDeltaVar", "relDeltaCV", "relDeltaNoise", "relDeltaFano",
                            "meanF") 

library(reshape2)
for(replNum in 1:numReps) {
  
  # read sim results
  GPFsumm <- read.table(file = paste("GPFsel_rep", replNum, ".summary", sep = ""),
                       sep = " ", header = FALSE)
  colnames(GPFsumm) <- c("generation", meanGenotColnames, meanPhenotColnames, varPhenotColnames, "meanF")
  generations <- unique(GPFsumm$generation)
  GPFsummLong <- melt(GPFsumm[, c("generation", meanGenotColnames)], 
                       id.vars = c("generation"), 
                       variable.name = "gene",
                       value.name = "meanG")
  GPFsummLong$gene <- gsub("meanG_", "", GPFsummLong$gene)
  GPFsummLong_meanP <- melt(GPFsumm[, c("generation", meanPhenotColnames)],
                       id.vars = c("generation"), 
                       variable.name = "gene",
                       value.name = "meanP")
  GPFsummLong_varP <- melt(GPFsumm[, c("generation", varPhenotColnames)],
                             id.vars = c("generation"), 
                             variable.name = "gene",
                             value.name = "varP")

  # make df
  oneReplStats <- data.frame(matrix(nrow = num_nodes*(length(generations)), ncol = 14))
  colnames(oneReplStats) <- c("scen", "rep", "generation", "gene", 
                              "meanG", "meanP", "varP", 
                              "CVP", "noiseP", "FanoP",
                              "relDeltaVar", "relDeltaCV", "relDeltaNoise", "relDeltaFano")
  oneReplStats$scen <- "sel"
  oneReplStats$rep <- replNum
  oneReplStats$generation <- GPFsummLong$generation
  oneReplStats$gene <- GPFsummLong$gene
  oneReplStats$meanG <- GPFsummLong$meanG
  oneReplStats$meanP <- GPFsummLong_meanP$meanP
  oneReplStats$varP <- GPFsummLong_varP$varP
  
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
  
  #oneReplStats$meanF <- 
  
  # add to big df
  allReplStats <- rbind(allReplStats, oneReplStats)
}

for(replNum in 1:numReps) {
  
  # read sim results
  GPFsumm <- read.table(file = paste("GPFneutr_rep", replNum, ".summary", sep = ""),
                        sep = " ", header = FALSE)
  colnames(GPFsumm) <- c("generation", meanGenotColnames, meanPhenotColnames, varPhenotColnames)
  generations <- unique(GPFsumm$generation)
  GPFsummLong <- melt(GPFsumm[, c("generation", meanGenotColnames)], 
                      id.vars = c("generation"), 
                      variable.name = "gene",
                      value.name = "meanG")
  GPFsummLong$gene <- gsub("meanG_", "", GPFsummLong$gene)
  GPFsummLong_meanP <- melt(GPFsumm[, c("generation", meanPhenotColnames)],
                            id.vars = c("generation"), 
                            variable.name = "gene",
                            value.name = "meanP")
  GPFsummLong_varP <- melt(GPFsumm[, c("generation", varPhenotColnames)],
                           id.vars = c("generation"), 
                           variable.name = "gene",
                           value.name = "varP")

  # make df
  oneReplStats <- data.frame(matrix(nrow = num_nodes*(length(generations)), ncol = 14))
  colnames(oneReplStats) <- c("scen", "rep", "generation", "gene", 
                              "meanG", "meanP", "varP", 
                              "CVP", "noiseP", "FanoP",
                              "relDeltaVar", "relDeltaCV", "relDeltaNoise", "relDeltaFano")
  oneReplStats$scen <- "neu"
  oneReplStats$rep <- replNum
  oneReplStats$generation <- GPFsummLong$generation
  oneReplStats$gene <- GPFsummLong$gene
  oneReplStats$meanG <- GPFsummLong$meanG
  oneReplStats$meanP <- GPFsummLong_meanP$meanP
  oneReplStats$varP <- GPFsummLong_varP$varP
  
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
library(dplyr)

allReplStats <- allReplStats %>%
  group_by(scen, gene, generation) %>%
  summarise(meanG = mean(meanG),
            meanP = mean(meanP),
            varP = mean(varP),
            CVP = mean(CVP),
            noiseP = mean(noiseP),
            FanoP = mean(FanoP),
            relDeltaVar = mean(relDeltaVar),
            relDeltaCV = mean(relDeltaCV),
            relDeltaNoise = mean(relDeltaNoise),
            relDeltaFano = mean(relDeltaFano))


g_static_area <- num_out_gens*starting_noise
p_static_area <- num_out_gens

selPress <- allReplStats %>%
  group_by(scen, gene) %>%
  summarise(s_g_area = sum(meanG - starting_noise)/g_static_area,
            s_g_area_abs = abs(s_g_area),
            s_p_area_relDeltaVar = sum(relDeltaVar)/p_static_area,
            s_p_area_relDeltaCV = sum(relDeltaCV)/p_static_area,
            s_p_area_relDeltaNoise = sum(relDeltaNoise)/p_static_area,
            s_p_area_relDeltaFano = sum(relDeltaFano)/p_static_area)

selPress$gene <- factor(selPress$gene, levels = node_names)
selPress$scen <- as.factor(selPress$scen)

# keep only first and last gen, and cast to long format 
firstLastgen <- allReplStats[allReplStats$generation %in% c(1, num_generations), ]
firstLastgen <- melt(firstLastgen, id.vars = c("scen", "gene", "generation"))
firstLastgen <- dcast(firstLastgen, id.vars = c("scen", "gene"), scen + gene ~ variable + generation)

sel <- firstLastgen[firstLastgen$scen == "sel", ]
sel <- sel[match(node_names, sel$gene), ]
neu <- firstLastgen[firstLastgen$scen == "neu", ]
neu <- neu[match(node_names, neu$gene), ]
firstLastgen <- rbind(sel, neu)
firstLastgen <- inner_join(firstLastgen, selPress, by = c("scen", "gene"))

# write summary of this network results
write.table(firstLastgen,
            file = "allReplEvolResults.txt",
            sep = "\t", row.names = FALSE)
