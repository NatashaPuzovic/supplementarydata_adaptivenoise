#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# This is a script to concatenate the results of all network replicates

# Author: Natasha Puzovic (puzovic@evolbio.mpg.de)
# Date: June 2020
# Modified: Sept 2020


# test if all the arguments are present: if not, return an error
if (length(args) != 4) {
  stop("Please specify results folder; network topology type, e.g. ER, BA, WS; density; number of networks.", call. = FALSE)
} else if (length(args) == 4) {
  resultsFolder = args[1]
  topo = args[2]
  dens = args[3]
  numNets = args[4]
}

setwd(paste(resultsFolder))

# make df
allNetsResults <- data.frame(matrix(nrow = 0, ncol = 0))

for(netNum in 1:numNets) {

  print(paste(paste("net_", netNum, sep = "")))
  netFolder = paste("net_", netNum, sep = "")
  setwd(netFolder)

  if(!file.exists("allReplEvolResults_fluctSel.txt"))
  {
    print(paste(paste("MISSING RESULTS net_", netNum, sep = "")))
    setwd("..")
    next
  }
  
  # read preprocessed results of one net
  allReplEvolResults <- read.table(file = "allReplEvolResults_fluctSel.txt", sep = "\t", header = TRUE)
  #allReplEvolResults <- read.table(file = "allReplEvolResults_dirSelNeg.txt", sep = "\t", header = TRUE)
  #allReplEvolResults <- read.table(file = "allReplEvolResults_fluctSel.txt", sep = "\t", header = TRUE)
  
  # read net data
  netMetrics <- read.table(file = "netMetrics.txt", sep = " ", header = TRUE)
  rownames(netMetrics) <- NULL

  # add net number & net data
  allReplEvolNet = as.data.frame(cbind(rep(topo, dim(allReplEvolResults)[1]),
                                       rep(netNum, dim(allReplEvolResults)[1]),
                                       allReplEvolResults,
                                       netMetrics[, -which(names(netMetrics) == "node_names")]))
  colnames(allReplEvolNet)[1:2] = c("topo", "net")

  # add this net results to allnet df
  allNetsResults = rbind(allNetsResults, allReplEvolNet)
  setwd("..")

}

# read params
params <- read.table(file = "net_1/prmsS_fluctSel_rep1",
#params <- read.table(file = "net_1/prmsS_fluctSel_rep1",
#params <- read.table(file = "net_1/prmsS_fluctSel_rep1",
                     sep = " ", header = FALSE, col.names = c("param", "val"))

# add params of this sim set to results
allNetsResults$num_generations =  params$val[params$param == "NUM_GENERATIONS"]
allNetsResults$num_nodes = params$val[params$param == "NUM_NODES"]
allNetsResults$dens = as.numeric(dens)
allNetsResults$pop_size = params$val[params$param == "POP_SIZE"]


# write summary of all network results (whole sim batch)
write.table(allNetsResults,
            #file = "allNetsResults_dirSelPos.txt",
            #file = "allNetsResults_dirSelNeg.txt",
            file = "allNetsResults_fluctSel.txt",
            sep = "\t", row.names = FALSE)

######################################################### end ############################################################
