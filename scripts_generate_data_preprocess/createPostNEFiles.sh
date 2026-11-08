#!/bin/bash

# create params for mutation-selection balance (needs netEst.pop)
Rscript scripts_and_params/createPrmsV_forMutSelBalance.R 2000

# create first gen file for mutation-selection balance (needs netEst.pop)
Rscript scripts_and_params/createFirstGenPopFile_forMutSelBalance.R 2000
