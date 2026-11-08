#!/bin/bash
clustRunFolder=$1
cd "${clustRunFolder}"

# preprocess replicate results for all networks
numNets=$2

# run Rscript in each net folder
for ((netNum = 1; netNum <= $numNets; netNum++))
do
	echo "working on net $netNum"
	cd "net_$netNum"
	Rscript ../../../scripts_generate_data_preprocess/selStrength_dirSelNeg.R
	cd ..
done

# output:
# in each net folder - allReplEvolResults.txt
