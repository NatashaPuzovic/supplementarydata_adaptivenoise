#!/bin/bash
baseFolder=$1
mkdir sim_data/${baseFolder}
mkdir sim_data/${baseFolder}/scripts_and_params

generatedNetsFolder=$2
numNets=$3

# run Rscript in each net folder
for ((netNum = 1; netNum <= $numNets; netNum++))
do
	mkdir sim_data/${baseFolder}/"net_$netNum"
	cp sim_data/generated_net_topos/${generatedNetsFolder}/"rmtrxFoundingNet_$netNum.txt" sim_data/${baseFolder}/"net_$netNum"/rmtrxFoundingNet
	echo "made net folder $netNum"
done
