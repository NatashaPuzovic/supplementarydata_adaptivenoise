#!/bin/bash
clustRunFolder=$1
numNets=$2

cd "${clustRunFolder}"

# write net metrics for all networks


# run Rscript in each net folder
for ((netNum = 1; netNum <= $numNets; netNum++))
do
	echo "working on net $netNum"
	cd "net_$netNum"
	Rscript ../../../scripts_generate_data_preprocess/getNetMetricsOneNet.R
	cd ..
done

# output:
# in each net folder - netMetrics.txt
