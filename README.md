# Supplementary data: evolution of expression noise in gene regulatory networks
This repository contains the code to reproduce the statistical analysis of the following publication (https://doi.org/10.1101/2023.12.04.569843):
"Network-level selection constrains the adaptive evolution of gene expression noise in regulatory genes in changing environments"

Evolutionary simulations were used to study the evolution of gene-specific expression noise and mean expression level in gene regulatory networks.

## Installation and data download
Clone this repo:

```
$ git clone git@github.com:NatashaPuzovic/supplementarydata_expressionnoise.git
$ cd supplementarydata_adaptivenoise
```

Download the zipped raw data from zenodo (https://zenodo.org/records/24049819) in this folder and unzip.

```
$ unzip results.zip
```

## Reproducing the results
The full code to reproduce the results in the maintext, supplementary information, and additional tests is in the scripts_analysis folder. The code to generate the input files for the simulations and preprocess the raw simulation output is in the scripts_generate_data_preprocess folder.

Running the scripts will output figures in the ./plots folder.

## Raw simulation results
Only preprocessed simulation results are available in the zenodo repository due to the large size of the raw simulations results. Please contact me directly if you want the raw simulation results (~1TB). 

## Simulator
The code used for running the simulations presented in the paper is not included here, but is available in this repo: https://github.com/NatashaPuzovic/netlings. The version of the program used in this study is tagged as v2.1.1 Release.