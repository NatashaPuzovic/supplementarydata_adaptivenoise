#!/bin/bash

###### calculate network metrics
# in /results
bash ../scripts_generate_data_preprocess/getNetMetrics.sh 20230213_ER_40node 2000;

###### summarize GGPVF data from evolved populations
# in /results
# summarize the results of all replicates per net sample
bash ../scripts_generate_data_preprocess/calculateSelStrength_dirSelPos.sh 20230213_ER_40node 2000; # dir sel pos
bash ../scripts_generate_data_preprocess/calculateSelStrength_dirSelNeg.sh 20230213_ER_40node 2000; # dir sel neg
bash ../scripts_generate_data_preprocess/calculateSelStrength_fluctSel.sh 20230213_ER_40node 2000; # fluct sel

###### combined network metrics and summarized GPF
# in /results
Rscript ../scripts_generate_data_preprocess/combineNetAndEvol_dirSelPos.R 20230213_ER_40node ER 0.05 2000; # dir sel pos
Rscript ../scripts_generate_data_preprocess/combineNetAndEvol_dirSelNeg.R 20230213_ER_40node ER 0.05 2000; # dir sel neg
Rscript ../scripts_generate_data_preprocess/combineNetAndEvol_fluctSel.R 20230213_ER_40node ER 0.05 2000; # fluct sel
