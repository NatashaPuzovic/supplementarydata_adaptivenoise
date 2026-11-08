#!/bin/bash

###### calculate network metrics
# in /results
bash ../scripts_generate_data_preprocess/getNetMetrics.sh 20230213_ER_40node 2000;

###### summarize GGPVF data from evolved populations
# in /results
# summarize the results of all replicates per net sample
bash ../scripts_generate_data_preprocess/calculateSelStrength_dirSelPos.sh 20230213_ER_40node 2000;
bash ../scripts_generate_data_preprocess/calculateSelStrength_dirSelPos_mutSBasal.sh 20230213_ER_40node 2000;
bash ../scripts_generate_data_preprocess/calculateSelStrength_dirSelPos_mutNoise.sh 20230410_ER_40node 2000;
time bash ../scripts_generate_data_preprocess/calculateSelStrength_dirSelPos_mutNoise.sh 20230411_ER_40node 2000;
bash ../scripts_generate_data_preprocess/calculateSelStrength_dirSelNeg.sh 20230213_ER_40node 2000;
bash ../scripts_generate_data_preprocess/calculateSelStrength_fluctSel.sh 20230213_ER_40node 2000;

###### combined network metrics and summarized GPF
# in /results
Rscript ../scripts_generate_data_preprocess/combineNetAndEvol_dirSelPos.R 20230213_ER_40node ER 0.05 2000; 
Rscript ../scripts_generate_data_preprocess/combineNetAndEvol_dirSelPos.R 20230401_ER_40node ER 0.05 2000;
Rscript ../scripts_generate_data_preprocess/combineNetAndEvol_dirSelPos.R 20230409_ER_40node ER 0.05 2000; #hardcoded 20230213 for netmetrics info
Rscript ../scripts_generate_data_preprocess/combineNetAndEvol_dirSelPos.R 20230410_ER_40node ER 0.05 2000; #hardcoded 20230213 for netmetrics info
Rscript ../scripts_generate_data_preprocess/combineNetAndEvol_dirSelPos.R 20230411_ER_40node ER 0.05 2000; #hardcoded 20230213 for netmetrics info
Rscript ../scripts_generate_data_preprocess/combineNetAndEvol_dirSelNeg.R 20230213_ER_40node ER 0.05 2000;
Rscript ../scripts_generate_data_preprocess/combineNetAndEvol_fluctSel.R 20230213_ER_40node ER 0.05 2000;
Rscript ../scripts_generate_data_preprocess/combineNetAndEvol_dirSelPos.R 20230401_ER_40node ER 0.05 2000;
