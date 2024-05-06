CONTENTS OF THIS FOLDER ——————————————

BM_JM_tutorial.R : A step-by-step implementation of BM-JM and the associated estimation and dynamic prediction procedures described in "Bayesian multivariate joint modeling of longitudinal, recurrent, and competing risk terminal events in patients with chronic kidney disease".

BM_JM_prepare_jags.R : Function to prepare the training data set to be passed onto JAGS.

BM_JM_jags.R : R file containing the JAGS model of the proposed BM-JM.

BM_JM_prepare_prediction.R : R file containing all the functions needed in BM-JM dynamic prediction procedures.

BM_JM_measure_prediction_accuracy.R : R file containing functions and codes to calculate dynamic AUC and BS from BM-JM predictions.

dat.long.RData, dat.recurrent.RData, dat.terminal.RData: Simulated multivariate training data set as described in "Bayesian multivariate joint modeling of longitudinal, recurrent, and competing risk terminal events in patients with chronic kidney disease".

pred.subject.RData: Simulated multivariate testing data set as described in "Bayesian multivariate joint modeling of longitudinal, recurrent, and competing risk terminal events in patients with chronic kidney disease".

INTRODUCTION ——————————————

The contents of this folder allow for implementation of the BM-JM estimation and inference, as well as dynamic prediction and predictive accuracy measurement procedures described in "Bayesian multivariate joint modeling of longitudinal, recurrent, and competing risk terminal events in patients with chronic kidney disease". Users can implemente the aforementioned procedures on the attached simulated sample training and testing data set. Detailed instructions on how to perform the aforementioned procedures and visualize results are included in BM_JM_tutorial.R.

REQUIREMENTS ——————————————

The included R programs require R 4.1.0 (R Core Team, 2021) and the packages listed in BM_JM_tutorial.R. 

INSTALLATION ——————————————

Load the R program files into the global environment and install required packages using commands in BM_JM_tutorial.R, as well as pre-install JAGS software.
