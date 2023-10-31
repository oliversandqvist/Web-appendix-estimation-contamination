# Web appendix

## Introduction
This repository contains code and output to reproduce results in the paper [Estimation for multistate models subject to reporting delays
and incomplete event adjudication]() by Kristian Buchardt, Christian Furrer, and Oliver Lunding Sandqvist. The subdirectory *Data application* contains the data set LEC-DK19 collected by a large Danish insurer. A total of 416,483 insured are included across five tables concerning disabilities (**disab**), reactivations (**reac**), disability delays (**delay**), disability adjudications (**adjDisab**), and reactivation adjudications (**adjReac**). 

## Code structure

* [Data application/LECDK19.RData](<Data application/LECDK19.RData>): The LEC-DK19 data set.
* [Data application/estimate.R](<Data application/estimate.R>): Load the LEC-DK19 data, run the proposed and naive estimation procedure, plot the results, run the bootstrap procedure, and save the results. This code is used to generate results in Section 6 and Appendix G of the paper. The output is the following.
    * [Data application/Results/paramEst.Rda](<Data application/Results/paramEst.Rda>): Estimated parameter values.
    * [Data application/Results/paramBoot.Rda](<Data application/Results/paramBoot.Rda>): Estimated bootstrap parameter values.
    * [Data application/Figures/AdjRates.png](<Data application/Figures/AdjRates.png>): Plot of fitted adjudication hazards.
    * [Data application/Figures/DelayRates.png](<Data application/Figures/DelayRates.png>): Plot of fitted reverse time hazards for disability reporting delays.
    * [Data application/Figures/FittedRates.png](<Data application/Figures/FittedRates.png>): Plot of fitted biometric hazards.
* [Numerical study/run_sim.R](<Numerical study/run_sim.R>): Simulates data, repeatedly run the proposed, Poisson approximation and naive estimation procedure, and save the results. Dependencies and output are as follows.
    * [Numerical study/Simulation/semi_markov_sim.R](<Numerical study/Simulation/semi_markov_sim.R>): Functions to simulate biometric multistate model.
    * [Numerical study/Simulation/semi_markov_sim_adj.R](<Numerical study/Simulation/semi_markov_sim_adj.R>): Functions to simulate adjudication multistate model.
    * [Numerical study/Results/dfParamEst.Rda](<Numerical study/Results/dfParamEst.Rda>): Parameter estimates.
    * [Numerical study/Misc/ParamTableAndPlot.R](<Numerical study/Misc/ParamTableAndPlot.R>): Script to calculate mean of parameters, standard deviations, and plot histograms.
    * [Numerical study/Figures/Histogram.png](<Numerical study/Figures/Histogram.png>): Histogram of parameter estimates.
* [Numerical study/run_sim_onlyTheta7.R](<Numerical study/run_sim_onlyTheta7.R>) and [Numerical study/run_sim_boot.R](<Numerical study/run_sim_boot.R>): Repeatedly simulates data and runs the estimation procedure for θ<sub>7</sub> and similarly simulates bootstrap samples and estimates θ<sub>7</sub>. Dependencies and output are as follows.
    * [Numerical study/Simulation/semi_markov_sim.R](<Numerical study/Simulation/semi_markov_sim.R>) and [Numerical study/Simulation/semi_markov_sim_adj.R](<Numerical study/Simulation/semi_markov_sim_adj.R>): See above.
    * [Numerical study/Results/dfParamTheta7.Rda](<Numerical study/Results/dfParamTheta7.Rda>) and [Numerical study/Results/dfParamBoot.Rda](<Numerical study/Results/dfParamBoot.Rda>): Parameter and bootstrap estimates.
    * [Numerical study/Misc/ComputeCIBoot.R](<Numerical study/Misc/ComputeCIBoot.R>): Script to calculate bootstrap confidence intervals.

## Data description

The variables included in LEC-DK19 are:

**disab**   
[,1] Gender (gender)  
[,2] Age (age)  
[,3] Occurrence of disability event (disabOcc)  
[,4] Exposure (expo)  
[,5] Duration at the time of analysis since the disability event occurred (durDisab)  
[,6] Duration at the time of analysis since the disability event was reported (durDisabReport)  
[,7] Date at the beginning of the observation (dateStart)  
[,8] Adjudication state at the time of analysis (adjState)  
[,9] Indicator for whether the reported disability has been rejected before at the time of analysis (rejectedBefore)  
[,10] Id of the insured (id)  

**reac**  
[,1] Gender (gender)  
[,2] Age (age)  
[,3] Occurrence of reactivation event (reacOcc)  
[,4] Exposure (expo)  
[,5] Date at the beginning of the observation (dateStart)  
[,6] Duration in the disabled state (durDisab)  
[,7] Adjudication state at the time of analysis (adjState)  
[,8] Duration at the time of analysis since the reactivation event occurred (durReac)    
[,9] Id of the insured (id)  

**delay**  
[,1] Age (age)  
[,2] Gender (gender)  
[,3] Exposure (expo)  
[,4] Adjudication state at the time of analysis (adjState)  
[,5] Duration at the time of analysis since the disability event was reported (durDisabReport)   
[,6] Maximal reporting delay that could have been observed (repDelayMax)  
[,7] Indicator for whether the reported disability has been rejected before at the time of analysis (rejectedBefore)  
[,8] Duration in the disabled state (durDisab)  
[,9] Observed reporting delay (repDelayDisab)  
[,10] Id of the insured (id)  

**adjDisab**  
[,1] Age (age)  
[,2] Gender (gender)  
[,3] Exposure for reverse time hazard (expo)  
[,4] Adjudication state (adjState)  
[,5] Duration since the disability event occurred (durDisab)  
[,6] Duration since the disability event was reported (durDisabReport)  
[,7] Occurrence of awarding of disability benefits (awardOcc)   
[,8] Occurrence of rejection of disability claim (rejectOcc)   
[,9] Occurrence of reapplication of disability claim (reapplyOcc)   
[,10] Occurrence of death event (deadOcc)  
[,11]  Indicator for whether the reported disability has been rejected before (rejectedBefore)   
[,12] Id of the insured (id)  

**adjReac**  
[,1] Adjudication state (adjState)  
[,2] Gender (gender)   
[,3] Age (age)   
[,4] Duration since the reactivation event occurred (durReac)    
[,5] Duration since the disability event occurred (durDisab)   
[,6] Exposure (expo)  
[,7] Occurrence of awarding of disability benefits (awardOcc)   
[,8] Occurrence of reapplication of disability claim (reapplyOcc)   
[,9] Occurrence of rejection of disability claim (rejectOcc)   
[,10] Occurrence of death event (deadOcc)  
[,11] Id of the insured (id)  
