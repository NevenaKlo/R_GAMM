# R_GAMM
This repository contains the R script used to analyse data from a visual world eye-tracking experiment using a Generalized Additive Mixed Model.

GAMMs are a type of regression analysis, but the models use a different function compared to linear regression. GAMMs use a smooth function of the independent variable(s) and, as opposed to the fixed coefficient β of linear regression models, the function f can vary across the range of x (where x = e.g., time). 

The imported data set contains information about participants, trials, items, conditions, and fixations to Target. The aim of the analysis is to detect differences in patterns across conditions.

In this script, a Gaussian distribution is used to take into account autocorrelation. We thus include a parameter in the model determined by checking for autocorrelation in a plain model of the data (i.e., a model in which we do not include the relevant independent variables). In order to use a Gaussian distributin, the binary data (1 “looks to target” vs 0 “looks to competitor”) were logit-transformed. Logit is another way of expressing probability (as is “proportions of looks to Target”). The logit of a probability is the log odds of (p / 1-p).  

The package mgcv is used to fit a model which estimates the effect of the independent variable "condition" on the logit of looks to Target, with random smooths for participants and items.

The package itsadug is used to plot the results and run pairwise comparisons.

