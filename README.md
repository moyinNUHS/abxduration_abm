[![License: CC0-1.0](https://licensebuttons.net/l/zero/1.0/88x31.png)](http://creativecommons.org/publicdomain/zero/1.0/) ![R](https://img.shields.io/badge/R-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)

This folder contains codes to implement and assess three increasingly complex stochastic agent based models to determine if shortening treatment duration is an effective strategy in reducing the prevalence of AMR, and under what conditions does shortening duration most effectively reduces AMR. 

The three models are named: simple 3-state, co-carriage 5-state and population growth models.

The codes are organised as below: 

## MODELS
Contains codes for the models: `model_simple3state`, `model_cocarriage5state`, `model_populationgrowth`
And related functions that are required to run the models 

The types of output that can be calculated from each code file are labelled in the file names


## RUN_MODELS
Contains codes to run the models stored in the `/MODELS` folder 

Include checks for agreement between runs to decide if sample size ie size of latin hypercube is adequate, using Symmetric Best Measure of Agreement (SBMA) between the PRCC coefficients of two runs with different sample sizes.


## RUNS
Contains results of the runs from each model and scenario

## PLOTS
Codes to generate plots in the manuscript

## AA TESTS
Codes to run Vargha-Delaney A-Test to find the optimal number of iterations per simulation. 
The test compares distributions of simulation outputs under identical parameter values, which is a non-parametric measure of the difference between the distributions of model outputs suggesting if the outputs are consistent.

## UNIT TESTS
Contains code which tests each step of each model to ensure they give the expected outputs


