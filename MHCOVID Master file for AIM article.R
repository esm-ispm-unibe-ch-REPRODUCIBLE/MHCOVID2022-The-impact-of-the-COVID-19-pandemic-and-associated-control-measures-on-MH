
#************ MASTER FILE FOR MHCOVID ANALYSIS *************************************

# Set up environment------------------------------
library(meta)
library(tidyr)
library(tibble)
library(stringr)
library(readxl)
library(grid)
library(dplyr)
library(rjags)
library(R2jags)
library(Matrix)
library(knitr)
library(kableExtra)
library(rms) 
library(meta)
rm(list=ls())

# Load the data-----------------------------------
outcomes<- read_excel("outcomes.xlsx", na = "NA" )
metadata <- read_excel("metadata.xlsx", na = "NA" )

# Select only longitudinal studies, creates SMDs, and transforms variables that are needed in regressions
source("Script to process the MHCOVID data.R") 
# Output: saves the MHCOVIDdataset.csv (internally called "dataset")

##**************************************************************************************
# Analysis of change in scores -  pre versus during 
##**************************************************************************************

source("Script to perform Bayesian meta-analysis.R") 
# Creates Figure 2 (saved as .pdf file in your working directory) and funnel plots (see graphs window)
# The objects "prepostresults" contain the results from the Bayesian meta-analysis for depression and anxiety. 

source("Script to perform Bayesian meta-regressions.R") #for depression & anxiety and for psychological distress. 
# The objects saved with names "prepostmetaregNAMEOFVARIABLE" are the JAGS results from meta-regression for depression and anxiety. 

      
##**************************************************************************************
# Bayesian dose-response meta-analysis for Anxiety and Depression using study-specific coefficients 
# and fixed (common) effect between depression and anxiety. Creates the panels of Figure 3. 
##**************************************************************************************
# Note:
# The script runs 20,000 simulation for each of the four variables and will take a while to execute.
# The results in the article are based on 100,000 simulations; if you want to run 80,000 more simulations
# open the script and un-comment the commented lines. 

source("Script for dose-response.R")
# The script creates: 
# - the dose-response panels in figure 3 (saved as .pdf files in your working directory)
# - report about the estimated coefficients and heterogeneity parameters (saved as "beta and tau X.txt" files in your working directory)
# - the data to plot in the dose-effect figures (saved as .csv files in your working directory)
# - JAGS objects named XFEAD where X is the dose variable name (e.g. logDeathFEAD has the JAGS output from the 
# dose-response model with the log-cumulative number of deaths as dose variable)

