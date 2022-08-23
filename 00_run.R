#-------------------------------------------------------------------------------
# File  : 00_run.R
#
# Performs data processing, re-creates figures 2a and 2b from the manuscript,
# provides example usage of Bayesian meta-regression functions, and performs
# dose-response analysis.
#-------------------------------------------------------------------------------
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

# Load the data ----------------------------------------------------------------
outcomes<- read_excel("data/outcomes.xlsx", na = "NA" )
metadata <- read_excel("data/metadata.xlsx", na = "NA" )

#-------------------------------------------------------------------------------
# Select only longitudinal studies, creates SMDs, and transforms variables that
# are needed in regressions.
#
# Processed data are saved to `mhcovid_dataset.csv`
#-------------------------------------------------------------------------------
source("01_process_data.R")

#-------------------------------------------------------------------------------
# Analysis of change in scores -- pre versus during
#-------------------------------------------------------------------------------
source("02_run_meta-analysis.R")
# Creates Figure 2 (saved as .pdf file in your working directory) and funnel plots (see graphs window)
# The objects "prepostresults" contain the results from the Bayesian meta-analysis for depression and anxiety.

source("03_run_bayesian_meta-regression.R")
# The objects saved with names "prepostmetaregNAMEOFVARIABLE" are the JAGS results from meta-regression for depression and anxiety.

#-------------------------------------------------------------------------------
# Bayesian dose-response meta-analysis for Anxiety and Depression using study-specific coefficients
# and fixed (common) effect between depression and anxiety. Creates the panels of Figure 3.
#-------------------------------------------------------------------------------
# Note:
# The script runs 20,000 simulation for each of the four variables and will take a while to execute.
# The results in the article are based on 100,000 simulations; if you want to run 80,000 more simulations
# open the script and un-comment the commented lines.

source("04_run_dose_response_models.R")
# The script creates:
# * the dose-response panels in figure 3 (saved as .pdf files in your working directory)
# * report about the estimated coefficients and heterogeneity parameters (saved as "beta and tau X.txt" files in your working directory)
# * the data to plot in the dose-effect figures (saved as .csv files in your working directory)
# * JAGS objects named XFEAD where X is the dose variable name (e.g. logDeathFEAD has the JAGS output from the
# * dose-response model with the log-cumulative number of deaths as dose variable)

