#-------------------------------------------------------------------------------
# File  : 03_run_bayesian_meta-regression.R
#
# Performs Bayesian meta-regression for Depression and Anxiety, and Psychological
# distress. Example usage is provided for a single covariate in both cases.
# N.B. the two analyses use different data preparation functions and models, as
# Psychological distress is a single condition, while Depression and Anxiety are
# assessed together.
#-------------------------------------------------------------------------------
# Example usage: Psychological distress ----------------------------------------
# Prepare the JAGS data for the Age covariate
data.prepost.jags.age <- MakeSingleConditionJagsData(data = data.prepost.psych,
                                                     studyid = record_id,
                                                     y = d,
                                                     var,
                                                     covlabel="cov",
                                                     corr=0.6,
                                                     covariate = pop0_central_age)
# Perform the Bayesian meta-regression
metareg.prepost.age <- jags(data = data.prepost.jags.age,
                            inits = NULL,
                            parameters.to.save = c("muMH", "tauMH", "gamma"),
                            model.file = prepostmetaregMHmodelRE1C,
                            n.chains = 3,
                            n.iter = 10000,
                            n.burnin = 1000,
                            DIC = FALSE,
                            n.thin = 2)
# Show the meta-regression output
print(metareg.prepost.age)
#-------------------------------------------------------------------------------
# The analysis can be performed as above for the following covariates, simply by
# replacing the `covariate` argument in the call to `MakeSingleConditionJagsData`
#
# pop0_central_age: central age within the study (mean or median, according to available data)
# pop0_percent_female: percentage of female participants in the study
# gdp_div_10000_minus_min: gross domestic product per capita in 2019, normalised by subtracting the lowest value
# gini_minus_min: GINI inequality index, normalised by subtracting the lowest value
# study_design_01: study design, 0 is longitudinal, 1 is cross-sectional at multiple timepoints
# rob_info_bias_num: risk of information bias
# rob_non_bias_num: risk of non-response bias
# rob_is_target_pop_num: risk of bias due to invited particpants matching the study's target population
# is_usa: binary indicating if the study was conducted in the Unites States
# is_china: binary indicating if the study was conducted in China
# days_after_first: number of days since the first case in the study country
# stringency: OxCGRT stringency measure
# log_confirmed_cumulative_100k: cumulative number of confirmed cases per 100,000 inhabitants (log-transformed)
# log_deaths_cumulative_100k: cumulative number of deaths per 100,000 inhabitants (log-transformed)
#-------------------------------------------------------------------------------
# TODO: check with Georgia what this is
# TODO: if required, move to 01_process_data.R
data.prepost.psych <- data.prepost.psych %>%
  mutate(OldBmeasurement = sqrt(-exact_days_after_first)) %>%
  mutate(OldBmeasurement = ifelse(is.na(OldBmeasurement), 0, OldBmeasurement))
# Example usage: Depression and Anxiety ----------------------------------------
# Prepare the JAGS data for the Age covariate
data.prepost.jags.age <- MakeMultipleConditionJagsData(data = data.prepost.daa,
                                                       studyid = record_id,
                                                       y = d,
                                                       var,
                                                       covlabel="cov",
                                                       corr = 0.6,
                                                       covariate = pop0_central_age)
metareg.prepost.age <- jags(data = data.prepost.jags.age,
                            inits = NULL,
                            parameters.to.save = c("muMH", "tauMH", "gamma"),
                            model.file = prepostmetaregMHmodelRE,
                            n.chains = 3,
                            n.iter = 10000,
                            n.burnin = 1000,
                            DIC = FALSE,
                            n.thin = 2)
# Show the meta-regression output
print(metareg.prepost.age)
