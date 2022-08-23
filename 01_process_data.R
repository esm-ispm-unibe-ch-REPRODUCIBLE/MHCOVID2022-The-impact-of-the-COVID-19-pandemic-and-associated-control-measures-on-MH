#-------------------------------------------------------------------------------
# File  : 01_process_data.R
#
# Processes the clean (meta-) data files, adding new variables and calculating
# values in order to facilitate subsequent analyses.
#-------------------------------------------------------------------------------
# Rename components of the OxCGRT stringency index with human-readable names
outcomes <- outcomes %>% rename(
  ESI = EconomicSupportIndex,
  CHIndex = ContainmentHealthIndex,
  School_closing = C1,
  Workplace_closing = C2,
  Stay_home_req = C6,
  Restrictions_movement = C7,
  Facial_cover = H6
)

# Ensure timepoint is a date
outcomes$timepoint <- as.Date(outcomes$timepoint,format = "%m/%d/%y")

# Create population and subgroup variables -------------------------------------
outcomes$study_population[outcomes$study_population=="Women,Adults"]<-"Adults"
outcomes$population[outcomes$population=="Women,Adults"]<-"Women"

# Indicate if the analysed population matches the study population (as opposed to subgroup)
outcomes <- outcomes %>%
  mutate(stpop = study_population == population)
# Outcome population to analyse
outcomes <- outcomes %>%
  mutate(out_pop = ifelse(stpop == 0 & (population == "Men" | population == "Women"),
                          paste(population, study_population),
                          population))
# Create sex variable
outcomes$sex<-"Mixed"
outcomes$sex[grepl("Women",outcomes$out_pop)]<-"Women"
outcomes$sex[grepl("Men",outcomes$out_pop)]<-"Men"


# Create a dataset with only continuous, longitudinal studies ------------------
data <- outcomes %>%
  mutate(is_continuous_outcome = !is.na(sd)) %>%  # continuous outcomes
  filter(is_continuous_outcome == 1) %>%
  mutate(is_main_analysis = population == study_population) %>%
  filter(is_main_analysis == 1) %>%  # studies for main analysis, not subgroups
  filter(is_longitudinal == 1) %>%  # studies marked as longitudinal
  group_by(record_id) %>%
  filter(n_distinct(timepoint) > 1) %>%  # studies with 2 or more unique timepoints
  left_join(metadata, by = c("record_id")) %>%  # attach metadata
  select(!ends_with(".y")) %>%  # remove metadata columns that were already present
  filter(record_id != 66150) %>%  # this study has different instruments per timepoint
  arrange(record_id, condition, scale, population)

# Calculate standardised mean differences --------------------------------------
source("calculate_smd_vectorised.R")
data <- data %>%
  filter(study_design != 0) %>%
  group_by(record_id, scale, population) %>%
  filter(n() > 1) %>%
  ungroup()

smd <- calculate_smd(data) %>%
  ungroup() %>%
  select(record_id, s2, d, var, starts_with("cov"))

# Keep useful columns, and join SMD data
data <- data %>%
  rename(country = country.x,
         rob_info_bias = rob_info_bias.x,
         num_invited = num_invited.x,
         num_assessed = num_assessed.x,
         response_rate = response_rate.x,
         rob_is_target_pop = rob_is_target_pop.x,
         rob_non_bias = rob_non_bias.x) %>%
  select(record_id, doi,author_1,year,population,
         study_population, country, condition, scale, timepoint,
         n_timepoints, is_longitudinal, is_prepandemic, sample_size,
         is_binary, score, sd, y, sey, stringency, days_after_first,
         days_after_pandemic,confirmed_cumulative, deaths_cumulative,
         confirmed_avg, deaths_avg, cumulative_stringency,country_population_2019,
         confirmed_per_100000, deaths_per_100000,confirmed_avg_per_100000, deaths_avg_per_100000,
         gdp_per_capita_2019, gini_2019, effect_direction, rob_info_bias, num_invited,
         num_assessed, response_rate, rob_is_target_pop, rob_non_bias,
         stpop, out_pop, sex, is_main_analysis, is_continuous_outcome, pop0_name,
         tag_num_timepoints, study_design, method_recruitment, method_collection,
         multi_country, pop0_central_age, pop0_sample_size, pop0_mean_med,
         pop0_sd, pop0_min_age, pop0_max_age, pop0_percent_female,
         pop_percent_phys_con, pop_percent_psych_con, pop_percent_covid19,
         pop_ethnicity, pop_ethnicity_other, ESI, CHIndex, School_closing,
         Workplace_closing, Stay_home_req, Restrictions_movement, Facial_cover)
data.smd <- cbind(data, select(smd, -c("record_id")))
data.smd <- data.smd %>%
  arrange(record_id, condition, scale, timepoint) %>%
  filter(!(record_id == 121277 & (condition == "Social anxiety" | condition == "Panic/Somatic symptoms")))

# Create covariates ------------------------------------------------------------
anxiety.conditions <- c("Social anxiety", "Generalized Anxiety Disorder",
                        "Panic/somatic symptoms", "Panic disorder")
dataset <- data.smd %>%
  mutate(effect_direction = ifelse(is.na(effect_direction), 1, effect_direction)) %>%
  mutate(d = effect_direction * d) %>%
  mutate(scaled_score = score / s2) %>%
  mutate(study_design_01 = study_design - 1) %>%
  # Transform existing covariates for specific analyses
  mutate(exact_days_after_first = days_after_first) %>%
  mutate(days_after_first = ifelse(days_after_first < 0, 0, days_after_first)) %>%
  mutate(condition = ifelse(condition %in% anxiety.conditions, "Anxiety", condition)) %>%
  mutate(gini_minus_min = gini_2019 - min(gini_2019)) %>%
  mutate(gdp_div_10000_minus_min = gdp_per_capita_2019 / 10000 - min(gdp_per_capita_2019 / 10000)) %>%
  mutate(log_cumulative_stringency = log(cumulative_stringency + 1)) %>%
  mutate(sqrt_cumulative_stringency = sqrt(cumulative_stringency)) %>%
  mutate(logconfirmed_cumulative100k = log(confirmed_cumulative / country_population_2019 * 100000 + 1)) %>%
  mutate(logdeaths_cumulative100k = log(deaths_cumulative / country_population_2019 * 100000 + 1)) %>%
  mutate(sqrtconfirmed_cumulative100k = sqrt(confirmed_cumulative / country_population_2019 * 100000)) %>%
  mutate(sqrtdeaths_cumulative100k = sqrt(deaths_cumulative / country_population_2019 * 100000)) %>%
  mutate(days_after_pandemic = ifelse(is.na(days_after_pandemic), days_after_first, days_after_pandemic)) %>%
  # Indicators for study's country
  mutate(is_china = country == "China") %>%
  mutate(is_usa = country == "United States") %>%
  # For stringency measures, pre-pandemic timepoints should be 0 and not NA
  mutate(ESI = ifelse(is_prepandemic == 1, 0, ESI)) %>%
  mutate(CHIndex = ifelse(is_prepandemic == 1, 0, CHIndex)) %>%
  mutate(School_closing = ifelse(is_prepandemic == 1, 0, School_closing)) %>%
  mutate(Workplace_closing = ifelse(is_prepandemic == 1, 0, Workplace_closing)) %>%
  mutate(Restrictions_movement = ifelse(is_prepandemic == 1, 0, Restrictions_movement)) %>%
  mutate(Stay_home_req = ifelse(is_prepandemic == 1, 0, Stay_home_req)) %>%
  mutate(Facial_cover = ifelse(is_prepandemic == 1, 0, Facial_cover))

# Recode RoB responses as numerics
rob.map <- data.frame(value = c(0:2), row.names = c("Low risk", "Unclear risk", "High risk"))
rob <- dataset %>%
  select(c(record_id, condition, sample_size, rob_info_bias, rob_is_target_pop, rob_non_bias)) %>%
  mutate(rob_info_bias_num = rob.map[rob_info_bias, "value"]) %>%
  mutate(rob_is_target_pop_num = rob.map[rob_is_target_pop, "value"]) %>%
  mutate(rob_non_bias_num = rob.map[rob_non_bias, "value"]) %>%
  ungroup() %>%
  group_by(record_id, condition) %>%
  # set the risk of bias as the maximum entered for any condition within a study
  group_modify(~ mutate(.x, rob_info_bias_num = max(rob_info_bias_num))) %>%
  group_modify(~ mutate(.x, rob_non_bias_num = max(rob_non_bias_num))) %>%
  group_modify(~ mutate(.x, rob_is_target_pop_num = max(rob_is_target_pop_num))) %>%
  # get mean sample size within study
  group_modify(~ mutate(.x, mean_study_sample_size = mean(sample_size))) %>%
  ungroup() %>%
  select(c(ends_with("num"), mean_study_sample_size))
# Bind back to original data.frame
dataset <- cbind(dataset, rob)

# Correct some ASCII character encoding errors (likely caused by Excel)
dataset$author_1[dataset$author_1=="Br√§scher"]<-"Braescher"
dataset$author_1[dataset$author_1=="Gim√©nez-Das√≠"]<-"Gimenez-Dasi"
dataset$author_1[dataset$author_1=="Sch√§fer"]<-"Schaefer"
dataset$author_1[dataset$author_1=="Sch√ºtzwohl"]<-"Schuetzwohl"
dataset$author_1[dataset$author_1=="S√∏nderskov"]<-"Sonderskov"
dataset <- dataset %>%
  mutate(author_year = paste(author_1, year))

# Write clean dataset to disk --------------------------------------------------
write.csv(dataset, "mhcovid_dataset.csv" )
