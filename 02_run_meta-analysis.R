#-------------------------------------------------------------------------------
# File  : 02_run_meta-analysis.R
#
# Takes clean, processed outcome data and performs meta-analysis. Creates funnel
# and forest plots of the results, creating figures 2a and 2b from the paper.
#-------------------------------------------------------------------------------
# Perform frequentist pre-during analysis for all conditions
data <- read.csv("data/dataset_alex.csv")
data.prepost <- data %>%
  group_by(record_id, condition) %>%
  mutate(sum_prepandemic = sum(is_prepandemic)) %>%
  filter(sum_prepandemic > 0) %>%
  # TODO: calculate these in 01_process_data.R
  mutate(se = sqrt(var)) %>%
  mutate(log_mean_study_sample_size = log(mean_study_sample_size))

# Subset the data by conditions, to create separate plots
data.prepost.daa <- data.prepost %>%
  filter(condition %in% c("Anxiety", "Depression"))
data.prepost.psych <- data.prepost %>%
  filter(condition == "Psychological distress")
conditions <- c("Psychological distress", "Mental Wellbeing",
                "Sleep disturbance", "Internet Gaming Disorder",
                "Problematic smartphone-application use", "Life satisfaction",
                "Problematic social media use", "Somatoform disorder")
data.prepost.other <- data.prepost %>%
  filter(condition %in% conditions)

# Meta-analysis for Depression and Anxiety
meta.prepost.daa <- metagen(TE = d,
                            seTE = se,
                            studlab = author_year,
                            data = data.prepost.daa,
                            subset = (is_prepandemic == 0),
                            subgroup = condition,
                            overall = TRUE,
                            fixed = FALSE,
                            prediction = TRUE,
                            prediction.subgroup = TRUE)
# Meta-analysis for other conditions
meta.prepost.other <- metagen(TE = d,
                              seTE = se,
                              studlab = author_year,
                              data = data.prepost.other,
                              subset = (is_prepandemic == 0),
                              subgroup = condition,
                              overall = FALSE)
# Meta-analysis for psychological distress
meta.prepost.psych <- metagen(TE = d,
                              seTE = se,
                              studlab = author_year,
                              data = data.prepost.psych,
                              subset = (is_prepandemic == 0),
                              subgroup = condition,
                              overall = FALSE)

# Draw funnel plot for anxiety and depression
funnel(meta.prepost.daa, contour = c(0.9, 0.95, 0.99), xlab="SMD anxiety and depression")
legend(-0.50, 0.02,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "< 0.01"),
       fill = c("darkgray", "gray", "lightgray"))
# Draw funnel plot for psychological distress
funnel(meta.prepost.psych,contour = c(0.9, 0.95, 0.99),xlab="SMD psychological distress")
legend(-0.50, 0.02,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "< 0.01"),
       fill = c("darkgray", "gray", "lightgray"))

# Create forest plot of pre-post pandemic meta-analysis in other conditions (Figure 2b)
pdf(paste("Figure 2b Pre-post forest plot other conditions", ".pdf", sep=""), width = 11, height = 11)
cols.forest <- c("studlab", "country", "population", "sample_size",
                 "pop0_percent_female", "rob_info_bias", "rob_non_bias",
                 "days_after_first")
labs.forest <- c("Study", "Country", "Population", "Sample size ", "Fraction females ",
                 "Information bias ", "Non-reponse bias ", "Days in pandemic")
forest(meta.prepost.other,
       test.overall = FALSE,
       comb.fixed = FALSE,
       print.byvar = TRUE,
       overall  = FALSE,
       weight.study = "same",
       title = "Post-Pre SMD",
       pooled.totals = FALSE,
       xlab = paste("SMD"),
       smlab = "SMD",
       print.pval.Q = FALSE,
       plotwidth = "6cm",
       print.I2 = FALSE,
       lwd = 1,
       fontsize = 7,
       fs.test.effect.subgroup = 0,
       subgroup = FALSE,
       fs.hetstat = 7,
       just = "left",
       addrow = FALSE,
       leftcols = cols.forest,
       leftlabs = labs.forest,
       colgap = "0.01cm",
       col.by = "black",
       hetlab = "",
       calcwidth.subgroup = FALSE,
       overall.hetstat = FALSE,
       hetstat = FALSE,
       text.random.w = sapply(meta.prepost.other$k.w, paste, "observations"),
       sort = days_after_first)
dev.off()

# Perform Bayesian pre-during analysis -----------------------------------------
source("Make pre-post Bayesian meta-analysis functions.R")

data.prepost.jags <- MakePrepostJagsData(data = data.prepost.daa,
                                         studyid = record_id,
                                         y = d,
                                         var = var,
                                         covlabel = "cov",
                                         corr = 0.6)
data.prepost.jags$Y1 <- data.prepost.jags$Y
prepost.jags.results <- jags(data=data.prepost.jags,
                             inits=NULL,
                             parameters.to.save = c("mu","muMH", "predSMD", "tau", "tauMH","predSMDMH"),
                             model.file = prepostMHmodelRE,
                             n.chains = 3,
                             n.iter = 50000,
                             n.burnin = 10000,
                             DIC = FALSE,
                             n.thin = 3)

# Extract useful metrics from JAGS output
summary.t <- as.data.frame(t(prepost.jags.results$BUGSoutput$summary))
mus <- as.data.frame(t(select(summary.t, starts_with("mu"))))
taus <- as.data.frame(t(select(summary.t, starts_with("tau"))))
pred.smds <- as.data.frame(t(select(summary.t, starts_with("predSMD"))))
prepost.results <- rbind(mus, pred.smds, taus) %>%
  select(`50%`, `2.5%`, `97.5%`) %>%
  rename("Posterior Median"=`50%`, "lowCI"=`2.5%`, "highCI" =`97.5%`)
rm(summary.t, mus, taus, pred.smds)

# Modify metagen output object such that forest is able to correctly plot it
meta.prepost.daa$TE.random.w     <- prepostresults[1:2, 1]
meta.prepost.daa$lower.random.w  <- prepostresults[1:2, 2]
meta.prepost.daa$upper.random.w  <- prepostresults[1:3, 3]
meta.prepost.daa$TE.random       <- prepostresults[3, 1]
meta.prepost.daa$lower.random    <- prepostresults[3, 2]
meta.prepost.daa$upper.random    <- prepostresults[3, 3]
meta.prepost.daa$lower.predict.w <- prepostresults[4:5, 2]
meta.prepost.daa$upper.predict.w <- prepostresults[4:5, 3]
meta.prepost.daa$lower.predict   <- prepostresults[6, 2]
meta.prepost.daa$upper.predict   <- prepostresults[6, 3]
meta.prepost.daa$tau.w           <- prepostresults[7:8, 1]
meta.prepost.daa$lower.tau.w     <- prepostresults[7:8, 2]
meta.prepost.daa$upper.tau.w     <- prepostresults[7:8, 3]
meta.prepost.daa$tau             <- prepostresults[9, 1]
meta.prepost.daa$lower.tau       <- prepostresults[9, 2]
meta.prepost.daa$upper.tau       <- prepostresults[9, 3]

pdf(paste("Figure 2a Pre-post forest plot Anxiety and Depression",".pdf",sep=""), width = 11, height = 9)
cols.forest <- c("studlab", "country", "population", "sample_size",
                 "pop0_percent_female", "rob_info_bias", "rob_non_bias",
                 "days_after_first")
labs.forest <- c("Study", "Country", "Population", "Sample size ", "Fraction females ",
                 "Information bias ", "Non-reponse bias ", "Days in pandemic")
forest(meta.prepost.daa,
       test.overall = FALSE,
       prediction = TRUE,
       comb.common = FALSE,
       print.byvar = TRUE,
       overall  = TRUE,
       title = "Post-Pre SMD",
       pooled.total = FALSE,
       xlab = paste("SMD"),
       smlab = "SMD",
       print.pval.Q = FALSE,
       plotwidth = "6cm",
       print.I2 = FALSE,
       lwd = 1,
       fontsize = 7,
       fs.test.effect.subgroup = 0,
       fs.hetstat = 7,
       just = "left",
       addrow = FALSE,
       leftcols = cols.forest,
       leftlabs = labs.forest,
       colgap = "0.01cm",
       col.by = "black",
       hetlab = "",
       calcwidth.subgroup = TRUE,
       print.tau = TRUE,
       print.tau.ci = TRUE,
       text.random.w = sapply(meta.prepost.daa$k.w, paste, "observations"),
       sort = days_after_first)

dev.off()
# TODO: clean up *all* unused data objects
rm(data.prepost, data.prepost.jags,prepost.jags.results, meta.prepost.daa, meta.prepost.other, meta.prepost.psych)
