#-------------------------------------------------------------------------------
# File  : 04_run_dose_response_models.R
#
# Perform dose-response anlyses for Depression and Anxiety outcomes.
#-------------------------------------------------------------------------------
source("Make dose-effect meta-analysis functions.R")
set.seed(100371)
data.daa <- filter(data, (condition == "Depression" | condition == "Anxiety"))

# Dose-response analysis: days since first recorded case -----------------------
DaysFEAD <- doseeffectPOOLEDJAGSresults.fun(
  name.dose = "days_after_first",
  dataset = data.daa,
  model = pooledMHmodelFE,
  n.iter = 20000,
  corr = 0.6,
  graphsfile = "Days since first case FE",
  textfile = "beta and tau days since first case FE")

DaysToPlot <- as.data.frame(DaysFEAD$plotdata)
write.csv(DaysToPlot,"DaysToPlot2.csv")

# Dose-response analysis: stringency -------------------------------------------
StrFEAD<-doseeffectPOOLEDJAGSresults.fun(
  name.dose = "stringency" ,
  dataset = data.daa,
  model = pooledMHmodelFE,
  n.iter = 20000,
  corr = 0.6,
  graphsfile = "Stringency FE",
  textfile = "beta and tau stringency FE")

StringencyToPlot<-as.data.frame(StrFEAD$plotdata)
write.csv(StringencyToPlot,"StringencyToPlot.csv")

# Dose-response analysis: cumulative cases per 100K (log-transformed) ----------
logCasesFEAD<-doseeffectPOOLEDJAGSresults.fun(
  name.dose = "logconfirmed_cumulative100k",
  dataset = data.daa,
  model = pooledMHmodelFE,
  n.iter = 20000,
  corr = 0.6,
  graphsfile = "log(ccases) FE",
  textfile = "beta and tau log(ccases) FE",
  step.pred = 0.1,
  ylimits = c(-1,0.5))

logCasesToPlot<-as.data.frame(logCasesFEAD$plotdata)
write.csv(logCasesToPlot,"logCasesToPlot.csv")

# Dose-response analysis: cumulative deaths per 100K (log-transformed) ---------
logDeathsFEAD<-doseeffectPOOLEDJAGSresults.fun(
  name.dose = "logdeaths_cumulative100k",
  dataset = data.daa,
  model = pooledMHmodelFE,
  n.iter = 20000,
  corr = 0.6,
  graphsfile = "log(cdeaths) FE",
  textfile = "beta and tau log(cdeaths) FE",
  step.pred = 0.1,
  ylimits = c(-1,0.5))

logDeathsToPlot<-as.data.frame(logDeathsFEAD$plotdata)
write.csv(logDeathsToPlot,"logDeathsToPlot.csv")
