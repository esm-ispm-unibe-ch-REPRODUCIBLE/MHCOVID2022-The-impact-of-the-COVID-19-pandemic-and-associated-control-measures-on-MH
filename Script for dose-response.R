
datasetOriginal=dataset #copy, now datasetOriginal has all conditions, and dataset only depression and anxiety

source("Make dose-effect meta-analysis functions.R")
dataset=filter(datasetOriginal, condition=="Depression"|condition=="Anxiety")

set.seed(100371)
DaysFEAD<-doseeffectPOOLEDJAGSresults.fun(
  name.dose="days_after_first",
  dataset=dataset,
  model=pooledMHmodelFE,
  n.iter=20000,
  corr=0.6, 
  graphsfile="Days since first case FE",
  textfile="beta and tau days since first case FE")
#DaysFEAD<-doseeffectPOOLEDJAGSresults.fun(name.dose="days_after_first",dataset=dataset,model=pooledMHmodelFE,n.iter=80000,corr=0.6, graphsfile="Days since first case FE", textfile="beta and tau days since first case FE",allJAGSresults=DaysFEAD$allJAGSresults,updateJAGS=T)
DaysToPlot<-as.data.frame(DaysFEAD$plotdata) 
write.csv(DaysToPlot,"DaysToPlot2.csv")

StrFEAD<-doseeffectPOOLEDJAGSresults.fun(
  name.dose="stringency" ,
  dataset=dataset,
  model=pooledMHmodelFE,
  n.iter=20000,
  corr=0.6,
  graphsfile="Stringency FE",
  textfile="beta and tau stringency FE")
#StrFEAD<-doseeffectPOOLEDJAGSresults.fun(name.dose="stringency" ,dataset=dataset,model=pooledMHmodelFE,n.iter=80000,corr=0.6, graphsfile="Stringency FE", textfile="beta and tau stringency FE",allJAGSresults=StrFEAD$allJAGSresults,updateJAGS=T)
StringencyToPlot<-as.data.frame(StrFEAD$plotdata) 
write.csv(StringencyToPlot,"StringencyToPlot.csv")

logCasesFEAD<-doseeffectPOOLEDJAGSresults.fun(
  name.dose="logconfirmed_cumulative100k",
  dataset=dataset,
  model=pooledMHmodelFE,
  n.iter=20000,
  corr=0.6, 
  graphsfile="log(ccases) FE", 
  textfile="beta and tau log(ccases) FE",
  step.pred=0.1,
  ylimits=c(-1,0.5))
#logCasesFEAD<-doseeffectPOOLEDJAGSresults.fun(name.dose="logconfirmed_cumulative100k" ,dataset=dataset,model=pooledMHmodelFE,n.iter=80000,corr=0.6, graphsfile="log(ccases) FE", textfile="beta and tau log(ccases) FE",allJAGSresults=logCasesFEAD$allJAGSresults,updateJAGS=T,step.pred=0.1,ylimits=c(-1,0.5))
logCasesToPlot<-as.data.frame(logCasesFEAD$plotdata) 
write.csv(logCasesToPlot,"logCasesToPlot.csv")

logDeathsFEAD<-doseeffectPOOLEDJAGSresults.fun(
  name.dose="logdeaths_cumulative100k",
  dataset=dataset,
  model=pooledMHmodelFE,
  n.iter=20000,
  corr=0.6,
  graphsfile="log(cdeaths) FE",
  textfile="beta and tau log(cdeaths) FE",
  step.pred=0.1,
  ylimits=c(-1,0.5) )
# logDeathsFEAD<-doseeffectPOOLEDJAGSresults.fun(name.dose="logdeaths_cumulative100k",dataset=dataset,model=pooledMHmodelFE,n.iter=80000,corr=0.6, graphsfile="log(cdeaths) FE", textfile="beta and tau log(cdeaths) FE",step.pred=0.1,ylimits=c(-1,0.5),updateJAGS=T,allJAGSresults=logDeathsFEAD$allJAGSresults)
logDeathsToPlot<-as.data.frame(logDeathsFEAD$plotdata)
write.csv(logDeathsToPlot,"logDeathsToPlot.csv")

