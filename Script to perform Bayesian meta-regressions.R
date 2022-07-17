
###### Perform BAYESIAN pre-post analysis for Psychological Distress ------
# #Note that a different data making function and model need to be used, because I we have only one condition
# #use: makeprepostJAGSdata1C and prepostMHmodelRE1C
#to see results from the meta-regressions, just print the produced objects, named e.g. prepostmetaregAGE,prepostmetaregSEX etc

prepostmetaregJAGSdata=makeprepostmetaregJAGSdata1C(data=dataprepostPDB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=pop0_central_age)
prepostmetaregAGE=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE1C,
                       n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)

prepostmetaregJAGSdata=makeprepostmetaregJAGSdata1C(data=dataprepostPDB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=pop0_percent_female)
prepostmetaregSEX=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE1C,
                       n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)

prepostmetaregJAGSdata=makeprepostmetaregJAGSdata1C(data=dataprepostPDB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=gdpDiv10000minusmin)
prepostmetaregGDP=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE1C,
                       n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)

prepostmetaregJAGSdata=makeprepostmetaregJAGSdata1C(data=dataprepostPDB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=giniMinusmin)
prepostmetaregGINI=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE1C,
                        n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)

prepostmetaregJAGSdata=makeprepostmetaregJAGSdata1C(data=dataprepostPDB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=study_design01)
prepostmetaregDESIGN=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE1C,
                          n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)

prepostmetaregJAGSdata=makeprepostmetaregJAGSdata1C(data=dataprepostPDB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=RoBRepresdich)
prepostmetaregRoBPrepres=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE1C,
                              n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)

prepostmetaregJAGSdata=makeprepostmetaregJAGSdata1C(data=dataprepostPDB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=RoBInfoBiasdich)
prepostmetaregRoBInfo=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE1C,
                           n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata1C(data=dataprepostPDB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=RoBNonRespdich)
prepostmetaregRoBNonResp=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE1C,
                              n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)

prepostmetaregJAGSdata=makeprepostmetaregJAGSdata1C(data=dataprepostPDB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=USA)
prepostmetaregUSA=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE1C,
                       n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)

prepostmetaregJAGSdata=makeprepostmetaregJAGSdata1C(data=dataprepostPDB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=CHINA)
prepostmetaregCHINA=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE1C,
                         n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)

dataprepostPDB$sqrtN<-sqrt(dataprepostPDB$sample_size)-sqrt(min(dataprepostPDB$sample_size))
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata1C(data=dataprepostPDB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=sqrtN)
prepostmetaregsqrtN=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE1C,
                         n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)

dataprepostPDB$OldBmeasurement<-sqrt(-dataprepostPDB$exact_days_after_first)
dataprepostPDB$OldBmeasurement[is.na(dataprepostPDB$OldBmeasurement)]<-0

prepostmetaregJAGSdata=makeprepostmetaregJAGSdata1C(data=dataprepostPDB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=OldBmeasurement)
prepostmetaregOldBmeasuremen=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE1C,
                                  n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)

#days      
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata1C(data=dataprepostPDB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=days_after_first)
prepostmetaregDays=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE1C,
                        n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)
#stringency      
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata1C(data=dataprepostPDB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=stringency)
prepostmetaregStrin=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE1C,
                         n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)

#cases      
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata1C(data=dataprepostPDB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=logconfirmed_cumulative100k)
prepostmetaregLogCumCases=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE1C,
                               n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)
#deaths       
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata1C(data=dataprepostPDB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=logdeaths_cumulative100k)
prepostmetaregLogCumDeaths=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE1C,
                                n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)









###### Perform BAYESIAN pre-post meta-regression analysis for Anxiety, Depression  ------

#to see results from the meta-regressions, just print the produced objects, named e.g. prepostmetaregAGE,prepostmetaregSEX etc


#age
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata(data=dataprepostB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=pop0_central_age)
prepostmetaregAGE=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE,
                       n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)
#female %
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata(data=dataprepostB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=pop0_percent_female)
prepostmetaregSEX=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE,
                       n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)
#GDP      
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata(data=dataprepostB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=gdpDiv10000minusmin)
prepostmetaregGDP=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE,
                       n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)
#GINI      
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata(data=dataprepostB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=giniMinusmin)
prepostmetaregGINI=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE,
                        n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)
#Design      
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata(data=dataprepostB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=study_design01)
prepostmetaregDESIGN=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE,
                          n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)
#RoB representative sample     
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata(data=dataprepostB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=RoBRepresdich)
prepostmetaregRoBPrepres=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE,
                              n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)
#RoB information bias     
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata(data=dataprepostB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=RoBInfoStudy)
prepostmetaregRoBInfo=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE,
                           n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)
#RoB non-response bias 
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata(data=dataprepostB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=RoBNonRespdich)
prepostmetaregRoBNonResp=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE,
                              n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)
#USA vs other      
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata(data=dataprepostB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=USA)
prepostmetaregUSA=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE,
                       n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)
#China vs other
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata(data=dataprepostB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=CHINA)
prepostmetaregCHINA=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE,
                         n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)
#sample size      
dataprepostB$sqrtN<-sqrt(dataprepostB$sample_size)-sqrt(min(dataprepostB$sample_size))
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata(data=dataprepostB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=sqrtN)
prepostmetaregsqrtN=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE,
                         n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)
#baseline measurement timing      
dataprepostB$OldBmeasurement<-sqrt(-dataprepostB$exact_days_after_first)
dataprepostB$OldBmeasurement[is.na(dataprepostB$OldBmeasurement)]<-0

prepostmetaregJAGSdata=makeprepostmetaregJAGSdata(data=dataprepostB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=OldBmeasurement)
prepostmetaregOldBmeasuremen=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE,
                                  n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)

#days      
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata(data=dataprepostB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=days_after_first)
prepostmetaregDays=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE,
                        n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)
#stringency      
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata(data=dataprepostB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=stringency)
prepostmetaregStrin=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE,
                         n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)

#cases      
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata(data=dataprepostB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=logconfirmed_cumulative100k)
prepostmetaregLogCumCases=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE,
                               n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)
#deaths       
prepostmetaregJAGSdata=makeprepostmetaregJAGSdata(data=dataprepostB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6, covariate=logdeaths_cumulative100k)
prepostmetaregLogCumDeaths=jags(data=prepostmetaregJAGSdata,inits=NULL,parameters.to.save = c("muMH", "tauMH", "gamma"),model.file = prepostmetaregMHmodelRE,
                                n.chains=3,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 2)

rm(dataprepostPDB,dataprepostB)
