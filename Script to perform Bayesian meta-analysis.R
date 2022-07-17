###### Perform Frequentist pre-during analysis for all conditions ------
a=dataset %>% group_by(record_id,condition) %>% summarise(sum(is_prepandemic))
keep<-filter(a,`sum(is_prepandemic)`>0) %>% select(record_id,condition) 
dataprepost=(inner_join(dataset,keep,by=c("record_id", "condition")))
dataprepost$se=sqrt(dataprepost$var)
dataprepost$days_after_first<-as.integer(dataprepost$days_after_first)
dataprepostB=dataprepost[dataprepost$condition %in% c("Anxiety","Depression"),]
dataprepostPDB=dataprepost[dataprepost$condition %in% c("Psychological distress"),]
rm(a,keep)

#split the conditions in two groups because they are too many to plot
d1=dataprepost[dataprepost$condition %in% c("Anxiety","Depression"),]
d1$logss<-log(d1$meanStudySS)
d2=dataprepost[dataprepost$condition %in% c( "Psychological distress", "Mental Wellbeing","Life satisfaction","Sleep disturbance","Internet Gaming Disorder","Problematic smartphone-application use","Problematic social media use", "Somatoform disorder" )  ,]
d3=dataprepost[dataprepost$condition %in% c( "Psychological distress")  ,]

#meta-analysis for "Anxiety","Depression"
prepostmeta1=metagen(d,se,studlab=authoryear,subset = is_prepandemic==0 , subgroup=condition, data=d1, overall =T, fixed=F, prediction=T,prediction.subgroup =T)
#meta-analysis for other conditions
prepostmeta2=metagen(d,se,studlab=authoryear,subset = is_prepandemic==0, subgroup=condition, data=d2, overall =F)
#meta-analysis for psychological distress
prepostmeta3=metagen(d,se,studlab=authoryear,subset = is_prepandemic==0, subgroup=condition, data=d3, overall =F)

#funnel plot for anxiety and depression
funnel(prepostmeta1,contour = c(0.9, 0.95, 0.99),xlab="SMD anxiety and depression")
legend(-0.50, 0.02,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "< 0.01"),
       fill = c("darkgray", "gray", "lightgray"))
#funnel plot for psychological distress
funnel(prepostmeta3,contour = c(0.9, 0.95, 0.99),xlab="SMD psychological distress")
legend(-0.50, 0.02,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "< 0.01"),
       fill = c("darkgray", "gray", "lightgray"))

#Create Figure 2b
pdf(paste("Figure 2b Pre-post forest plot other conditions",".pdf",sep=""), width = 10, height = 10)
forest(prepostmeta2,
       test.overall=F, comb.fixed=F,print.byvar=T, overall =F,weight.study="same",
       title="Post-Pre SMD",pooled.totals=F, xlab=paste("SMD"),smlab="SMD",print.pval.Q=F,
       plotwidth="6cm",print.I2=F,lwd=1,fontsize=7,fs.test.effect.subgroup=0, subgroup=F,
       fs.hetstat=7, just="left",addrow=F,
       leftcols = c("studlab","country","population","sample_size","pop0_percent_female", "rob_info_bias.x","rob_non_bias.x","days_after_first"),
       leftlabs=c("Study","Country", "Population","Sample size","Fraction females","Information bias","Non-reponse bias","Days in pandemic"),
       colgap="0.01cm",col.by="black",hetlab="",calcwidth.subgroup=F,overall.hetstat=F,hetstat=F,
       text.random.w=sapply(prepostmeta2$k.w, paste, "observations"),sort=days_after_first)
dev.off()


###### Perform BAYESIAN pre-during analysis for  Anxiety, Depression ------
source("Make pre-post Bayesian meta-analysis functions.R")

prepostJAGSdata<-makeprepostJAGSdata(data=dataprepostB,studyid=record_id, y=d,var,covlabel="cov", corr=0.6)
prepostJAGSdata$Y1<-prepostJAGSdata$Y
prepostJAGSresults<-jags(data=prepostJAGSdata,inits=NULL,parameters.to.save = c("mu","muMH", "predSMD", "tau", "tauMH","predSMDMH"),model.file = prepostMHmodelRE,
                         n.chains=3,n.iter = 50000,n.burnin = 10000,DIC=F,n.thin = 3)

mus=as.data.frame(t(dplyr::select(as.data.frame(t(prepostJAGSresults$BUGSoutput$summary)),starts_with("mu"))))
taus=as.data.frame(t(dplyr::select(as.data.frame(t(prepostJAGSresults$BUGSoutput$summary)),starts_with("tau"))))
predSMDs=as.data.frame(t(dplyr::select(as.data.frame(t(prepostJAGSresults$BUGSoutput$summary)),starts_with("predSMD"))))
prepostresults=rbind(mus,predSMDs,taus) %>% select(`50%`,`2.5%`,`97.5%`) %>% rename("Posterior Median"=`50%`, "lowCI"=`2.5%`,"highCI" =`97.5%` )
rm(taus,mus,predSMDs,d1,d2)

###Trick forest.meta to show the  the Bayesian results and create Figure 2a----

prepostmeta1B<-prepostmeta1

prepostmeta1B$TE.random <- prepostresults[3,1]
prepostmeta1B$lower.random <- prepostresults[3,2]
prepostmeta1B$upper.random <- prepostresults[3,3]
prepostmeta1B$lower.predict <- prepostresults[6,2]
prepostmeta1B$upper.predict <- prepostresults[6,3]

prepostmeta1B$TE.random.w<- prepostresults[1:2,1]
prepostmeta1B$lower.random.w <- prepostresults[1:2,2]
prepostmeta1B$upper.random.w <- prepostresults[1:3,3]
prepostmeta1B$lower.predict.w <- prepostresults[4:5,2]
prepostmeta1B$upper.predict.w <- prepostresults[4:5,3]

prepostmeta1B$tau.w<-prepostresults[7:8,1]
prepostmeta1B$tau<-prepostresults[9,1]

prepostmeta1B$lower.tau.w<-prepostresults[7:8,2]
prepostmeta1B$upper.tau.w<-prepostresults[7:8,3]
prepostmeta1B$lower.tau<-prepostresults[9,2]
prepostmeta1B$upper.tau<-prepostresults[9,3]

pdf(paste("Figure 2a Pre-post forest plot Anxiety and Depression",".pdf",sep=""), width = 11, height = 9)
forest(prepostmeta1B,test.overall=F, prediction=T, comb.common=F,print.byvar=T, overall =T,
       title="Post-Pre SMD",pooled.total=F, xlab=paste("SMD"),smlab="SMD",print.pval.Q=F,
       plotwidth="6cm",print.I2=F,lwd=1,fontsize=7,fs.test.effect.subgroup=0,
       fs.hetstat=7, just="left",addrow=F,
       leftcols = c("studlab","country","population","sample_size","pop0_percent_female", "rob_info_bias.x","rob_non_bias.x","days_after_first"),
       leftlabs=c("Study","Country", "Population","Sample size","Fraction females","Information bias","Non-reponse bias","Days in pandemic"),
       colgap="0.01cm",col.by="black",hetlab="",calcwidth.subgroup=T,print.tau=T,print.tau.ci=T,
       text.random.w=sapply(prepostmeta1$k.w, paste, "observations"), sort=days_after_first)

dev.off()

rm(d3,dataprepost, prepostJAGSdata,prepostJAGSresults, prepostmeta1,prepostmeta1B,prepostmeta2,prepostmeta3)
