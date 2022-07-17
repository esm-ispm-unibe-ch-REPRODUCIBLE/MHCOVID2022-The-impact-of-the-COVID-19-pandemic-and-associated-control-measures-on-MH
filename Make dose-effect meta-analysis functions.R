

##                <<<<<<<<<<<<<<<< Pooled analysis for all conditions >>>>>>>>>>>>>>>>>>>>>>>>>>>>  ------

  #### 1 FUNCTION TO CREATE THE JAGS DATA---- 
makeJAGSDoseEffectDataPOOLED=function(data,studyid,dose, y,var,covlabel="cov", corr=0.6, knots=c(0.20,0.50,0.80), step.pred=1)
{#data=dataset
  data$studyid = eval(substitute(studyid), data) # data$studyid=data$record_id
  data$dose=eval(substitute(dose), data)#data$dose=data$days_after_first
  data$y = eval(substitute(y), data)#data$y=data$d
  data$var = eval(substitute(var), data)
  #!!!! put here the direction
  
  data=data[order(data$studyid),]
  data=data %>% filter(!is.na(var) & !is.na(dose) & !is.na(y))
  ns=length(unique(data$studyid))
  nc=length(unique(data$condition))
  nts <- as.numeric(table(data$studyid))#***combination of number of timepoints and conditions by study
  max.nts <- max(nts)
  maxcov=dim(combn(max.nts-1,2))[2]
  data$condid<-as.numeric(as.factor(data$condition))
  data$dose2<-as.numeric(unlist(rms::rcs(unlist(data$dose),knots)[,2]))
  data=data %>% group_by(record_id, condition) %>% group_modify(~ mutate(.x,tid=1:nrow(.x)))
  ntbycon=data %>% count(record_id, condid, condition)
  ncond=(ungroup(data %>% count(record_id, condition)) %>% count(record_id))$n #**number of reported conditions in each study
  data=data %>% group_by(record_id) %>% group_modify(~ add_count(.x)) %>% group_modify(~ add_row(.x, d=rep(NA, maxcov - nrow(.x))))
  
  #Create the outcome data Y matrix
  Y=data %>% select(d,tid,condid,condition) %>% filter(tid>1)
  dimY=max(table(Y$record_id))
  Y=ungroup(Y%>% group_by(record_id)%>% group_modify(~ add_row(.x, d=rep(NA, dimY - nrow(.x)))) )
  #condmat=matrix(Y$condid,ncol=dimY,byrow=T)#**
  Y=matrix(Y$d,ncol=dimY,byrow=T)#**
  
  #Create the two doses matrices
  dimdose=max(table(data$record_id))
  dose1=matrix(data$dose,ncol=dimdose,byrow=T)#*
  dose2=matrix(data$dose2,ncol=dimdose,byrow=T)#*
  
  #Create the condition by timepoints matrix
  dimntbycon=max(ncond)
  ntbycon=ungroup(ntbycon%>% group_by(record_id)%>% group_modify(~ add_row(.x, n=rep(NA, dimntbycon - nrow(.x)))) )
  condmat=matrix(ntbycon$condid,ncol=dimntbycon,byrow=T)#**
  ntbycon=matrix(ntbycon$n,ncol=dimntbycon,byrow=T)#**
  
  #create the precision matrices
  datasplit<-split(dplyr::select(as.data.frame(data),var,starts_with(covlabel),tid,condid),data$studyid )
  defaultW <- getOption("warn")  #temporarily suppress warnings
  options(warn = -1)
  
  create.covmat2=function(mat,corr=0.6){
    #corr is the correlation between SMDs from different conditions
    mat=mat[mat$tid>1,]
    dims<-dim(mat)[1]
    matsplit<-split(dplyr::select(mat,var,starts_with(covlabel)),mat$condid )
    create.covmat1=function(mat){
      dims<-dim(mat)[1]
      newmat<-matrix(0,nrow=dims, ncol=dims)
      newmat[lower.tri(newmat)]=unlist(mat[1,-1])
      newmat<-newmat+t(newmat)
      diag(newmat)<-mat$var
      newmat}
    percondVlist=bdiag(lapply(matsplit,create.covmat1))
    Covmati=sqrt(bdiag(diag(percondVlist))%*%diag(percondVlist))*corr*(percondVlist==0)+percondVlist
    Covmati
    # print(j)
    # j=j+1
    # mat=datasplit[[j]]
   }#makes the variance-covariances in a study with multiple timepoints and conditions. Assumes correlation 0.6 between condition SMDs
  
  Vlist=lapply(datasplit,create.covmat2)
  options(warn = defaultW)#bring warnings back
  
  prec<-lapply(Vlist,solve)
  prec<-lapply(prec,as.matrix)
  precmat<-plyr::ldply(prec,rbind)
  
  #create data needed for predictions
  new.dose=seq(round(min(data$dose,na.rm=T)),round(max(data$dose,na.rm=T)),by=step.pred)
  f.new.dose=as.vector(unlist(rms::rcs(new.dose,knots)[,2]))
  nt.new=length(new.dose)
  
  #create  position matrices for the outcomes Y and the mean (k1 and k2)
  
  k2=t(apply(ntbycon-1,1,cumsum))
  # print(k2)
  u=cbind(rep(1,ns),ntbycon-1)
  u<-t(apply(u, 1, cumsum))
  # print((u))
  k1=u[,1:dim(k2)[2]]
  k1[is.na(k2)]<-NA#**
  # cat(paste("k1="))
  # print(k1)
  
  #create  position matrices for the dose (k3 and k4)
  k3<-k1+col(k1)-1#**
  # print(paste("k3="))
  # print(k3)
  k4=t(apply(ntbycon,1,cumsum))#**
  
  JAGSdata <- list(Y = Y, dose=dose1, dose2=dose2,nts = nts, ncond=ncond, condmat=condmat, k1=k1,k2=k2,k3=k3,k4=k4, nc=nc, ns = ns, prec = precmat[,-1],nt.new=nt.new,f.new.dose=f.new.dose,new.dose=new.dose)
  return(JAGSdata)
}

    #### 2  JAGS MODEL random effects pooling conditions---------------
   pooledMHmodelRE<-function(){
  b[1] <-1
  
  for (i in 1:ns) { ## for each study
    # multivariate normal likelihood of the response 
    Y[i,1:(nts[i]-ncond[i])]  ~ dmnorm(mean[i,1:(nts[i]-ncond[i])], prec[(b[i]):(b[i]+nts[i]-ncond[i]-1),1:(nts[i]-ncond[i])])
    # window to change the index in precision matrix
    b[i+1] <- b[i]+ nts[i]-ncond[i]
    
    for (s in 1:ncond[i]) {
      mean[i,k1[i,s]:k2[i,s]]<-b1[i,condmat[i,s]]*(dose[i,(1+k3[i,s]):k4[i,s]]-dose[i, k3[i,s]])+ b2[i,condmat[i,s]]*(dose2[i,(1+k3[i,s]):k4[i,s]]-dose2[i, k3[i,s]])
    }
  }
  
  for(i in 1:ns){
    for(c in 1:nc){
      b1[i,c]~dnorm(beta1[c],prec.beta)
      b2[i,c]~dnorm(beta2[c],prec.beta)}}
  
  for(c in 1:nc){
    # beta1[c]<-beta1.pooled
    # beta2[c]<-beta2.pooled
    beta1[c]~dnorm(beta1.pooled,prec.betac)
    beta2[c]~dnorm(beta2.pooled,prec.betac)
  }
  
  beta1.pooled ~ dnorm(0,10)
  beta2.pooled ~ dnorm(0,10)
  
  prec.beta<-1/variance
  variance<-tau*tau
  tau~ dnorm(0,1)%_%T(0,2)
  prec.betac<-1/variancec
  variancec<-tauc*tauc
  tauc~ dnorm(0,1)%_%T(0,2)
  
  for( j in 1:nt.new){
    for(c in 1:nc){
      predSMD[j,c] <- (beta1[c]*new.dose[j]+ beta2[c]*f.new.dose[j])
    }}
  
  for( j in 1:nt.new){
    predSMDpooled[j] ~ dnorm(beta1.pooled*new.dose[j]+ beta2.pooled*f.new.dose[j],0.5*prec.betac)
    SMDpooled<- (beta1.pooled*new.dose[j]+ beta2.pooled*f.new.dose[j])
  }
  
}

    #### 2  JAGS MODEL fixed effects pooling conditions---------------
    pooledMHmodelFE<-function(){
  b[1] <-1
  
  for (i in 1:ns) { ## for each study
    # multivariate normal likelihood of the response 
    Y[i,1:(nts[i]-ncond[i])]  ~ dmnorm(mean[i,1:(nts[i]-ncond[i])], prec[(b[i]):(b[i]+nts[i]-ncond[i]-1),1:(nts[i]-ncond[i])])
    # window to change the index in precision matrix
    b[i+1] <- b[i]+ nts[i]-ncond[i]
    
    for (s in 1:ncond[i]) {
      mean[i,k1[i,s]:k2[i,s]]<-b1[i,condmat[i,s]]*(dose[i,(1+k3[i,s]):k4[i,s]]-dose[i, k3[i,s]])+ b2[i,condmat[i,s]]*(dose2[i,(1+k3[i,s]):k4[i,s]]-dose2[i, k3[i,s]])
    }
  }
  
  for(i in 1:ns){
    for(c in 1:nc){
      b1[i,c]~dnorm(beta1[c],prec.beta)
      b2[i,c]~dnorm(beta2[c],prec.beta.f)}}
  
  for(c in 1:nc){
    beta1[c]<-beta1.pooled
    beta2[c]<-beta2.pooled
  }
  
  beta1.pooled ~ dnorm(0,1)
  beta2.pooled ~ dnorm(0,1)
  
  prec.beta<-1/variance
  variance<-tau*tau
  tau~ dnorm(0,1)%_%T(0,2)

  prec.beta.f<-1/variance.f
  variance.f<-tau.f*tau.f
  tau.f~ dnorm(0,1)%_%T(0,2)
  
  coef1~dnorm(beta1.pooled,prec.beta)
  coef2~dnorm(beta2.pooled,prec.beta.f)
  for(j in 1:nt.new){
      meanSMDpooled[j]<-beta1.pooled*new.dose[j]+ beta2.pooled*f.new.dose[j]
      #predSMDpooled[j] ~ dnorm(meanSMDpooled[j],0.5*prec.beta)
    predSMDpooled[j]<-coef1*new.dose[j]+ coef2*f.new.dose[j]
    }

}

    #### 2  JAGS MODEL fixed effects for different condition---------------
    pooledMHmodelFEpercond<-function(){
      b[1] <-1
      
      for (i in 1:ns) { ## for each study
        # multivariate normal likelihood of the response 
        Y[i,1:(nts[i]-ncond[i])]  ~ dmnorm(mean[i,1:(nts[i]-ncond[i])], prec[(b[i]):(b[i]+nts[i]-ncond[i]-1),1:(nts[i]-ncond[i])])
        # window to change the index in precision matrix
        b[i+1] <- b[i]+ nts[i]-ncond[i]
        
        for (s in 1:ncond[i]) {
          mean[i,k1[i,s]:k2[i,s]]<-b1[i,condmat[i,s]]*(dose[i,(1+k3[i,s]):k4[i,s]]-dose[i, k3[i,s]])+ b2[i,condmat[i,s]]*(dose2[i,(1+k3[i,s]):k4[i,s]]-dose2[i, k3[i,s]])
        }
      }
      
      for(i in 1:ns){
        for(c in 1:nc){
          b1[i,c]~dnorm(beta1[c],prec.beta)
          b2[i,c]~dnorm(beta2[c],prec.beta)}}
      
      for(c in 1:nc){
        beta1[c] ~ dnorm(0,1)
        beta2[c] ~ dnorm(0,1)
      }
      
      
      prec.beta<-1/variance
      variance<-tau*tau
      tau~ dnorm(0,1)%_%T(0,2)
      
      for( j in 1:nt.new){
        for(c in 1:nc){
          predSMD[j,c] <- (beta1[c]*new.dose[j]+ beta2[c]*f.new.dose[j])
        }}
      
      
    }
    
    #### 2  JAGS MODEL fixed effects linear only---------------
    pooledMHmodelFElinear<-function(){
      b[1] <-1
      
      for (i in 1:ns) { ## for each study
        # multivariate normal likelihood of the response 
        Y[i,1:(nts[i]-ncond[i])]  ~ dmnorm(mean[i,1:(nts[i]-ncond[i])], prec[(b[i]):(b[i]+nts[i]-ncond[i]-1),1:(nts[i]-ncond[i])])
        # window to change the index in precision matrix
        b[i+1] <- b[i]+ nts[i]-ncond[i]
        
        for (s in 1:ncond[i]) {
          mean[i,k1[i,s]:k2[i,s]]<-b1[i,condmat[i,s]]*(dose[i,(1+k3[i,s]):k4[i,s]]-dose[i, k3[i,s]])
        }
      }
      
      for(i in 1:ns){
        for(c in 1:nc){
          b1[i,c]~dnorm(beta1[c],prec.beta)
          }}
      
      for(c in 1:nc){
        beta1[c]<-beta1.pooled
      }
      
      beta1.pooled ~ dnorm(0,1)
    
      prec.beta<-1/variance
      variance<-tau*tau
      tau~ dnorm(0,1)%_%T(0,2)
      
      for( j in 1:nt.new){
        for(c in 1:nc){
          predSMD[j,c] <- (beta1[c]*new.dose[j]+ beta2[c]*f.new.dose[j])
        }}
      
      
    }
    
    
    #### 3  Function to plot results ---------------
plotMHall.fun<-function(allJAGSdata=allJAGSdata, allJAGSresults=allJAGSresults, ylimits, names.conditions=names.conditions, graphsfile="trend and traceplots", textfile="beta and tau", name.dose="dosename", append=T, traceplot=F,percondition=F){

  ### TABLES
  #betas and taus
  sink(file = paste(textfile,".txt",sep=""), append = append)
  betas=as.data.frame(t(dplyr::select(as.data.frame(t(allJAGSresults$BUGSoutput$summary)),starts_with("beta"))))
  taus=as.data.frame(t(dplyr::select(as.data.frame(t(allJAGSresults$BUGSoutput$summary)),starts_with("tau"))))
  print(kable(betas,  format = "markdown", digits = 3) )
  print(kable(taus,  format = "markdown", digits = 3) )
  sink()
  
  ### PLOTS
  pdf(paste(graphsfile,".pdf",sep=""), width = 10, height = 8)
  #plot the dose-effect lines
  
  predictedSMD=as.data.frame(t(dplyr::select(as.data.frame(t(allJAGSresults$BUGSoutput$summary)),starts_with("predSMD"))))
  predictedSMD=dplyr::select(predictedSMD,mean,`2.5%`,`97.5%`)
  meanSMD=as.data.frame(t(dplyr::select(as.data.frame(t(allJAGSresults$BUGSoutput$summary)),starts_with("meanSMD"))))
  meanSMD=dplyr::select(meanSMD,mean,`2.5%`,`97.5%`)
  predSMD<-cbind.data.frame(predictedSMD,meanSMD)
  predSMD$new.dose=allJAGSdata$new.dose
  names(predSMD)<-c("Pred","2.5%PrI","97.5%PrI", "Mean","2.5%CrI","97.5%CrI","new.dose")
  
  
  par(mfrow=c(1,1))
  if(percondition==T){predSMD$condition=c(sapply(c(names.conditions),rep,length(allJAGSdata$new.dose)))
  plot(allJAGSdata$new.dose,rep(0,length(allJAGSdata$new.dose)), type='n', xlab=name.dose,ylab=paste(unique(predSMD$condition), "SMD"),ylim=ylimits)}
  if(percondition==F){predSMD$condition=c(sapply(c("all"),rep,length(allJAGSdata$new.dose)))
  predsmd<-filter(predSMD,condition=="all")
  plot(predsmd$new.dose,predsmd$Mean,type='l', xlab=name.dose,ylab= "SMD",ylim=ylimits,col="red",lwd=2)
  lines(predsmd$new.dose,predsmd$`2.5%PrI`,lty=3,col="blue")
  lines(predsmd$new.dose,predsmd$`97.5%PrI`,lty=3,col="blue")
  lines(predsmd$new.dose,predsmd$`2.5%CrI`,lty=2)
  lines(predsmd$new.dose,predsmd$`97.5%CrI`,lty=2)
  }
  rug(allJAGSdata$dose)
  j=1
  for(i in c(names.conditions)){
    predsmd<-filter(predSMD,condition==i)
    lines(predsmd$new.dose,predsmd$mean,type='l', xlab=name.dose,ylab=paste(unique(predsmd$condition), "SMD"),ylim=ylimits,col=j,lwd=1)
    j<-j+1}
  
  #traceplots
  if(traceplot==T){
    traceplot(allJAGSresults,var="beta1.pooled")
    traceplot(allJAGSresults,var="beta2.pooled")
    traceplot(allJAGSresults,var="tau")
    }
  dev.off()
  
  out=predSMD
  return(out)
}

    #### 4  Wrap function to fit pooled dose-effects and show results ---------------
# name.dose="Facial_cover"
# dataset=dataset
# model=pooledMHmodelFE
# n.iter=10000
# corr=0.6
# graphsfile="days_after_first"
# textfile="beta and tau days_after_first FE"
# step.pred=1
# percondition=F
#quant.rcs=c(0.2,0.6,0.8)


doseeffectPOOLEDJAGSresults.fun<-function(name.dose,dataset=dataset,model=pooledMHmodelFE,n.iter=1000,corr=0.6, graphsfile, textfile,updateJAGS=F,allJAGSresults, knots,ylimits=c(-1,0.5),traceplot=F, percondition=F,step.pred=1,quant.rcs=c(0.2,0.5,0.8)){
data=select(dataset, record_id,d,var,condition,starts_with("cov"), starts_with(name.dose) )  
  #data=filter(data, condition=="Depression"|condition=="Anxiety" )
  data$dose=as.numeric(unlist(dplyr::select(data,starts_with(name.dose))))
  data=filter(data, !is.na(dose) )  
  names.conditions=sort(unique(data$condition))
  if(missing(knots)) {myknots=unlist(quantile(unlist(data$dose),quant.rcs))}
  #myknots=c(5,40,60)
  
  sink(file = paste(textfile,".txt",sep=""), append = F)
  cat(paste(" \n ____ Analysis of", length(unique(data$record_id)), "studies in",length(names.conditions), " conditions ____ \n")[1])
  cat(paste("\n____ Knots in",name.dose , " ____\n"))
  print(myknots)
  print(kable(margin.table(table(data$condition,data$record_id),1),  format = "markdown", digits = 3, caption = "Available timepoints per condition",col.names = c("Condition", "Nr of Timepoints") ))
  sink()
  
  allJAGSdata=makeJAGSDoseEffectDataPOOLED(data,studyid=record_id,dose=data$dose, y=d,var=var,covlabel="cov", corr=corr, knots=myknots,step.pred=step.pred)
  if(updateJAGS==F){   allJAGSresults=jags(data = allJAGSdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled',"predSMDpooled","meanSMDpooled", "tau","tau.f"),model.file = model,
                      n.chains=3,n.iter = 1000,n.burnin = 100,DIC=F,n.thin = 2)
}
  allJAGSresults=update(allJAGSresults, n.iter=n.iter, n.burnin = round(n.iter/5))
  plotdata=plotMHall.fun(allJAGSdata=allJAGSdata, allJAGSresults=allJAGSresults, ylim=ylimits,names.conditions=names.conditions,  graphsfile=graphsfile, textfile=textfile, name.dose=name.dose, append=T,traceplot=traceplot,percondition=percondition)
  out<-list(allJAGSresults=allJAGSresults,allJAGSdata=allJAGSdata,plotdata=plotdata)
  out
  }


##

