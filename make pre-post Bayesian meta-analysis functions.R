#***** META-ANALYSES for more than one condition - random effects
#### 1 FUNCTION TO CREATE THE JAGS DATA---- 
makeprepostJAGSdata=function(data,studyid, y,var,covlabel="cov", corr=0.6)
{
  data$studyid = eval(substitute(studyid), data) 
  data$y = eval(substitute(y), data)
  data$var = eval(substitute(var), data)
  
  # data=dataprepostB
  # data$studyid=data$record_id
  # data$y=data$d
  # covlabel="cov"

  data=data[order(data$studyid),]
  data=data %>% filter(!is.na(var)  & !is.na(y))
  ns=length(unique(data$studyid))
  nc=length(unique(data$condition))
  nts <- as.numeric(table(data$studyid))#***combination of number of timepoints and conditions by study
  max.nts <- max(nts)
  maxcov=dim(combn(max.nts-1,2))[2]
  data$condid<-as.numeric(as.factor(data$condition))
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
  
  
  #create  position matrices for the outcomes Y and the mean (k1 and k2)
  
  k2=t(apply(ntbycon-1,1,cumsum))
  u=cbind(rep(1,ns),ntbycon-1)
  u<-t(apply(u, 1, cumsum))
  k1=u[,1:dim(k2)[2]]
  k1[is.na(k2)]<-NA#**
  k3<-k1+col(k1)-1#**
  k4=t(apply(ntbycon,1,cumsum))#**
  #Create the two doses matrices
  dimdose=max(table(data$record_id))
  dose=matrix(rep(1, length(data$record_id)),ncol=dimdose,byrow=T)#*

  ### Print
  ConditionsKey<-cbind(1:nc,sort(unique(data$condition)))
  print(kable( ConditionsKey,  format = "markdown", digits = 3,col.names = c("Condition code", "Condition")) )
  
  prepostJAGSdata <- list(Y = Y, nts = nts, ncond=ncond, condmat=condmat, k1=k1,k2=k2, dose=dose, k3=k3,k4=k4,nc=nc, ns = ns, prec = precmat[,-1])
  return( prepostJAGSdata)
}

#### 2  JAGS MODEL random effects for conditions---------------
prepostMHmodelRE<-function(){
  # This part is for estimating results from each condition separately
  b[1] <-1
  for (i in 1:ns) { ## for each study
    # multivariate normal likelihood of the response 
    Y[i,1:(nts[i]-ncond[i])]  ~ dmnorm(mean[i,1:(nts[i]-ncond[i])], prec[(b[i]):(b[i]+nts[i]-ncond[i]-1),1:(nts[i]-ncond[i])])
    # window to change the index in precision matrix
    b[i+1] <- b[i]+ nts[i]-ncond[i]
    
    for (s in 1:ncond[i]) {
      mean[i,k1[i,s]:k2[i,s]]<-delta[i,condmat[i,s]]*dose[i,(1+k3[i,s]):k4[i,s]]
    }
  }
  
  # This part is for pooling the conditions together
  b1[1] <-1
  for (i in 1:ns) { ## for each study
    # multivariate normal likelihood of the response 
    Y1[i,1:(nts[i]-ncond[i])]  ~ dmnorm(meanMH[i,1:(nts[i]-ncond[i])], prec[(b1[i]):(b1[i]+nts[i]-ncond[i]-1),1:(nts[i]-ncond[i])])
    # window to change the index in precision matrix
    b1[i+1] <- b1[i]+ nts[i]-ncond[i]
    
    for (s in 1:ncond[i]) {
      meanMH[i,k1[i,s]:k2[i,s]]<-deltaMH[i,condmat[i,s]]*dose[i,(1+k3[i,s]):k4[i,s]]
    }
  }
  
  for(i in 1:ns){
    for(c in 1:nc){
      delta[i,c]~dnorm(mu[c],precc[c])
    }}
  
  for(i in 1:ns){
    for(c in 1:nc){
      deltaMH[i,c]~dnorm(muMH,precMH)
    }}
  
  #Priors for the two conditions
    for(c in 1:nc){
      mu[c] ~ dnorm(0,0.001)
      precc[c]<-1/variance[c]
      variance[c]<-tau[c]*tau[c]
      tau[c]~dnorm(0,1)%_%T(0,2)
    }
  
  #Priors for the pooled estimate
  muMH ~ dnorm(0,0.00001)
  precMH<-1/varianceMH
  varianceMH<-tauMH*tauMH
  tauMH~dnorm(0,1)%_%T(0,2)
  
  
 #prediction intervals
    for(c in 1:nc){
      predSMD[c]~dnorm(mu[c],precc[c])
    }
  predSMDMH~dnorm(muMH,precMH)
}


##### Meta-regressions 

#***** META-REGRESSIONS
#### 1 FUNCTION TO CREATE THE META-REGRESSION JAGS DATA---- 
makeprepostmetaregJAGSdata=function(data,studyid, y,var,covlabel="cov", corr=0.6, covariate)
{
  data$studyid = eval(substitute(studyid), data) 
  data$y = eval(substitute(y), data)
  data$var = eval(substitute(var), data)
  data$covariate = eval(substitute(covariate), data)
  
  # data=dataprepostB
  # data$studyid=data$record_id
  # data$y=data$d
  # covlabel="cov"
  
  data=data[order(data$studyid),]
  data=data %>% filter(!is.na(var)  & !is.na(y) & !is.na(covariate))
  ns=length(unique(data$studyid))
  nc=length(unique(data$condition))
  nts <- as.numeric(table(data$studyid))#***combination of number of timepoints and conditions by study
  max.nts <- max(nts)
  maxcov=dim(combn(max.nts-1,2))[2]
  data$condid<-as.numeric(as.factor(data$condition))
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
  
  #Create the covariate matrix
  dimdose=max(table(data$record_id))
  covariate=matrix(data$covariate,ncol=dimdose,byrow=T)#*
  
  
  #create  position matrices for the outcomes Y and the mean (k1 and k2)
  
  k2=t(apply(ntbycon-1,1,cumsum))
  u=cbind(rep(1,ns),ntbycon-1)
  u<-t(apply(u, 1, cumsum))
  k1=u[,1:dim(k2)[2]]
  k1[is.na(k2)]<-NA#**
  k3<-k1+col(k1)-1#**
  k4=t(apply(ntbycon,1,cumsum))#**
  #Create the two doses matrices
  dimdose=max(table(data$record_id))
  dose=matrix(rep(1, length(data$record_id)),ncol=dimdose,byrow=T)#*
  
  ### Print
  ConditionsKey<-cbind(1:nc,sort(unique(data$condition)))
  print(kable( ConditionsKey,  format = "markdown", digits = 3,col.names = c("Condition code", "Condition")) )
  
  prepostJAGSdata <- list(Y = Y, nts = nts, ncond=ncond, condmat=condmat, k1=k1,k2=k2, dose=dose, k3=k3,k4=k4,nc=nc, ns = ns, prec = precmat[,-1], covariate=covariate)
  return( prepostJAGSdata)
}

#### 2  JAGS  META-REGRESSION  MODEL random effects for conditions---------------
prepostmetaregMHmodelRE<-function(){
  #this model does not assume any difference between the conditions
  b[1] <-1
  for (i in 1:ns) { ## for each study
    # multivariate normal likelihood of the response 
    Y[i,1:(nts[i]-ncond[i])]  ~ dmnorm(mean[i,1:(nts[i]-ncond[i])], prec[(b[i]):(b[i]+nts[i]-ncond[i]-1),1:(nts[i]-ncond[i])])
    # window to change the index in precision matrix
    b[i+1] <- b[i]+ nts[i]-ncond[i]
    
    for (s in 1:ncond[i]) {
      mean[i,k1[i,s]:k2[i,s]]<-delta[i,condmat[i,s]]*dose[i,(1+k3[i,s]):k4[i,s]]+gamma*covariate[i,k1[i,s]:k2[i,s]]
    }
  }
  
  for(i in 1:ns){
    for(c in 1:nc){
      delta[i,c]~dnorm(muMH,precMH)
    }}

    #Priors
  muMH ~ dnorm(0,0.00001)
  precMH<-1/varianceMH
  varianceMH<-tauMH*tauMH
  tauMH~dnorm(0,1)%_%T(0,2)
  gamma~dnorm(0,0.00001)
  
  
  #prediction intervals
    predSMD~dnorm(muMH,precMH)

 
}


#***** META-ANALYSES for one condition - random effects- I use it for the Psychological distress analysis
#### 1 FUNCTION TO CREATE THE JAGS DATA---- 
makeprepostJAGSdata1C=function(data,studyid, y,var,covlabel="cov", corr=0.6)
{
  data$studyid = eval(substitute(studyid), data) 
  data$y = eval(substitute(y), data)
  data$var = eval(substitute(var), data)
  
  # data=dataprepostPDB
  # data$studyid=data$record_id
  # data$y=data$d
  # covlabel="cov"
  
  data=data[order(data$studyid),]
  data=data %>% filter(!is.na(var)  & !is.na(y))
  ns=length(unique(data$studyid))
  nc=length(unique(data$condition))
  nts <- as.numeric(table(data$studyid))#***combination of number of timepoints and conditions by study
  max.nts <- max(nts)
  maxcov=dim(combn(max.nts-1,2))[2]
  data$condid<-as.numeric(as.factor(data$condition))
  data=data %>% group_by(record_id, condition) %>% group_modify(~ mutate(.x,tid=1:nrow(.x)))
  ntbycon=data %>% count(record_id, condid, condition)
  ncond=(ungroup(data %>% count(record_id, condition)) %>% count(record_id))$n #**number of reported conditions in each study
  data=data %>% group_by(record_id) %>% group_modify(~ add_count(.x)) #%>% group_modify(~ add_row(.x, d=rep(NA, maxcov - nrow(.x))))
  
  #Create the outcome data Y matrix
  Y=data %>% select(d,tid,condid,condition) %>% filter(tid>1)
  dimY=max(table(Y$record_id))
  Y=ungroup(Y%>% group_by(record_id)%>% group_modify(~ add_row(.x, d=rep(NA, dimY - nrow(.x)))) )
  #condmat=matrix(Y$condid,ncol=dimY,byrow=T)#**
  Y=matrix(Y$d,ncol=dimY,byrow=T)#**
  
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
  }#makes the variance-covariances in a study with multiple timepoints and conditions. Assumes correlation 0.6 between condition SMDs
  
  Vlist=lapply(datasplit,create.covmat2)
  options(warn = defaultW)#bring warnings back
  
  prec<-lapply(Vlist,solve)
  prec<-lapply(prec,as.matrix)
  precmat<-plyr::ldply(prec,rbind)

  #Create a pseudo-dose matrix to change the dimensions as needed
  dimdose=max(table(data$record_id))
  dose=matrix(1,ncol=dimdose,nrow=ns)#*
  
  prepostJAGSdata <- list(Y = Y, nts = nts, dose=dose, ns = ns, prec = precmat[,-1])
  return( prepostJAGSdata)
}

#### 2  JAGS MODEL random effects for one condition---------------
prepostMHmodelRE1C<-function(){
    b1[1] <-1
    for (i in 1:ns) { ## for each study
      # multivariate normal likelihood of the response 
      Y[i,1:(nts[i]-1)]  ~ dmnorm(meanMH[i,1:(nts[i]-1)], prec[(b1[i]):(b1[i]+nts[i]-2),1:(nts[i]-1)])
      # window to change the index in precision matrix
      b1[i+1] <- b1[i]+ nts[i]-1
      meanMH[i,1:(nts[i]-1)]<-deltaMH[i]*dose[i,1:(nts[i]-1)]
    }
    
    for(i in 1:ns){
      deltaMH[i]~dnorm(muMH,precMH)
    }
    
    #Priors for the pooled estimate
    muMH ~ dnorm(0,0.00001)
    precMH<-1/varianceMH
    varianceMH<-tauMH*tauMH
    tauMH~dnorm(0,1)%_%T(0,2)
    
    #prediction intervals
    predSMDMH~dnorm(muMH,precMH)
  }



##### Meta-regressions 

#***** META-REGRESSIONS
#### 1 FUNCTION TO CREATE THE META-REGRESSION JAGS DATA---- 
makeprepostmetaregJAGSdata1C=function(data,studyid, y,var,covlabel="cov", corr=0.6, covariate)
{ 
  # data=dataprepostPDB
  # data$studyid=data$record_id
  # data$y=data$d
  # covlabel="cov"
  # data$covariate=data$pop0_central_age

  
  data$studyid = eval(substitute(studyid), data) 
  data$y = eval(substitute(y), data)
  data$var = eval(substitute(var), data)
  data$covariate = eval(substitute(covariate), data)

    data=data[order(data$studyid),]
    data=data %>% filter(!is.na(var)  & !is.na(y) & !is.na(covariate))
    ns=length(unique(data$studyid))
    nc=length(unique(data$condition))
    nts <- as.numeric(table(data$studyid))#***combination of number of timepoints and conditions by study
    max.nts <- max(nts)
    maxcov=dim(combn(max.nts-1,2))[2]
    data$condid<-as.numeric(as.factor(data$condition))
    data=data %>% group_by(record_id, condition) %>% group_modify(~ mutate(.x,tid=1:nrow(.x)))
    ntbycon=data %>% count(record_id, condid, condition)
    ncond=(ungroup(data %>% count(record_id, condition)) %>% count(record_id))$n #**number of reported conditions in each study
    data=data %>% group_by(record_id) %>% group_modify(~ add_count(.x)) #%>% group_modify(~ add_row(.x, d=rep(NA, maxcov - nrow(.x))))
    
    #Create the outcome data Y matrix
    Y=data %>% select(d,tid,condid,condition) %>% filter(tid>1)
    dimY=max(table(Y$record_id))
    Y=ungroup(Y%>% group_by(record_id)%>% group_modify(~ add_row(.x, d=rep(NA, dimY - nrow(.x)))) )
    #condmat=matrix(Y$condid,ncol=dimY,byrow=T)#**
    Y=matrix(Y$d,ncol=dimY,byrow=T)#**
    
    #Create the covariate matrix
    xvar=data %>% select(covariate,tid,condid,condition) %>% filter(tid>1)
    xvar=ungroup(xvar%>% group_by(record_id)%>% group_modify(~ add_row(.x, covariate=rep(NA, dimY - nrow(.x)))) )
    #condmat=matrix(Y$condid,ncol=dimY,byrow=T)#**
    covariate=matrix(xvar$covariate,ncol=dimY,byrow=T)#**
    
    
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
    }#makes the variance-covariances in a study with multiple timepoints and conditions. Assumes correlation 0.6 between condition SMDs
    
    Vlist=lapply(datasplit,create.covmat2)
    options(warn = defaultW)#bring warnings back
    
    prec<-lapply(Vlist,solve)
    prec<-lapply(prec,as.matrix)
    precmat<-plyr::ldply(prec,rbind)
    
    #Create a pseudo-dose matrix to change the dimensions as needed
    dimdose=max(table(data$record_id))
    dose=matrix(1,ncol=dimdose,nrow=ns)#*
  
    
    
    prepostJAGSdata <- list(Y = Y, nts = nts, dose=dose, ns = ns, prec = precmat[,-1], covariate=covariate)

  return(prepostJAGSdata)
}

#### 2  JAGS  META-REGRESSION  MODEL random effects for conditions---------------
prepostmetaregMHmodelRE1C<-function(){
        b1[1] <-1
        for (i in 1:ns) { ## for each study
          # multivariate normal likelihood of the response 
          Y[i,1:(nts[i]-1)]  ~ dmnorm(meanMH[i,1:(nts[i]-1)], prec[(b1[i]):(b1[i]+nts[i]-2),1:(nts[i]-1)])
          # window to change the index in precision matrix
          b1[i+1] <- b1[i]+ nts[i]-1
          meanMH[i,1:(nts[i]-1)]<-deltaMH[i]*dose[i,1:(nts[i]-1)]+gamma*covariate[i,1:(nts[i]-1)]
        }
        
        for(i in 1:ns){
          deltaMH[i]~dnorm(muMH,precMH)
        }
        
        #Priors for the pooled estimate
        muMH ~ dnorm(0,0.00001)
        gamma ~ dnorm(0,0.00001)
        precMH<-1/varianceMH
        varianceMH<-tauMH*tauMH
        tauMH~dnorm(0,1)%_%T(0,2)
        
        #prediction intervals
        predSMDMH~dnorm(muMH,precMH)
      }
