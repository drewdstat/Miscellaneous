#' Functions for simulating power for linear regressions or Weighted Quantile 
#' Sum (WQS) regressions 
#' 
#' Power calculations are commonly done to roughly estimate statistical power to 
#' detect a true association in planned statistical analyses given expected 
#' sample sizes, effect sizes, and other parameters. This package provides 3 
#' separate functions for simulating power for linear or logistic regressions, 
#' as well as for WQS regressions (Carrico et al. 2015), which I frequently used 
#' as mixture exposure regressions. The first function, \code{pwrfrompub}, takes 
#' a number of parameters typically reported in a publication and estimates 
#' power based on simulations involving those parameters. It's common that the 
#' only a priori information one might have regarding an analysis of interest is 
#' the often sparse information reported in a typical paper, and so this function
#' simulates a set of simple regressions based on what is typically reported.
#' 
#' @details
#' \code{pwrfrompub} utilizes estimates from a publication, namely regression 
#' coefficients, their confidence intervals, sample sizes, and approximate means 
#' and variances of the exposure of interest. These estimates are used to 
#' calculate an approximation of the residual standard deviation, using the 
#' following formulae (see 
#' \href{https://stats.stackexchange.com/questions/44838/how-are-the-standard-errors-of-coefficients-calculated-in-a-regression/44841#44841}{this Stackexchange thread} 
#' and \href{https://stats.stackexchange.com/questions/85943/how-to-derive-the-standard-error-of-linear-regression-coefficient}{this thread}):
#' \deqn{
#' SE_\hat{\beta} = {\frac {LCI_{\beta} - \beta}{- qt(0.975, df=N-2)}} \\
#' RSS = SE_\hat{\beta}^2 * SE_{exposure}^2 * (N-1) * (N-2) \\
#' SD_{residuals} = \sqrt{\frac{RSS}{N-1}}
#' }
#' These formulae are based on a univariable regression, but they should also 
#' be approximately correct with multivariable regression when the predictor 
#' variables are uncorrelated. With increasing positive correlation between 
#' predictors, the above formulae will overestimate the residual SD, meaning 
#' that power estimates made using those formulae will be more conservative than 
#' what may actually happen in the future analysis.
#' 

pwrfrompub<-function(niter,nobs_poss,betas_poss,nobs_study,beta_study,lci_study,Xmean,Xsd,
                     vc=NULL,method=c("Mine","Other"),powergoal=0.8,nparams_study=NULL,
                     family="gaussian",OR=F,logout=F,beta0=0,residsigma=NULL,wqs=F,wts,q=5,interp=T,
                     Xsd_study=NULL,binexp=F,p0=NULL){
  library(dplyr)
  library(paramtest)
  if(!is.null(nparams_study)){
    df<-nobs_study-nparams_study
  } else {
    df<-nobs_study-2
  }
  if(is.null(Xsd_study)){Xsd_study<-Xsd}
  betase_study<-(lci_study-beta_study)/-qt((0.95/2+0.5),df)
  RSS<-betase_study^2*Xsd_study^2*(nobs_study-1)*df #assumes no covariance between predictors
  if(is.null(residsigma)){
    residualsd<-sqrt(RSS/(nobs_study-1))
    #as positive correlation between predictors increases, this estimate of residSD overestimates,
    #  leading to more conservative (lower) power estimates
  } else {
    residualsd<-residsigma
  }
  
  Resultsmat<-expand.grid(N=nobs_poss,EffectSize=betas_poss) #50 combos
  Resultsmat[,c("Power")]<-NA
  lmpowersim<-function(niter,nobs,residsd,expmean,expsd,betas,vc=NULL){
    pvals<-as.data.frame(matrix(NA,niter,length(betas)))
    seeds<-c()
    if(is.null(vc)){
      vc<-diag(x=expsd^2,length(expsd))
    }
    for(i in 1:niter){
      newseed<-sample(1:1E6,1)
      seeds[i]<-newseed
      set.seed(newseed)
      X<-mvtnorm::rmvnorm(nobs,expmean,vc)
      set.seed(newseed-1)
      eps<-rnorm(nobs,0,residsd)
      y<-cbind(1,X)%*%betas+eps
      moddat<-as.data.frame(cbind(y,X))
      names(moddat)<-c("y",paste0("x",1:length(expmean)))
      mod<-lm(y~.,moddat)
      pvals[i,]<-summary(mod)$coef[,4]
    }
    names(pvals)<-paste0("beta",0:(length(betas)-1))
    return(list(Pvals=pvals,Seeds=seeds))
  }
  lmpowersim_single<-function(simNum,nobs,residsd,expmean,expsd,b1,b0=0,fam="gaussian",binX=F){
    if(binX){
      X<-sample(c(0,1),nobs,replace=T,prob=c(1-expmean,expmean))
      if(all(X==0)){X<-sample(c(0,1),nobs,replace=T,prob=c(1-expmean,expmean))}
    } else {
      X<-rnorm(nobs,expmean,expsd)
    }
    
    eps<-rnorm(nobs,0,residsd)
    if(fam=="gaussian"){
      y<-cbind(1,X)%*%c(b0,b1)+eps
    } else if(fam=="binomial"){
      eta<-cbind(1,X)%*%c(b0,b1)#+eps
      prob_y<-1/(1+exp(-eta))
      y<-rbinom(nobs,size=1,prob=prob_y)
    } else if(class(fam)=="family"&fam$link=="log"){
      zeta0<-((1-p0)*exp(b0))/(1-exp(b0)*p0)
      zeta1<-((1-p0)*exp(b1))/(1-exp(b1)*p0)
      eta<-cbind(1,X)%*%c(log(zeta0),log(zeta1))
      prob_y<-1/(1+exp(-eta))
      y<-rbinom(nobs,size=1,prob=prob_y)
    } else { 
      stop("family must be gaussian, binomial, or binomial(link=log)")
    }
    
    moddat<-as.data.frame(cbind(y,X))
    names(moddat)<-c("y",paste0("x",1:length(expmean)))
    if(fam=="gaussian"){
      model<-lm(y~.,moddat)
    } else if(fam=="binomial"){
      model<-glm(y~.,data=moddat,family="binomial")
    } else if(class(fam)=="family"&fam$link=="log"){
      model<-glm(y~.,data=moddat,family=fam)
    }
    
    est <- coef(summary(model))['x1', 'Estimate']
    se <- coef(summary(model))['x1', 'Std. Error']
    p <- coef(summary(model))['x1', grep("Pr(>",colnames(coef(summary(model))),fixed=T)]
    return(c(xm=mean(X), xsd=sd(X), ym=mean(y), ysd=sd(y), est=est, se=se, p=p,
             sig=(p < .05)))
  }
  lm_test <- function(simNum, N, b1, b0=0, xm=0, xsd=1) {
    x <- rnorm(N, xm, xsd)
    y <- rnorm(N, b0 + b1*x, sqrt(1 - b1^2))  # var. approx. 1 after accounting
    # for explained variance by x
    model <- lm(y ~ x)
    
    # pull output from model
    est <- coef(summary(model))['x', 'Estimate']
    se <- coef(summary(model))['x', 'Std. Error']
    p <- coef(summary(model))['x', 'Pr(>|t|)']
    
    return(c(xm=mean(x), xsd=sd(x), ym=mean(y), ysd=sd(y), est=est, se=se, p=p,
             sig=(p < .05)))
  }
  wqspowersim<-function(simNum,nobs,residsd,expmeans,VC,w,qq,b1,b0=0,fam="gaussian"){
    #library(gWQS)
    X<-mvtnorm::rmvnorm(nobs,expmeans,VC)
    quantile_transform<-function(mixture,nq){
      for(i in 1:ncol(mixture)){
        newbreaks<-unique(quantile(mixture[,i],probs=seq(0,1,by=1/nq),na.rm=T))
        mixture[,i]<-cut(mixture[,i],breaks=newbreaks,labels=F,include.lowest=T)-1
      }
      return(mixture)
    }
    Xq<-quantile_transform(X,nq=qq)
    WQS<-Xq%*%w
    eps<-rnorm(nobs,0,residsd)
    if(fam=="gaussian"){
      y<-cbind(1,WQS)%*%c(b0,b1)+eps
    } else if(fam=="binomial"){
      eta<-cbind(1,WQS)%*%c(b0,b1)#+eps
      prob_y<-1/(1+exp(-eta))
      y<-rbinom(nobs,size=1,prob=prob_y)
      #prob_y<-exp(eta)/(1+exp(eta))
      #runis<-runif(nobs,0,1)
      #y<-ifelse(runis<prob_y,1,0)
    } else {
      stop("family must be gaussian or binomial")
    }
    
    moddat<-as.data.frame(cbind(y,WQS))
    names(moddat)<-c("y","wqs")
    if(fam=="gaussian"){
      model<-lm(y~.,moddat)
    } else if(fam=="binomial"){
      model<-glm(y~.,moddat,family="binomial")
    }
    
    est <- coef(summary(model))['wqs', 'Estimate']
    se <- coef(summary(model))['wqs', 'Std. Error']
    p <- coef(summary(model))['wqs', grep("Pr(>",colnames(coef(summary(model))),fixed=T)]
    return(c(xm=mean(X), xsd=sd(X), ym=mean(y), ysd=sd(y), est=est, se=se, p=p,
             sig=(p < .05)))
  }
  if(OR){betas_poss<-log(betas_poss)}
  if(logout){betas_poss<-log((betas_poss/100)+1)}
  if(wqs){
    power_lm<-grid_search(wqspowersim,params=list(nobs=nobs_poss,b1=betas_poss),n.iter=niter,
                          output="data.frame",parallel="snow",ncpus=parallel::detectCores(),
                          residsd=residualsd,expmeans=Xmean,VC=vc,w=wts,fam=family,b0=beta0,qq=q)
    pwrsumm<-results(power_lm) %>%
      group_by(nobs.test,b1.test) %>%
      summarise(power=mean(sig))
  } else if(method=="Mine"){
    power_lm<-grid_search(lmpowersim_single,params=list(nobs=nobs_poss,b1=betas_poss),n.iter=niter,
                          output="data.frame",parallel="snow",ncpus=parallel::detectCores(),
                          residsd=residualsd,expmean=Xmean,expsd=Xsd,fam=family,b0=beta0,binX=binexp)
    pwrsumm<-results(power_lm) %>%
      group_by(nobs.test,b1.test) %>%
      summarise(power=mean(sig))
  } else if(method=="Other"){
    power_lm <- grid_search(lm_test, params=list(N=nobs_poss,b1=betas_poss),n.iter=niter,
                            output='data.frame',xm=Xmean,xsd=Xsd,b0=beta0,parallel='snow',
                            ncpus=parallel::detectCores())
    pwrsumm<-results(power_lm) %>%
      group_by(N.test,b1.test) %>%
      summarise(power=mean(sig))
  } else {
    stop("Option 'method' must be either 'Mine' or 'Other'.")
  }
  names(pwrsumm)<-c("N","Beta","Power")
  pwrsumm$PwrLabel<-as.factor(ifelse(round(pwrsumm$Power,2)>=powergoal,"Y","N"))
  if(all(pwrsumm$PwrLabel=="N")) colvals<-"white" else if(all(pwrsumm$PwrLabel=="Y")) colvals<-"black" else
    colvals<-c("white","black")
  if(OR){
    pwrsumm$OR<-as.factor(exp(pwrsumm$Beta))
    g1<-ggplot(pwrsumm,aes(x=N,y=OR))+theme_classic()+
      geom_tile(aes(fill=Power),alpha=0.5)+
      geom_text(aes(label=round(Power,2),color=PwrLabel))+
      #scale_y_continuous(breaks=unique(pwrsumm$OR))+
      scale_color_manual(values=colvals,guide="none")+
      ylab("Odds Ratio")
    pwrsumm$OR<-NULL
  } else if(logout){
    pwrsumm$PercBeta<-as.factor(round((exp(pwrsumm$Beta)-1)*100,2))
    g1<-ggplot(pwrsumm,aes(x=N,y=PercBeta))+theme_classic()+
      geom_tile(aes(fill=Power),alpha=0.5)+
      geom_text(aes(label=round(Power,2),color=PwrLabel))+
      #scale_y_continuous(breaks=unique(pwrsumm$PercBeta))+
      scale_color_manual(values=colvals,guide="none")+
      ylab("Percent Change")
    pwrsumm$PercBeta<-NULL
  } else {
    pwrsumm$BetaFac<-as.factor(pwrsumm$Beta)
    g1<-ggplot(pwrsumm,aes(x=N,y=BetaFac))+theme_classic()+
      geom_tile(aes(fill=Power),alpha=0.5)+
      geom_text(aes(label=round(Power,2),color=PwrLabel))+
      #scale_y_continuous(breaks=unique(pwrsumm$Beta))+
      scale_color_manual(values=colvals,guide="none")+
      ylab("Effect Size")
    pwrsumm$BetaFac<-NULL
  }
  pwrsumm$PwrLabel<-NULL
  
  if(interp){
    library(akima)
    s=interp(pwrsumm$N,pwrsumm$Beta,pwrsumm$Power,
             xo=c(seq(min(pwrsumm$N),median(pwrsumm$N),length=25),
                  seq(median(pwrsumm$N)+(max(pwrsumm$N)-min(pwrsumm$N))/25,max(pwrsumm$N),length=25)),
             yo=seq(min(pwrsumm$Beta),max(pwrsumm$Beta),length=100))
    s1 <- interp2xyz(s)
  } else {
    s1<-NULL
  }
  
  #dim(s1)
  #s1[round(s1[,"z"],3)==0.80,c("x","y","z")]
  return(list(Summary=pwrsumm,GGplot=g1,Interp=s1))
}

pwrfrommod<-function(niter,nobs_poss,betas_poss,model,powergoal=0.8,OR=F,logout=F){
  residualsd<-sd(model$residuals)
  modX<-model.matrix(model$formula,model$model)[,-1]
  modXmeans<-colMeans(modX,na.rm=T)
  modVC<-cov(modX)
  uniques<-sapply(as.data.frame(modX),function(x) unique(x))
  is.binary<-sapply(uniques,function(x) all(x%in%c(0,1)))
  if("family"%in%names(model)){
    family<-as.character(model$family)[1]
  } else {family<-"gaussian"}
  if(OR){betas_poss<-log(betas_poss)}
  if(logout){betas_poss<-log((betas_poss/100)+1)}
  lmpowersim<-function(simNum,nobs,residsd,expmeans,vc,betas,b1,binarybool,fam="gaussian"){
    X<-mvtnorm::rmvnorm(nobs,expmeans,vc)
    if(any(binarybool)){
      for(i in which(binarybool)){
        X[,i]<-ifelse(X[,i]>quantile(X[,i],probs=1-expmeans[i]),1,0)
      }; rm(i)
    }
    betas[2]<-b1
    eps<-rnorm(nobs,0,residsd)
    if(fam=="gaussian"){
      y<-cbind(1,X)%*%betas+eps
    } else if(fam=="binomial"){
      eta<-cbind(1,X)%*%betas#+eps
      prob_y<-1/(1+exp(-eta))
      y<-rbinom(nobs,size=1,prob=prob_y)
    } else {
      stop("family must be gaussian or binomial")
    }
    
    moddat<-as.data.frame(cbind(y,X))
    names(moddat)<-c("y",paste0("x",1:length(expmeans)))
    if(fam=="gaussian"){
      model<-lm(y~.,moddat)
    } else if(fam=="binomial"){
      model<-glm(y~.,moddat,family="binomial")
    }
    
    est <- coef(summary(model))['x1', 'Estimate']
    se <- coef(summary(model))['x1', 'Std. Error']
    p <- coef(summary(model))['x1', grep("Pr(>",colnames(coef(summary(model))),fixed=T)]
    return(c(xm=mean(X), xsd=sd(X), ym=mean(y), ysd=sd(y), est=est, se=se, p=p,
             sig=(p < .05)))
  }
  power_lm<-grid_search(lmpowersim,params=list(nobs=nobs_poss,b1=betas_poss),n.iter=niter,
                        output="data.frame",parallel="snow",ncpus=parallel::detectCores(),
                        residsd=residualsd,expmeans=modXmeans,vc=modVC,betas=model$coefficients,
                        fam=family,binarybool=is.binary)
  pwrsumm<-results(power_lm) %>%
    group_by(nobs.test,b1.test) %>%
    summarise(power=mean(sig))
  
  names(pwrsumm)<-c("N","Beta","Power")
  pwrsumm$PwrLabel<-as.factor(ifelse(round(pwrsumm$Power,2)>=powergoal,"Y","N"))
  if(all(pwrsumm$PwrLabel=="N")) colvals<-"white" else if(all(pwrsumm$PwrLabel=="Y")) colvals<-"black" else
    colvals<-c("white","black")
  if(OR){
    pwrsumm$OR<-exp(pwrsumm$Beta)
    pwrsumm$OR<-as.factor(pwrsumm$OR)
    g1<-ggplot(pwrsumm,aes(x=N,y=OR))+theme_classic()+
      geom_tile(aes(fill=Power),alpha=0.5)+
      geom_text(aes(label=round(Power,2),color=PwrLabel))+
      #scale_y_continuous(breaks=unique(pwrsumm$OR))+
      scale_color_manual(values=colvals,guide="none")+
      ylab("Odds Ratio")
    pwrsumm$OR<-NULL
  } else if(logout){
    pwrsumm$PercBeta<-(exp(pwrsumm$Beta)-1)*100
    pwrsumm$PercBeta<-as.factor(pwrsumm$PercBeta)
    g1<-ggplot(pwrsumm,aes(x=N,y=PercBeta))+theme_classic()+
      geom_tile(aes(fill=Power),alpha=0.5)+
      geom_text(aes(label=round(Power,2),color=PwrLabel))+
      #scale_y_continuous(breaks=unique(pwrsumm$PercBeta))+
      scale_color_manual(values=colvals,guide="none")+
      ylab("Percent Change")
    pwrsumm$PercBeta<-NULL
  } else {
    pwrsumm$BetaFac<-as.factor(pwrsumm$Beta)
    g1<-ggplot(pwrsumm,aes(x=N,y=BetaFac))+theme_classic()+
      geom_tile(aes(fill=Power),alpha=0.5)+
      geom_text(aes(label=round(Power,2),color=PwrLabel))+
      #scale_y_continuous(breaks=unique(pwrsumm$Beta))+
      scale_color_manual(values=colvals,guide="none")+
      ylab("Effect Size")
    pwrsumm$BetaFac<-NULL
  }
  pwrsumm$PwrLabel<-NULL
  
  library(akima)
  s=interp(pwrsumm$N,pwrsumm$Beta,pwrsumm$Power,
           xo=c(seq(min(pwrsumm$N),median(pwrsumm$N),length=25),
                seq(median(pwrsumm$N)+(max(pwrsumm$N)-min(pwrsumm$N))/25,max(pwrsumm$N),length=25)),
           yo=seq(min(pwrsumm$Beta),max(pwrsumm$Beta),length=100))
  s1 <- interp2xyz(s)
  return(list(Summary=pwrsumm,GGplot=g1,Interp=s1))
}

pwrfromtheory<-function(niter=1000,nobs_poss=seq(500,2000,500),
                        betas_poss=seq(0.1,0.5,0.2),residsd_poss=1:3,Xmean=0,Xsd=1,
                        powergoal=0.8,family="gaussian",OR=F,
                        logout=F,beta0=0,wqs=F,wts,q=5,interp=T,
                        binexp=F,p0=NULL,ix=F,b2_poss=seq(0.1,0.5,0.2),
                        b3_poss=seq(0.1,0.5,0.2)){
  library(dplyr)
  library(paramtest)
  Resultsmat<-expand.grid(N=nobs_poss,EffectSize=betas_poss,ResidualSD=residsd_poss)
  Resultsmat[,c("Power")]<-NA
  
  lmpowersim_single<-function(simNum,nobs,residsd,expmean,expsd,b1,b0=0,fam="gaussian",binX=F){
    if(binX){
      X<-sample(c(0,1),nobs,replace=T,prob=c(1-expmean,expmean))
      if(all(X==0)){X<-sample(c(0,1),nobs,replace=T,prob=c(1-expmean,expmean))}
    } else {
      X<-rnorm(nobs,expmean,expsd)
    }
    
    eps<-rnorm(nobs,0,residsd)
    if(fam=="gaussian"){
      y<-cbind(1,X)%*%c(b0,b1)+eps
    } else if(fam=="binomial"){
      eta<-cbind(1,X)%*%c(b0,b1)#+eps
      prob_y<-1/(1+exp(-eta))
      y<-rbinom(nobs,size=1,prob=prob_y)
    } else if(class(fam)=="family"&fam$link=="log"){
      zeta0<-((1-p0)*exp(b0))/(1-exp(b0)*p0)
      zeta1<-((1-p0)*exp(b1))/(1-exp(b1)*p0)
      eta<-cbind(1,X)%*%c(log(zeta0),log(zeta1))
      prob_y<-1/(1+exp(-eta))
      y<-rbinom(nobs,size=1,prob=prob_y)
    } else { 
      stop("family must be gaussian, binomial, or binomial(link=log)")
    }
    
    moddat<-as.data.frame(cbind(y,X))
    names(moddat)<-c("y",paste0("x",1:length(expmean)))
    if(fam=="gaussian"){
      model<-lm(y~.,moddat)
    } else if(fam=="binomial"){
      model<-glm(y~.,data=moddat,family="binomial")
    } else if(class(fam)=="family"&fam$link=="log"){
      model<-glm(y~.,data=moddat,family=fam)
    }
    
    est <- coef(summary(model))['x1', 'Estimate']
    se <- coef(summary(model))['x1', 'Std. Error']
    p <- coef(summary(model))['x1', grep("Pr(>",colnames(coef(summary(model))),fixed=T)]
    return(c(xm=mean(X), xsd=sd(X), ym=mean(y), ysd=sd(y), est=est, se=se, p=p,
             sig=(p < .05)))
  }
  lmpowersim_contix<-function(simNum,nobs,residsd,X1mean,X1sd,X2mean,X2sd,b1,b2,b3,b0=0,fam="gaussian"){
    x1<-rnorm(nobs,X1mean,X1sd)
    x2<-rnorm(nobs,X2mean,X2sd)
    x1x2<-x1*x2
    X<-cbind(1,x1,x2,x1x2)
    eps<-rnorm(nobs,0,residsd)
    if(fam=="gaussian"){
      y<-X%*%c(b0,b1,b2,b3)+eps
    } else if(fam=="binomial"){
      eta<-X%*%c(b0,b1,b2,b3)#+eps
      prob_y<-1/(1+exp(-eta))
      y<-rbinom(nobs,size=1,prob=prob_y)
    } else if(class(fam)=="family"&fam$link=="log"){
      zeta0<-((1-p0)*exp(b0))/(1-exp(b0)*p0)
      zeta1<-((1-p0)*exp(b1))/(1-exp(b1)*p0)
      zeta2<-((1-p0)*exp(b2))/(1-exp(b2)*p0)
      zeta3<-((1-p0)*exp(b3))/(1-exp(b3)*p0)
      eta<-X%*%c(log(zeta0),log(zeta1),log(zeta2),log(zeta3))
      prob_y<-1/(1+exp(-eta))
      y<-rbinom(nobs,size=1,prob=prob_y)
    } else { 
      stop("family must be gaussian, binomial, or binomial(link=log)")
    }
    
    moddat<-as.data.frame(cbind(y,X[,c(2:3)]))
    names(moddat)<-c("y",paste0("x",1:2))
    if(fam=="gaussian"){
      model<-lm(y~x1*x2,moddat)
    } else if(fam=="binomial"){
      model<-glm(y~x1*x2,data=moddat,family="binomial")
    } else if(class(fam)=="family"&fam$link=="log"){
      model<-glm(y~x1*x2,data=moddat,family=fam)
    }
    
    est <- coef(summary(model))['x1:x2', 'Estimate']
    se <- coef(summary(model))['x1:x2', 'Std. Error']
    p <- coef(summary(model))['x1:x2', grep("Pr(>",colnames(coef(summary(model))),fixed=T)]
    return(c(x1m=mean(x1), x1sd=sd(x1), x2m=mean(x2), x2sd=sd(x2), ym=mean(y), ysd=sd(y), 
             est=est, se=se, p=p, sig=(p < .05)))
  }
  wqspowersim<-function(simNum,nobs,residsd,expmeans,VC,w,qq,b1,b0=0,fam="gaussian"){
    #library(gWQS)
    X<-mvtnorm::rmvnorm(nobs,expmeans,VC)
    quantile_transform<-function(mixture,nq){
      for(i in 1:ncol(mixture)){
        newbreaks<-unique(quantile(mixture[,i],probs=seq(0,1,by=1/nq),na.rm=T))
        mixture[,i]<-cut(mixture[,i],breaks=newbreaks,labels=F,include.lowest=T)-1
      }
      return(mixture)
    }
    Xq<-quantile_transform(X,nq=qq)
    WQS<-Xq%*%w
    eps<-rnorm(nobs,0,residsd)
    if(fam=="gaussian"){
      y<-cbind(1,WQS)%*%c(b0,b1)+eps
    } else if(fam=="binomial"){
      eta<-cbind(1,WQS)%*%c(b0,b1)#+eps
      prob_y<-1/(1+exp(-eta))
      y<-rbinom(nobs,size=1,prob=prob_y)
      #prob_y<-exp(eta)/(1+exp(eta))
      #runis<-runif(nobs,0,1)
      #y<-ifelse(runis<prob_y,1,0)
    } else {
      stop("family must be gaussian or binomial")
    }
    
    moddat<-as.data.frame(cbind(y,WQS))
    names(moddat)<-c("y","wqs")
    if(fam=="gaussian"){
      model<-lm(y~.,moddat)
    } else if(fam=="binomial"){
      model<-glm(y~.,moddat,family="binomial")
    }
    
    est <- coef(summary(model))['wqs', 'Estimate']
    se <- coef(summary(model))['wqs', 'Std. Error']
    p <- coef(summary(model))['wqs', grep("Pr(>",colnames(coef(summary(model))),fixed=T)]
    return(c(xm=mean(X), xsd=sd(X), ym=mean(y), ysd=sd(y), est=est, se=se, p=p,
             sig=(p < .05)))
  }
  if(OR){betas_poss<-log(betas_poss)}
  if(logout){betas_poss<-log((betas_poss/100)+1)}
  if(wqs){
    power_lm<-grid_search(wqspowersim,params=list(nobs=nobs_poss,b1=betas_poss,residsd=residsd_poss),
                          n.iter=niter,output="data.frame",parallel="snow",ncpus=parallel::detectCores(),
                          expmeans=Xmean,VC=vc,w=wts,fam=family,b0=beta0,qq=q)
    pwrsumm<-results(power_lm) %>%
      group_by(nobs.test,b1.test,residsd.test) %>%
      summarise(power=mean(sig))
    names(pwrsumm)<-c("N","Beta","ResidSD","Power")
  } else if(ix) {
    if(OR){b2_poss<-log(b2_poss)
      b3_poss<-log(b3_poss)
    }
    power_lm<-grid_search(lmpowersim_contix,
                          params=list(nobs=nobs_poss,b1=betas_poss,b2=b2_poss,b3=b3_poss,
                                      residsd=residsd_poss),n.iter=niter,
                          output="data.frame",parallel="snow",ncpus=parallel::detectCores(),
                          X1mean=Xmean,X1sd=Xsd,X2mean=Xmean,X2sd=Xsd,fam=family,b0=beta0)
    pwrsumm<-results(power_lm) %>%
      group_by(nobs.test,b1.test,b2.test,b3.test,residsd.test) %>%
      summarise(power=mean(sig))
    names(pwrsumm)<-c("N","Beta1","Beta2","Beta3","ResidSD","Power")
  } else {
    power_lm<-grid_search(lmpowersim_single,
                          params=list(nobs=nobs_poss,b1=betas_poss,residsd=residsd_poss),n.iter=niter,
                          output="data.frame",parallel="snow",ncpus=parallel::detectCores(),
                          expmean=Xmean,expsd=Xsd,fam=family,b0=beta0,binX=binexp)
    pwrsumm<-results(power_lm) %>%
      group_by(nobs.test,b1.test,residsd.test) %>%
      summarise(power=mean(sig))
    names(pwrsumm)<-c("N","Beta","ResidSD","Power")
  }
  pwrsumm$PwrLabel<-as.factor(ifelse(round(pwrsumm$Power,2)>=powergoal,"Y","N"))
  if(all(pwrsumm$PwrLabel=="N")) colvals<-"white" else if(all(pwrsumm$PwrLabel=="Y")) colvals<-"black" else
    colvals<-c("white","black")
  pwrsumm$ResidSDlabel<-as.factor(paste0("ResidSD=",pwrsumm$ResidSD))
  if(length(residsd_poss)==4) numplotcols<-2 else numplotcols<-min(c(3,length(residsd_poss)))
  if(!ix){
    if(OR){
      pwrsumm$OR<-as.factor(exp(pwrsumm$Beta))
      g1<-ggplot(pwrsumm,aes(x=N,y=OR))+theme_classic()+
        geom_tile(aes(fill=Power),alpha=0.5)+
        geom_text(aes(label=round(Power,2),color=PwrLabel))+
        facet_wrap(.~ResidSDlabel,ncol=numplotcols)+
        scale_color_manual(values=colvals,guide="none")+
        ylab("Odds Ratio")
      pwrsumm$OR<-NULL
    } else if(logout){
      pwrsumm$PercBeta<-as.factor(round((exp(pwrsumm$Beta)-1)*100,2))
      g1<-ggplot(pwrsumm,aes(x=N,y=PercBeta))+theme_classic()+
        geom_tile(aes(fill=Power),alpha=0.5)+
        geom_text(aes(label=round(Power,2),color=PwrLabel))+
        facet_wrap(.~ResidSDlabel,ncol=numplotcols)+
        scale_color_manual(values=colvals,guide="none")+
        ylab("Percent Change")
      pwrsumm$PercBeta<-NULL
    } else {
      pwrsumm$BetaFac<-as.factor(pwrsumm$Beta)
      g1<-ggplot(pwrsumm,aes(x=N,y=BetaFac))+theme_classic()+
        geom_tile(aes(fill=Power),alpha=0.5)+
        geom_text(aes(label=round(Power,2),color=PwrLabel))+
        facet_wrap(.~ResidSDlabel,ncol=numplotcols)+
        scale_color_manual(values=colvals,guide="none")+
        ylab("Effect Size")
      pwrsumm$BetaFac<-NULL
    }
  } else {g1<-NULL}
  pwrsumm$PwrLabel<-NULL
  
  if(interp&!ix){
    library(akima)
    slist<-list()
    for(ee in 1:length(residsd_poss)){
      pwrtemp<-pwrsumm[which(pwrsumm$ResidSD==residsd_poss[ee]),]
      s<-interp(pwrtemp$N,pwrtemp$Beta,pwrtemp$Power,
                xo=c(seq(min(pwrtemp$N),median(pwrtemp$N),length=25),
                     seq(median(pwrtemp$N)+(max(pwrtemp$N)-min(pwrtemp$N))/25,max(pwrtemp$N),length=25)),
                yo=seq(min(pwrtemp$Beta),max(pwrtemp$Beta),length=100))
      slist[[paste0("ResidSD=",residsd_poss[ee])]] <- interp2xyz(s)
    }
  } else {
    slist<-NULL
  }
  return(list(Summary=pwrsumm,GGplot=g1,Interp=slist))
}
