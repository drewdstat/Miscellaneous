Resmat_lm<-function(prednames,outnames,covnames=NULL,Data,logout=F,logpred=F,
                    logbasepred=10,logbaseout=exp(1),Outtitle="Outcome",
                    Predtitle="Exposure",MyMult=1,ixterm=NULL,Firth=F,binomlink="logit",binfunc=binomial,
                    robust=T,HCtype="HC0",predspline=F,predsplinedf=3,plotOR=F,plotPercChange=T,LOOCV=F,
                    facetcol=NULL,covix=NULL,ixpred=T,extradiag=T,leverage.test=F,leverage.cutoff=0.2,
                    leverage.meancutoff=NULL,post.power=F,effect.size=0.5,nsim=1E3,
                    mice=F,micevars=NULL,miceiter=10){
  library(caret)
  library(dplyr)
  library(ggplot2)
  library(kableExtra)
  library(logistf)
  library(sandwich)
  library(splines)
  library(emmeans)
  robustse <- function(x, HCtype="HC0") {
    suppressMessages(suppressWarnings(library(lmtest)))
    suppressMessages(suppressWarnings(library(sandwich)))
    mod1 <- coeftest(x, vcov = function(x) vcovHC(x,type=HCtype))
    if(any(class(x)=="glm")){
      cis<-coefci(x, vcov=function(x) vcovHC(x,type=HCtype))
    } else {
      cis<-coefci(x, vcov = function(x) vcovHC(x,type=HCtype),
                  df=x$df.residual)
    }
    mod1<-cbind(mod1,cis)
    return(mod1)
  }
  confint.qt<-function(beta,se,DF,IQR=1,level=0.95){
    ci.lower<-(beta*IQR)-((se*IQR)*qt(((level/2)+0.5), DF))
    ci.upper<-(beta*IQR)+((se*IQR)*qt(((level/2)+0.5), DF))
    CIs<-data.frame(CIL=ci.lower,CIU=ci.upper)
    return(CIs)
  }
  if(!is.null(ixterm)&predspline){
    stop("Code not yet set up to accommodate spline interactions.")
  }
  
  if(robust&Firth){
    stop("Code not yet set up to extract robust sandwich SEs from logistf objects.")
  }
  if(Firth&extradiag){
    stop("Code not yet set up to perform GLM diagnostics on logistf objects.")
  }
  
  predclasses<-sapply(Data[,prednames],class)
  if(any(predclasses=="character")){
    for(cl in which(predclasses=="character")){
      Data[,prednames[cl]]<-as.factor(Data[,prednames[cl]])
    }; rm(cl)
  }
  
  if(any(predclasses%in%c("factor"))){
    predlevellengths<-sapply(Data[,prednames],function(x) length(levels(x)))
    predlevelrows<-ifelse(predlevellengths%in%c(0,1),1,predlevellengths-1)
    predrowlength<-sum(predlevelrows)
  } else {
    predlevelrows<-rep(1,length(prednames))
    predrowlength<-length(prednames)
  }
  
  if(is.null(ixterm)&predspline!=F){
    Resultsmat<-as.data.frame(matrix(NA,predsplinedf*predrowlength*length(outnames),12))
  } else if(!is.null(ixterm)){
    nixlvl<-length(levels(Data[,ixterm]))
    if(class(Data[,ixterm])%in%c("factor","character")){
      if(class(Data[,ixterm])=="character"){Data[,ixterm]<-as.factor(Data[,ixterm])}
      if(is.null(ncol(combn(levels(Data[,ixterm]),2)))) {
        ixtermlength<-nixlvl+(nixlvl-1)+1
      } else {
        ixtermlength<-nixlvl+(nixlvl-1)+
          ncol(combn(levels(Data[,ixterm]),2))
      }
    } else {
      ixtermlength<-3
    }
    
    Resultsmat<-as.data.frame(matrix(NA,ixtermlength*predrowlength*length(outnames),
                                     12))
  } else {
    Resultsmat<-as.data.frame(matrix(NA,predrowlength*length(outnames),11))
  }
  
  lmlist<-list()
  
  if(extradiag){Extramat<-data.frame(Heterosk.p=rep(NA,nrow(Resultsmat)))}
  
  if(length(logout)==1){
    logout2<-rep(logout,length(outnames))
  } else {
    logout2<-logout
  }
  
  if(length(logpred)==1){
    logpred2<-rep(logpred,length(prednames))
  } else {
    logpred2<-logpred
  }
  
  for(i in 1:length(outnames)){
    for(j in 1:length(prednames)){
      if(is.null(ixterm)&predspline==T){
        rowstart<-(1+(j-1))+((i-1)*predsplinedf*predrowlength)+((predsplinedf-1)*(j-1)) 
        rowend<-(predsplinedf+(j-1))+((i-1)*predsplinedf*predrowlength)+ #length(prednames)
          ((predsplinedf-1)*(j-1))
        rownum<-rowstart:rowend
        lmnum<-(1+(j-1))+((i-1)*length(prednames))
      } else if(!is.null(ixterm)){
        rowstart<-(1+(j-1))+((i-1)*ixtermlength*predrowlength)+((ixtermlength-1)*(j-1))
        rowend<-(ixtermlength+(j-1))+((i-1)*ixtermlength*predrowlength)+ #length(prednames)
          ((ixtermlength-1)*(j-1))
        rownum<-rowstart:rowend
        lmnum<-(1+(j-1))+((i-1)*predrowlength) #length(prednames)
      } else if(any(sapply(Data[,prednames],FUN=function(x) class(x)=="factor"))) {
        rowstarts<-cumsum(c(1,predlevelrows[-length(predlevelrows)]))+((i-1)*predrowlength)
        rowends<-cumsum(predlevelrows)+((i-1)*predrowlength)
        rownum<-rowstarts[j]:rowends[j]
        lmnum<-(1+(j-1))+((i-1)*length(prednames))
      } else {
        rownum<-(1+(j-1))+((i-1)*length(prednames))
        lmnum<-rownum
      }
      
      if(logout2[i]==T){
        outexpr<-paste0("log(",outnames[i],",base=",logbaseout,")")
      } else {
        outexpr<-outnames[i]
      }
      
      if(logpred2[j]==T&predspline==T){
        predexpr<-paste0("ns(log(",prednames[j],",base=",logbasepred,"),",predsplinedf,")")
      } else if(logpred2[j]==T){
        predexpr<-paste0("log(",prednames[j],",base=",logbasepred,")")
      } else if(predspline==T){
        predexpr<-paste0("ns(",prednames[j],",",predsplinedf,")")
      } else {
        predexpr<-prednames[j]
      }
      
      if(is.null(covnames)){
        if(is.null(ixterm)){
          if(mice){
            ccvars<-c(outnames[i],prednames[j])
            if(is.null(micevars)){micevars<-ccvars[which(sapply(Data[,ccvars],function(x) any(is.na(x))))]}
            ccDat<-Data[complete.cases(Data[,ccvars[-which(ccvars%in%micevars)]]),]
            micedata<-mice(ccDat[,ccvars],m=miceiter,maxit=miceiter,printFlag=F)
          } else {
            ccDat<-Data[complete.cases(Data[,c(outnames[i],prednames[j])]),]
          }
          nobs<-dim(ccDat)[1]
          form1<-formula(paste0(outexpr,"~",predexpr))
        } else {
          if(mice){
            ccvars<-c(outnames[i],prednames[j],ixterm)
            if(is.null(micevars)){micevars<-ccvars[which(sapply(Data[,ccvars],function(x) any(is.na(x))))]}
            ccDat<-Data[complete.cases(Data[,ccvars[-which(ccvars%in%micevars)]]),]
            micedata<-mice(ccDat[,ccvars],m=miceiter,maxit=miceiter,printFlag=F)
          } else {
            ccDat<-Data[complete.cases(Data[,c(outnames[i],prednames[j],ixterm)]),]
          }
          nobs<-dim(ccDat)[1]
          form1<-formula(paste0(outexpr,"~",predexpr,"*",ixterm))
        }
      } else {
        if(is.list(covnames)){
          if(length(covnames)!=(length(prednames)*length(outnames))){
            stop(paste0("Error: covnames must have a length of ", 
                        length(prednames)*length(outnames)))
          }
          mycovnames<-covnames[[lmnum]]
        } else {
          mycovnames<-covnames
        }
        if(any(grepl("(",mycovnames,fixed=T))&any(grepl(":",mycovnames,fixed=T))){
          mycovnames_df<-gsub("ns(","",mycovnames,fixed=T)
          mycovnames_df<-gsub(",.*","",mycovnames_df)
          mycovnames_df<-mycovnames_df[-which(grepl(":",mycovnames_df,fixed=T))]
        } else if(any(grepl(":",mycovnames,fixed=T))){
          mycovnames_df<-mycovnames[-which(grepl(":",mycovnames,fixed=T))]
        } else if(any(grepl("(",mycovnames,fixed=T))){
          mycovnames_df<-gsub("ns(","",mycovnames,fixed=T)
          mycovnames_df<-gsub(",.*","",mycovnames_df)
        } else {
          mycovnames_df<-mycovnames
        }
        if(is.null(ixterm)&is.null(covix)){
          if(mice){
            ccvars<-c(outnames[i],prednames[j],mycovnames_df)
            if(is.null(micevars)){micevars<-ccvars[which(sapply(Data[,ccvars],function(x) any(is.na(x))))]}
            ccDat<-Data[complete.cases(Data[,ccvars[-which(ccvars%in%micevars)]]),]
            micedata<-mice(ccDat[,ccvars],m=miceiter,maxit=miceiter,printFlag=F)
          } else {
            ccDat<-Data[complete.cases(Data[,c(outnames[i],prednames[j],mycovnames_df)]),]
          }
          nobs<-dim(ccDat)[1]
          form1<-formula(paste0(outexpr,"~",predexpr,"+",paste(mycovnames,collapse="+")))
        } else if(is.null(covix)){
          if(mice){
            ccvars<-unique(c(outnames[i],prednames[j],mycovnames_df,ixterm))
            if(is.null(micevars)){micevars<-ccvars[which(sapply(Data[,ccvars],function(x) any(is.na(x))))]}
            ccDat<-Data[complete.cases(Data[,ccvars[-which(ccvars%in%micevars)]]),]
            micedata<-mice(ccDat[,ccvars],m=miceiter,maxit=miceiter,printFlag=F)
          } else {
            ccDat<-Data[complete.cases(Data[,c(outnames[i],prednames[j],mycovnames_df,ixterm)]),]
          }
          nobs<-dim(ccDat)[1]
          form1<-formula(paste0(outexpr,"~",predexpr,"*",ixterm,"+",
                                paste(mycovnames,collapse="+")))
        } else {
          if(mice){
            ccvars<-unique(c(outnames[i],prednames[j],mycovnames_df,ixterm))
            if(is.null(micevars)){micevars<-ccvars[which(sapply(Data[,ccvars],function(x) any(is.na(x))))]}
            ccDat<-Data[complete.cases(Data[,ccvars[-which(ccvars%in%micevars)]]),]
            micedata<-mice(ccDat[,ccvars],m=miceiter,maxit=miceiter,printFlag=F)
          } else {
            ccDat<-Data[complete.cases(Data[,c(outnames[i],prednames[j],mycovnames_df,ixterm)]),]
          }
          nobs<-dim(ccDat)[1]
          covixaux<-rep("",length(mycovnames))
          covixaux[which(mycovnames%in%covix)]<-paste0("*",ixterm)
          for(cv in 1:length(mycovnames)){
            mycovnames[cv]<-paste0(mycovnames[cv],covixaux[cv])
          }
          if(ixpred==T){
            form1<-formula(paste0(outexpr,"~",predexpr,"*",ixterm,"+",
                                  paste(mycovnames,collapse="+")))
          } else {
            form1<-formula(paste0(outexpr,"~",predexpr,"+",
                                  paste(mycovnames,collapse="+")))
          }
        }
      }
      if(length(unique(Data[which(!is.na(Data[,outnames[i]])),outnames[i]]))==2){
        if(binfunc(binomlink)$family=="poisson"){
          if(is.character(Data[,outnames[i]])) Data[,outnames[i]]<-as.factor(Data[,outnames[i]])
          if(is.factor(Data[,outnames[i]])){
            levels(Data[,outnames[i]])<-c("0","1")
            Data[,outnames[i]]<-as.numeric(as.character(Data[,outnames[i]]))
          } 
        }
        if(Firth==T){
          lm1<-logistf(form1,data=Data,family=binomial(link=binomlink),na.action=na.exclude)
          if(leverage.test==T){
            warning("Unclear how to get Cook's distance for a 'logistf' object. Ignoring leverage.test.")
          }
        } else {
          if(mice){
            lm1<-list()
            tempdatalist<-list()
            for(mim in 1:micedata$m){
              tempdatalist[[mim]]<-ccDat
              for(mip in 1:length(micevars)){
                tempdatalist[[mim]][match(rownames(micedata$imp[[micevars[mip]]]),rownames(ccDat)),
                                    micevars[mip]]<-
                  micedata$imp[[micevars[mip]]][,mim]
              }; rm(mip)
              lm1[[mim]]<-glm(form1,data=tempdatalist[[mim]],family=binfunc(link=binomlink),na.action=na.exclude)
            }; rm(mim)
          } else {
            lm1<-glm(form1,data=Data,family=binfunc(link=binomlink),na.action=na.exclude)
          }
          if(leverage.test==T){
            if(mice){
              for(mim in 1:length(lm1)){
                cooks<-cooks.distance(lm1[[mim]])
                if(!is.null(leverage.cutoff)){
                  threshold<-leverage.cutoff
                } else if(!is.null(leverage.meancutoff)){
                  threshold<-mean(cooks,na.rm=T)*leverage.meancutoff
                } else {threshold<-mean(cooks,na.rm=T)*4}
                lm1<-glm(form1,data=Data[-c(as.numeric(names(cooks[which(cooks>threshold)]))),],
                         family=binfunc(link=binomlink),na.action=na.exclude)
              }
            } else {
              cooks<-cooks.distance(lm1)
              if(!is.null(leverage.cutoff)){
                threshold<-leverage.cutoff
              } else if(!is.null(leverage.meancutoff)){
                threshold<-mean(cooks,na.rm=T)*leverage.meancutoff
              } else {threshold<-mean(cooks,na.rm=T)*4}
              lm1<-glm(form1,data=Data[-c(as.numeric(names(cooks[which(cooks>threshold)]))),],
                       family=binfunc(link=binomlink),na.action=na.exclude)
            }
          }
          if(robust==T){
            if(mice){
              warning("Unclear how to get pooled robust SEs from imputed GLMs; ignoring robust command")
            } else {
              robustcoef<-robustse(lm1,HCtype=HCtype)
            }
          }
        }
        
      } else {
        if(mice){
          lm1<-list()
          tempdatalist<-list()
          for(mim in 1:micedata$m){
            tempdatalist[[mim]]<-ccDat
            for(mip in 1:length(micevars)){
              tempdatalist[[mim]][match(rownames(micedata$imp[[micevars[mip]]]),rownames(ccDat)),
                                  micevars[mip]]<-micedata$imp[[micevars[mip]]][,mim]
            }; rm(mip)
            lm1[[mim]]<-lm(form1,data=tempdatalist[[mim]],na.action=na.exclude)
          }; rm(mim)
        } else {
          lm1<-lm(form1,data=Data,na.action=na.exclude)
        }
        
        if(leverage.test==T){
          if(mice){
            for(mim in 1:length(lm1)){
              cooks<-cooks.distance(lm1[[mim]])
              if(!is.null(leverage.cutoff)){
                threshold<-leverage.cutoff
              } else if(!is.null(leverage.meancutoff)){
                threshold<-mean(cooks,na.rm=T)*leverage.meancutoff
              } else {threshold<-mean(cooks,na.rm=T)*4}
              lm1<-lm(form1,data=Data[-c(as.numeric(names(cooks[which(cooks>threshold)]))),],
                       na.action=na.exclude)
              
            }
          } else {
            cooks<-cooks.distance(lm1)
            if(!is.null(leverage.cutoff)){
              threshold<-leverage.cutoff
            } else if(!is.null(leverage.meancutoff)){
              threshold<-mean(cooks,na.rm=T)*leverage.meancutoff
            } else {
              threshold<-mean(cooks,na.rm=T)*4
            }
            lm1<-lm(form1,data=Data[-c(as.numeric(names(cooks[which(cooks>threshold)]))),],
                    na.action=na.exclude)
          }
        }
        if(robust==T){
          print(summary(lm1)) #to remove
          if(mice){
            warning("Unclear how to do robust SEs with MICE, so ignoring 'robust' command")
          } else {
            robustcoef<-robustse(lm1,HCtype=HCtype)
          }
        }
      }
      lmlist[[lmnum]]<-lm1
      names(lmlist)[[lmnum]]<-paste0(outnames[i],"~",prednames[j])
      Resultsmat[rownum,1]<-outnames[i]
      Resultsmat[rownum,3]<-nobs
      if(mice){
        Resultsmat[rownum,2]<-as.character(summary(pool(lm1))$term[2:(1+length(rownum))])
      } else if(Firth){
        Resultsmat[rownum,2]<-lm1$terms[2:(1+length(rownum))]
      } else {
        Resultsmat[rownum,2]<-dimnames(summary(lm1)$coef)[[1]][2:(1+length(rownum))]
      }
      if(!is.null(ixterm)){
        Resultsmat[rownum,2]<-prednames[j]
        if(mice) {
          micesumm<-summary(pool(lm1),conf.int=T)
          if(class(Data[,ixterm])=="factor"){
            comb1<-combn(levels(Data[,ixterm]),2)
            Resultsmat[rownum,12]<-c(
              paste(micesumm$term[2],levels(Data[,ixterm]),sep="|"),
              paste0(comb1[2,1:(nixlvl-1)]," - ",comb1[1,1:(nixlvl-1)]),
              paste0(micesumm$term[2],"*",comb1[2,]," - ",comb1[1,])
            )
            miceixsumm<-as.data.frame(matrix(NA,ixtermlength,6))
            names(miceixsumm)<-c("beta","SE","DF","pval","CIL","CIU")
            miceixsumm[1,]<-micesumm[2,c(2:3,5:8)]
            miceixsumm[(1+nixlvl:((2*nixlvl)-2)),]<-
              micesumm[c(3:(2+(nixlvl-1))),c(2:3,5:8)]
            miceixsumm[(2*nixlvl):
                         ((2*nixlvl)+((nixlvl)-2)),]<-
              micesumm[grep(paste0(prednames[j],":"),micesumm$term),c(2:3,5:8)]
            newData<-ccDat
            for(rp in 1:(nixlvl-1)){
              newData[,ixterm]<-factor(newData[,ixterm],
                                       levels=levels(newData[,ixterm])[
                                         c(2:(length(levels(newData[,ixterm]))),1)])
              newlm1<-list()
              tempdatalist<-list()
              for(mim in 1:micedata$m){
                tempdatalist[[mim]]<-newData
                for(mip in 1:length(micevars)){
                  tempdatalist[[mim]][match(rownames(micedata$imp[[micevars[mip]]]),rownames(newData)),
                                      micevars[mip]]<-micedata$imp[[micevars[mip]]][,mim]
                }; rm(mip)
                if(any(class(lm1[[1]])=="glm")){
                  newlm1[[mim]]<-glm(form1,data=tempdatalist[[mim]],family=binfunc(link=binomlink),
                                     na.action=na.exclude)
                } else {
                  newlm1[[mim]]<-lm(form1,data=tempdatalist[[mim]],na.action=na.exclude)
                }
              }; rm(mim)
              
              newmicesumm<-summary(pool(newlm1),conf.int=T)
              miceixsumm[1+rp,]<-newmicesumm[2,c(2:3,5:8)]
              
              nix2take<-(nixlvl-1)-rp
              if(nix2take>0){
                newfillrows<-c(
                  (((3*nixlvl)-1)+((nix2take+1)*(rp-1))):
                    ((((3*nixlvl)-1)+((nix2take+1)*(rp-1)))+(nix2take-1))
                )
                if(is.null(covix)){
                  newtakerows<-c((nrow(newmicesumm)-nixlvl+2):
                                   (nrow(newmicesumm)-(nixlvl-1)+nix2take))
                } else {
                  newimpixrows<-which(grepl(paste0(prednames[j],":"),newmicesumm$term))
                  newtakerows<-newimpixrows[1:nix2take]
                }
                
                miceixsumm[newfillrows,]<- newmicesumm[newtakerows,c(2:3,5:8)]
              }
            }; rm(rp)
            betalist<-miceixsumm$beta
            selist<-miceixsumm$SE
            dflist<-miceixsumm$DF
            pvallist<-miceixsumm$pval
            cilist<-miceixsumm[,c("CIL","CIU")]
          } else {
            betalist<-c(micesumm[c(2:3,grep(paste0(prednames[j],":"),micesumm$term)),2])
            pvallist<-c(micesumm[c(2:3,grep(paste0(prednames[j],":"),micesumm$term)),6])
            selist<-c(micesumm[c(2:3,grep(paste0(prednames[j],":"),micesumm$term)),3])
            dflist<-c(micesumm[c(2:3,grep(paste0(prednames[j],":"),micesumm$term)),5])
            cilist<-cbind(micesumm[c(2:3,grep(paste0(prednames[j],":"),micesumm$term)),7],
                          micesumm[c(2:3,grep(paste0(prednames[j],":"),micesumm$term)),8])
            Resultsmat[rownum,12]<-
              as.character(micesumm$term[c(2:3,grep(paste0(prednames[j],":"),
                                                          micesumm$term))])
          }
        } else if(any(class(lm1)=="logistf")){
          if(class(Data[,ixterm])=="factor"){
            comb1<-combn(levels(Data[,ixterm]),2)
            if(is.null(ncol(comb1))) comb1<-as.matrix(comb1,nrow=length(comb1))
            Resultsmat[rownum,12]<-c(
              paste(names(lm1$coefficients)[2],levels(Data[,ixterm]),sep="|"),
              paste0(comb1[2,1:(nixlvl-1)]," - ",comb1[1,1:(nixlvl-1)]),
              paste0(lm1$terms[2],"*",comb1[1,]," - ",comb1[2,])
            )
            betalist<-rep(NA,ixtermlength)
            betalist[1]<-lm1$coefficients[2]
            betalist[(1+nixlvl:((2*nixlvl)-2))]<-
              lm1$coefficients[c(3:(2+(nixlvl-1)))]
            betalist[(2*nixlvl):
                       ((2*nixlvl)+((nixlvl)-2))]<-
              -1*lm1$coefficients[c(
                (length(lm1$coefficients)-(nixlvl-2)):
                  length(lm1$coefficients)
              )]
            
            selist<-rep(NA,ixtermlength)
            selist[1]<-sqrt(diag(lm1$var))[2]
            selist[(1+nixlvl:((2*nixlvl)-2))]<-
              sqrt(diag(lm1$var))[c(3:(2+(nixlvl-1)))]
            selist[(2*nixlvl):
                     ((2*nixlvl)+((nixlvl)-2))]<-
              sqrt(diag(lm1$var))[c(
                (length(lm1$coefficients)-(nixlvl-2)):
                  length(lm1$coefficients)
              )]
            
            cilist<-matrix(NA,ixtermlength,2)
            cilist[1,]<-c(lm1$ci.lower[2],lm1$ci.upper[2])
            cilist[c(1+nixlvl:((2*nixlvl)-2)),]<-
              cbind(lm1$ci.lower[c(3:(2+(nixlvl-1)))],
                    lm1$ci.upper[c(3:(2+(nixlvl-1)))])
            
            cilist[c((2*nixlvl):
                       ((2*nixlvl)+((nixlvl)-2))),]<-
              cbind(-1*lm1$ci.upper[c(
                (length(lm1$coefficients)-(nixlvl-2)):
                  length(lm1$coefficients)
              )],
              -1*lm1$ci.lower[c(
                (length(lm1$coefficients)-(nixlvl-2)):
                  length(lm1$coefficients)
              )])
            
            pvallist<-rep(NA,ixtermlength)
            pvallist[1]<-lm1$prob[2]
            pvallist[(1+nixlvl):((2*nixlvl)-1)]<-
              lm1$prob[c(3:(2+(nixlvl-1)))]
            pvallist[(2*nixlvl):
                       ((2*nixlvl)+(nixlvl-2))]<-
              lm1$prob[c(
                (length(lm1$prob)-(nixlvl-2)):
                  length(lm1$prob)
              )]
            
            newData<-Data
            for(rp in 1:(nixlvl-1)){
              newData[,ixterm]<-factor(newData[,ixterm],
                                       levels=levels(newData[,ixterm])[
                                         c(2:(length(levels(newData[,ixterm]))),1)])
              newlm1<-logistf(form1,data=newData,family=binomial(link=binomlink),na.action=na.exclude)
              betalist[1+rp]<-newlm1$coefficients[2]
              selist[1+rp]<-sqrt(diag(newlm1$var))[2]
              pvallist[1+rp]<-newlm1$prob[2]
              cilist[1+rp,]<-cbind(newlm1$ci.lower[2],
                                   newlm1$ci.upper[2])
              
              nix2take<-(nixlvl-1)-rp
              if(nix2take>0){
                newtakerows<-c((length(lm1$coefficients)-nixlvl+2):
                                 (length(lm1$coefficients)-(nixlvl-1)+nix2take))
                newfillrows<-c(
                  (((3*nixlvl)-1)+((nix2take+1)*(rp-1))):
                    ((((3*nixlvl)-1)+((nix2take+1)*(rp-1)))+(nix2take-1))
                )
                
                betalist[newfillrows]<- -1*newlm1$coefficients[newtakerows]
                selist[newfillrows]<-sqrt(diag(newlm1$var))[newtakerows]
                pvallist[newfillrows]<-newlm1$prob[newtakerows]
                cilist[newfillrows,]<-cbind(-1*newlm1$ci.upper[newtakerows],
                                            -1*newlm1$ci.lower[newtakerows])
              }
            }
          } else {
            betalist<-c(lm1$coefficients[c(2:3,length(lm1$coefficients))])
            pvallist<-c(lm1$prob[c(2:3,length(lm1$prob))])
            selist<-c(sqrt(diag(lm1$var))[c(2:3,nrow(lm1$var))])
            cilist<-cbind(lm1$ci.lower[c(2:3,length(lm1$prob))],
                          lm1$ci.upper[c(2:3,length(lm1$prob))])
          }
        } else if (robust==T) {
          if(class(Data[,ixterm])=="factor"){
            comb1<-combn(levels(Data[,ixterm]),2)
            Resultsmat[rownum,12]<-c(
              paste(dimnames(robustcoef)[[1]][2],levels(Data[,ixterm]),sep="|"),
              paste0(comb1[2,1:(nixlvl-1)]," - ",comb1[1,1:(nixlvl-1)]),
              paste0(dimnames(robustcoef)[[1]][2],"*",comb1[2,]," - ",comb1[1,])
            )
            betalist<-rep(NA,ixtermlength)
            betalist[1]<-robustcoef[2,1]
            betalist[(1+nixlvl:((2*nixlvl)-2))]<-
              robustcoef[c(3:(2+(nixlvl-1))),1]
            betalist[(2*nixlvl):
                       ((2*nixlvl)+((nixlvl)-2))]<-
              robustcoef[grep(paste0(prednames[j],":"),dimnames(robustcoef)[[1]]),1]#-1*
            #c((nrow(robustcoef)-(nixlvl-2)):nrow(robustcoef))
            selist<-rep(NA,ixtermlength)
            selist[1]<-robustcoef[2,2]
            selist[(1+nixlvl:((2*nixlvl)-2))]<-
              robustcoef[c(3:(2+(nixlvl-1))),2]
            selist[(2*nixlvl):
                     ((2*nixlvl)+((nixlvl)-2))]<-
              robustcoef[grep(paste0(prednames[j],":"),dimnames(robustcoef)[[1]]),2]
            #c((nrow(robustcoef)-(nixlvl-2)):nrow(robustcoef))
            cilist<-matrix(NA,ixtermlength,2)
            cilist[1,]<-c(robustcoef[2,5],robustcoef[2,6])
            cilist[1+nixlvl:((2*nixlvl)-2),]<-
              cbind(robustcoef[c(3:(2+(nixlvl-1))),5],
                    robustcoef[c(3:(2+(nixlvl-1))),6])
            
            cilist[c((2*nixlvl):((2*nixlvl)+((nixlvl)-2))),]<-
              cbind(robustcoef[grep(paste0(prednames[j],":"),dimnames(robustcoef)[[1]]),5],#-1*
                    robustcoef[grep(paste0(prednames[j],":"),dimnames(robustcoef)[[1]]),6])#-1*
            #c((nrow(robustcoef)-(nixlvl-2)):nrow(robustcoef))
            
            pvallist<-rep(NA,ixtermlength)
            pvallist[1]<-robustcoef[2,4]
            pvallist[(1+nixlvl:((2*nixlvl)-2))]<-
              robustcoef[c(3:(2+(nixlvl-1))),4]
            pvallist[(2*nixlvl):((2*nixlvl)+(nixlvl-2))]<-
              robustcoef[grep(paste0(prednames[j],":"),dimnames(robustcoef)[[1]]),4]
            #c((nrow(robustcoef)-(nixlvl-2)):nrow(robustcoef))
            
            newData<-Data
            for(rp in 1:(nixlvl-1)){
              newData[,ixterm]<-factor(newData[,ixterm],
                                       levels=levels(newData[,ixterm])[
                                         c(2:(length(levels(newData[,ixterm]))),1)])
              if(any(class(lm1)=="glm")){
                newlm1<-glm(form1,data=newData,
                            family=binfunc(link=binomlink),na.action=na.exclude)
              } else {
                newlm1<-lm(form1,data=newData,na.action=na.exclude)
              }
              
              newrobustcoef<-robustse(newlm1,HCtype=HCtype)
              betalist[1+rp]<-newrobustcoef[2,1]
              selist[1+rp]<-newrobustcoef[2,2]
              pvallist[1+rp]<-newrobustcoef[2,4]
              cilist[1+rp,]<-cbind(newrobustcoef[2,5],
                                   newrobustcoef[2,6])
              
              nix2take<-(nixlvl-1)-rp
              if(nix2take>0){
                newfillrows<-c(
                  (((3*nixlvl)-1)+((nix2take+1)*(rp-1))):
                    ((((3*nixlvl)-1)+((nix2take+1)*(rp-1)))+(nix2take-1))
                )
                if(is.null(covix)){
                  newtakerows<-c((nrow(newrobustcoef)-nixlvl+2):
                                   (nrow(newrobustcoef)-(nixlvl-1)+nix2take))
                } else {
                  newimpixrows<-which(grepl(paste0(prednames[j],":"),rownames(newrobustcoef)))
                  newtakerows<-newimpixrows[1:nix2take]
                }
                
                betalist[newfillrows]<- newrobustcoef[newtakerows,1]#-1*
                selist[newfillrows]<-newrobustcoef[newtakerows,2]
                pvallist[newfillrows]<-newrobustcoef[newtakerows,4]
                cilist[newfillrows,]<-cbind(newrobustcoef[newtakerows,5],#-1* ,6]
                                            newrobustcoef[newtakerows,6])#-1* ,5]
              }
            }
          } else {
            betalist<-c(robustcoef[c(2:3,grep(paste0(prednames[j],":"),dimnames(robustcoef)[[1]])),1])
            pvallist<-c(robustcoef[c(2:3,grep(paste0(prednames[j],":"),dimnames(robustcoef)[[1]])),4])
            selist<-c(robustcoef[c(2:3,grep(paste0(prednames[j],":"),dimnames(robustcoef)[[1]])),2])
            cilist<-cbind(robustcoef[c(2:3,grep(paste0(prednames[j],":"),dimnames(robustcoef)[[1]])),5],
                          robustcoef[c(2:3,grep(paste0(prednames[j],":"),dimnames(robustcoef)[[1]])),6])
            Resultsmat[rownum,12]<-
              dimnames(summary(lm1)$coef)[[1]][c(2:3,grep(paste0(prednames[j],":"),
                                                          dimnames(robustcoef)[[1]]))]
          }
        } else {
          if(class(Data[,ixterm])=="factor"){
            comb1<-combn(levels(Data[,ixterm]),2)
            Resultsmat[rownum,12]<-c(
              paste(rownames(summary(lm1)$coef)[2],levels(Data[,ixterm]),sep="|"),
              paste0(comb1[2,1:(nixlvl-1)]," - ",comb1[1,1:(nixlvl-1)]),
              paste0(rownames(summary(lm1)$coef)[2],"*",comb1[2,]," - ",comb1[1,])
            )
            betalist<-rep(NA,ixtermlength)
            betalist[1]<-lm1$coefficients[2]
            betalist[(1+nixlvl:((2*nixlvl)-2))]<-
              lm1$coefficients[c(3:(2+(nixlvl-1)))]
            betalist[(2*nixlvl):
                       ((2*nixlvl)+((nixlvl)-2))]<-
              lm1$coefficients[c( #-1*
                (length(lm1$coefficients)-(nixlvl-2)):
                  length(lm1$coefficients)
              )]
            
            selist<-rep(NA,ixtermlength)
            selist[1]<-summary(lm1)$coef[2,2]
            selist[(1+nixlvl:((2*nixlvl)-2))]<-
              summary(lm1)$coef[c(3:(2+(nixlvl-1))),2]
            selist[(2*nixlvl):
                     ((2*nixlvl)+((nixlvl)-2))]<-
              summary(lm1)$coef[c(
                (length(lm1$coefficients)-(nixlvl-2)):
                  length(lm1$coefficients)
              ),2]
            
            cilist<-matrix(NA,ixtermlength,2)
            lm1.ci<-suppressMessages(confint(lm1))
            cilist[1,]<-lm1.ci[2,]
            cilist[c(1+nixlvl:((2*nixlvl)-2)),]<-
              lm1.ci[c(3:(2+(nixlvl-1))),]
            cilist[c((2*nixlvl):
                       ((2*nixlvl)+((nixlvl)-2))),]<-
              lm1.ci[c(
                (length(lm1$coefficients)-(nixlvl-2)):
                  length(lm1$coefficients)
              ),]
            
            pvallist<-rep(NA,ixtermlength)
            pvallist[1]<-summary(lm1)$coef[2,4]
            pvallist[(1+nixlvl:((2*nixlvl)-2))]<-
              summary(lm1)$coef[c(3:(2+(nixlvl-1))),4]
            pvallist[(2*nixlvl):
                     ((2*nixlvl)+((nixlvl)-2))]<-
              summary(lm1)$coef[c(
                (length(lm1$coefficients)-(nixlvl-2)):
                  length(lm1$coefficients)
              ),4]
            
            newData<-Data
            for(rp in 1:(nixlvl-1)){
              newData[,ixterm]<-factor(newData[,ixterm],
                                       levels=levels(newData[,ixterm])[
                                         c(2:(length(levels(newData[,ixterm]))),1)])
              if(any(class(lm1)=="glm")){
                newlm1<-glm(form1,data=newData,family=binfunc(link=binomlink),na.action=na.exclude)
              } else {
                newlm1<-lm(form1,data=newData,na.action=na.exclude)
              }
              
              betalist[1+rp]<-newlm1$coefficients[2]
              selist[1+rp]<-summary(newlm1)$coef[2,2]
              pvallist[1+rp]<-summary(newlm1)$coef[2,4]
              cilist[1+rp,]<-suppressMessages(confint(newlm1))[2,]
              
              nix2take<-(nixlvl-1)-rp
              if(nix2take>0){
                newtakerows<-c((length(lm1$coefficients)-nixlvl+2):
                                 (length(lm1$coefficients)-(nixlvl-1)+nix2take))
                newfillrows<-c(
                  (((3*nixlvl)-1)+((nix2take+1)*(rp-1))):
                    ((((3*nixlvl)-1)+((nix2take+1)*(rp-1)))+(nix2take-1))
                )
                
                betalist[newfillrows]<- newlm1$coefficients[newtakerows]
                selist[newfillrows]<-summary(newlm1)$coef[newtakerows,2]
                pvallist[newfillrows]<-summary(newlm1)$coef[newtakerows,4]
                cilist[newfillrows,]<-suppressMessages(confint(newlm1))[newtakerows,]
              }
            }
          } else {
            Resultsmat[rownum,12]<-
              dimnames(summary(lm1)$coef)[[1]][c(2:3,grep(paste0(prednames[j],":"),
                                                          dimnames(summary(lm1)$coef)[[1]]))]
          }
        }
        
      } else if(predspline==T){
        Resultsmat[rownum,12]<-dimnames(summary(lm1)$coef)[[1]][2:4]
      }
      
      if(MyMult=="IQR"){
        if(logpred2[j]==TRUE){
          Mult<-as.numeric((IQR(Data[,prednames[j]],na.rm=T)/
                              summary(Data[,prednames[j]])[2])*100)
        } else {
          Mult<-IQR(Data[,prednames[j]],na.rm=T)
        }
      } else {
        Mult<-MyMult
      }
      
      if(!is.null(ixterm)){
        if(mice){
          if(any(class(lm1[[1]])=="glm")&binomlink!="probit"){
            Resultsmat[rownum,4]<-exp(betalist)
            Resultsmat[rownum,5:6]<-exp(cilist)
            Resultsmat[rownum,7]<-pvallist
          } else {
            Resultsmat[rownum,4]<-betalist
            Resultsmat[rownum,5:6]<-cilist
            Resultsmat[rownum,7]<-pvallist
          }
        } else if(any(class(lm1)=="logistf")){
          Resultsmat[rownum,4]<-exp(betalist)
          Resultsmat[rownum,5:6]<-exp(cilist)
          Resultsmat[rownum,7]<-pvallist
        } else if(robust==T) {
          if(any(class(lm1)=="glm")&binomlink!="probit"){
            Resultsmat[rownum,4]<-exp(betalist)
            Resultsmat[rownum,5:6]<-exp(cilist)
            Resultsmat[rownum,7]<-pvallist
          } else {
            Resultsmat[rownum,4]<-betalist
            Resultsmat[rownum,5:6]<-cilist
            Resultsmat[rownum,7]<-pvallist
          }
        } else {
          if(class(Data[,ixterm])=="factor"){
            if(any(class(lm1)=="glm")&binomlink!="probit"){
              Resultsmat[rownum,4]<-exp(betalist)
              Resultsmat[rownum,5:6]<-exp(cilist)
              Resultsmat[rownum,7]<-c(pvallist)
            } else {
              Resultsmat[rownum,4]<-betalist
              Resultsmat[rownum,5:6]<-cilist
              Resultsmat[rownum,7]<-pvallist
            }
          } else {
            mysumrows<-c(2:3,nrow(summary(lm1)$coef))
            betalist<-summary(lm1)$coef[mysumrows,1]
            selist<-summary(lm1)$coef[mysumrows,2]
            cilist<-suppressMessages(confint(lm1))[mysumrows,1:2]
            pvallist<-summary(lm1)$coef[mysumrows,4]
            if(any(class(lm1)=="glm")){
              Resultsmat[rownum,4:7]<-c(exp(betalist),exp(cilist),pvallist)
            } else {
              Resultsmat[rownum,4:7]<-c(betalist,cilist,pvallist)
            }
          }
        }
        
      } else {
        if(predspline==T){
          mysumrows<-c(2:4)
        } else if(class(Data[,prednames[j]])=="factor"){
          mysumrows<-c(2:(1+length(rownum)))
        } else {
          mysumrows<-2
        }
        
        if(mice){
          if(length(unique(Data[!is.na(Data[,outnames[i]]),outnames[i]]))==2){
            micesumm<-summary(pool(lm1),conf.int=T,exponentiate=T)
          } else {
            micesumm<-summary(pool(lm1),conf.int=T)
          }
          betalist<-micesumm[mysumrows,2]
          selist<-micesumm[mysumrows,3]
          cilist<-micesumm[mysumrows,7:8]
          dflist<-micesumm[mysumrows,"df"]
          pvallist<-micesumm[mysumrows,6]
        } else if(any(class(lm1)=="logistf")){
          betalist<-lm1$coefficients[2]
          selist<-sqrt(diag(lm1$var))[2]
          cilist<-cbind(lm1$ci.lower[2],lm1$ci.upper[2])
          pvallist<-lm1$prob[2]
        } else if(robust==T){
          betalist<-robustcoef[mysumrows,1] #used to be row = "2"
          selist<-robustcoef[mysumrows,2]
          cilist<-robustcoef[mysumrows,5:6]
          pvallist<-robustcoef[mysumrows,4]
        } else {
          betalist<-summary(lm1)$coef[mysumrows,1]
          selist<-summary(lm1)$coef[mysumrows,2]
          cilist<-confint.lm(lm1)[mysumrows,1:2]
          pvallist<-summary(lm1)$coef[mysumrows,4]
        }
        
        
        if(length(unique(Data[!is.na(Data[,outnames[i]]),outnames[i]]))==2&!mice){
          Resultsmat[rownum,4:7]<-c(exp(betalist),exp(cilist),pvallist)
        } else {
          Resultsmat[rownum,4:7]<-c(betalist,cilist,pvallist)
        }
      }
      
      if(mice){
        if(logpred2[j]==T){
          IQRBeta<-betalist*log(1+(Mult/100))
          IQRCIs<-confint.qt(betalist,selist,dflist,log(1+(Mult/100)))
        } else {
          IQRBeta<-betalist*Mult
          IQRCIs<-confint.qt(betalist,selist,dflist,Mult)
        }
      } else if(Firth==T|any(class(lm1)=="glm")&robust==F){
        if(Mult!=1){
          warning("The multiplying factor for beta was ignored because it's unclear how to get profile 
                  penalized log likelihood CIs for logistf or profile log likelihood CIs for glm for a 
                  change in predictor not equal to 1.")
        }
        IQRBeta<-betalist
        IQRCIs<-cilist
      } else {
        if(logpred2[j]==T){
          IQRBeta<-betalist*log(1+(Mult/100))
          IQRCIs<-confint.qt(betalist,selist,summary(lm1)$df[2],log(1+(Mult/100)))
        } else {
          IQRBeta<-betalist*Mult
          IQRCIs<-confint.qt(betalist,selist,summary(lm1)$df[2],Mult)
        }
      }
      
      if(logout2[i]==T|length(unique(Data[!is.na(Data[,outnames[i]]),outnames[i]]))==2){
        Resultsmat[rownum,8]<-(exp(IQRBeta)-1)*100
        Resultsmat[rownum,9:10]<-c((exp(IQRCIs)-1)*100)
      } else {
        Resultsmat[rownum,8]<-IQRBeta
        Resultsmat[rownum,9:10]<-IQRCIs
      }
      if(LOOCV==T){
        if(length(unique(Data[!is.na(Data[,outnames[i]]),outnames[i]]))!=2){
          t1<-train(form1,ccDat,trControl = trainControl(method="repeatedcv",number=10,
                                                         repeats=50),method="lm")
          Resultsmat[rownum,11]<-t1$results$RMSE
        } else if(any(class(lm1)=="glm")){
          t1<-train(formula(paste0("factor(",outexpr,")~",as.character(form1)[3])),ccDat,
                    trControl = trainControl(method="repeatedcv",number=10,repeats=50),
                    method="glm",family=binfunc(link=binomlink))
          Resultsmat[rownum,11]<-t1$results$Accuracy
        } else {
          Resultsmat[rownum,11]<-NA
        }
      } else {
        Resultsmat[rownum,11]<-NA
      }
      if(mice){
        Resultsmat[rownum,"AIC"]<-NA
      } else if(Firth){
        extractAIC.logistf<-function(fit, scale, k=2, ...){
          dev<- -2 * (fit$loglik['null']-fit$loglik['full'])
          AIC<- dev+k*fit$df
          edf<- fit$df
          return(AIC)
        }
        Resultsmat[rownum,"AIC"]<-extractAIC.logistf(lmlist[[lmnum]])
      } else {
        Resultsmat[rownum,"AIC"]<-AIC(lmlist[[lmnum]])
      }
      if(extradiag){
        if(mice){
          warning("Extra diagnostics not attainable for imputed models; ignoring extradiag")
        } else if(!any(class(lmlist[[lmnum]])=="glm")){
          cooks<-cooks.distance(lmlist[[lmnum]])
          Extramat[rownum,"Heterosk.p"]<-lmtest::bptest(lmlist[[lmnum]])$p.value
          Extramat[rownum,"N.Cooks.D_>_0.5"]<-length(which(cooks>0.5))
          Extramat[rownum,"MaxCooks.D"]<-max(cooks,na.rm=T)
        } else {
          cooks<-cooks.distance(lmlist[[lmnum]])
          # auxDat<-lmlist[[lmnum]]$model
          # if(any(grepl("(",names(auxDat),fixed=T))){names(auxDat)<-gsub(".*\\(|\\,.*|\\).*","",names(auxDat))}
          # print(names(auxDat))
          # print(form1)
          altmod<-glm(form1,data=lmlist[[lmnum]]$model,family=lmlist[[lmnum]]$family)
          simout<-DHARMa::simulateResiduals(altmod)
          quanttest<-DHARMa::testQuantiles(simout,plot=F)
          outliertest<-DHARMa::testOutliers(simout,plot=F,type='bootstrap')
          uniftest<-DHARMa::testUniformity(simout,plot=F)
          disptest<-DHARMa::testDispersion(simout,plot=F)
          Extramat[rownum,"Heterosk.p"]<-quanttest$p.value
          Extramat[rownum,"MinQuantileSplinePval"]<-min(quanttest$pvals)
          Extramat[rownum,"N.Cooks.D_>_0.5"]<-length(which(cooks>0.5))
          Extramat[rownum,"MaxCooks.D"]<-max(cooks,na.rm=T)
          Extramat[rownum,"UnifTestPval"]<-uniftest$p.value
          Extramat[rownum,"OutlierTestPval"]<-outliertest$p.value
          Extramat[rownum,"DispTestPval"]<-disptest$p.value
        }
      }
      if(post.power==T){
        suppressPackageStartupMessages(library(SimMultiCorrData))
        suppressPackageStartupMessages(library(pbapply))
        Xmat<-model.matrix(form1,ccDat)[,-1]
        Xmatbinary<-which(apply(Xmat,2,function(x) setequal(unique(x),c(0,1)))==T)
        Xmatcont<-which(apply(Xmat,2,function(x) setequal(unique(x),c(0,1)))==F)
        Xmat<-Xmat[,c(Xmatcont,Xmatbinary)]
        Xmatbinary<-which(apply(Xmat,2,function(x) setequal(unique(x),c(0,1)))==T)
        Xmatcont<-which(apply(Xmat,2,function(x) setequal(unique(x),c(0,1)))==F)
        M<-apply(Xmat,2,calc_moments)
        catmarginals<-as.list(1-M[1,Xmatbinary])
        support <- list() # default support will be generated inside simulation
        seeds<-sample(1:(nsim*10),nsim)
        powersimfunc<-function(seed=seeds){
          if(is.null(ixterm)){
            myvc<-cov2cor(cov(Xmat))
            invisible(
              Simmat<-rcorrvar(nrow(Xmat),length(Xmatcont),length(Xmatbinary),means = M[1,Xmatcont], 
                               vars =  (M[2,Xmatcont])^2, skews = M[3,Xmatcont], skurts = M[4,Xmatcont], 
                               fifths = M[5,Xmatcont], sixths = M[6,Xmatcont], marginal = catmarginals,
                               rho = myvc, seed = seed)
            )
            newX<-cbind(Simmat$continuous_variables,Simmat$ordinal_variables-1)
            names(newX)<-dimnames(Xmat)[[2]]
          } else {
            myvc<-cov2cor(cov(Xmat[,-length(Xmatcont)]))
            invisible(
              Simmat<-rcorrvar(nrow(Xmat),length(Xmatcont)-1,length(Xmatbinary),
                               means = M[1,Xmatcont[-length(Xmatcont)]], 
                               vars =  (M[2,Xmatcont[-length(Xmatcont)]])^2, 
                               skews = M[3,Xmatcont[-length(Xmatcont)]], 
                               skurts = M[4,Xmatcont[-length(Xmatcont)]], 
                               fifths = M[5,Xmatcont[-length(Xmatcont)]], 
                               sixths = M[6,Xmatcont[-length(Xmatcont)]], marginal = catmarginals,
                               rho = myvc, seed = seed)
            )
            newX<-cbind(Simmat$continuous_variables,Simmat$ordinal_variables-1)
            names(newX)<-dimnames(Xmat)[[2]][-length(Xmatcont)]
            newX[,dimnames(Xmat)[[2]][length(Xmatcont)]]<-
              newX[,dimnames(Xmat)[[2]][1]]*newX[,dimnames(Xmat)[[2]][length(Xmatcont)+1]]
            newX<-newX[,c(1:(length(Xmatcont)-1),ncol(newX),length(Xmatcont):(ncol(newX)-1))]
          }
          betas<-lm1$coef[names(newX)]+sapply(summary(lm1)$coef[names(newX),2],
                                              function(x) rnorm(1,0,x)) 
          #keep general beta structure but introduce noise with SD = beta SE
          if(is.null(ixterm)){
            betas[1]<-effect.size
          } else {
            betas[length(betas)]<-effect.size
          }
          my.y<-as.matrix(newX)%*%betas+lm1$coef[1]+rnorm(nrow(newX),sd(lm1$residuals)) #keep former residual SD
          newDat<-cbind(my.y,newX)
          names(newDat)[1]<-outnames[i]
          names(newDat)<-gsub(" ","_",names(newDat))
          names(newDat)<-gsub("$","",names(newDat),fixed=T)
          names(newDat)<-gsub(":","_",names(newDat),fixed=T)
          names(newDat)<-gsub("+","_",names(newDat),fixed=T)
          names(newDat)<-gsub("-","_",names(newDat),fixed=T)
          names(newDat)<-gsub("<","_",names(newDat),fixed=T)
          names(newDat)<-gsub(">","_",names(newDat),fixed=T)
          newlm<-lm(formula(paste0(outnames[i],"~",paste(names(newDat)[-1],collapse="+"))),newDat)
          if(robust==T){
            newcoefs<-robustse(newlm)
            rets<-newcoefs[2,c(4)]
          } else {
            rets<-summary(newlm)$coef[2,c(4)]
          }
          return(rets)
        }
        pboptions(type="timer")
        pvalvec<-pbsapply(seeds,powersimfunc)
        pow<-length(which(pvalvec<0.05))/length(pvalvec)
        Resultsmat[rownum,"Posthoc_Power"]<-pow
      }
    } #end j loop
  } #end i loop
  if(post.power==T){
    if(!is.null(ixterm)|predspline==T){
      if(any(logout)){
        names(Resultsmat)<-c(Outtitle,Predtitle,"N","Beta","LCI","UCI","p-value",
                             "PercBeta.Mult","PercLCI.Mult","PercUCI.Mult",
                             "K-fold RMSE","Variable","AIC","Posthoc_Power")
      } else if(length(unique(Data[!is.na(Data[,outnames[i]]),outnames[i]]))==2){
        names(Resultsmat)<-c(Outtitle,Predtitle,"N","OR","ORLCI","ORUCI","p-value",
                             "OR.Mult","ORLCI.Mult","ORUCI.Mult",
                             "K-fold Accuracy","Variable","AIC","Posthoc_Power")
      } else {
        names(Resultsmat)<-c(Outtitle,Predtitle,"N","Beta","LCI","UCI","p-value",
                             "Beta.Mult","LCI.Mult","UCI.Mult",
                             "K-fold RMSE","Variable","AIC","Posthoc_Power")
      }
      
    } else {
      if(any(logout)){
        names(Resultsmat)<-c(Outtitle,Predtitle,"N","Beta","LCI","UCI","p-value",
                             "PercBeta.Mult","PercLCI.Mult","PercUCI.Mult",
                             "K-fold RMSE","AIC","Posthoc_Power")
      } else if(length(unique(Data[!is.na(Data[,outnames[i]]),outnames[i]]))==2){
        names(Resultsmat)<-c(Outtitle,Predtitle,"N","OR","ORLCI","ORUCI","p-value",
                             "OR.Mult","ORLCI.Mult","ORUCI.Mult",
                             "K-fold Accuracy","AIC","Posthoc_Power")
      } else {
        names(Resultsmat)<-c(Outtitle,Predtitle,"N","Beta","LCI","UCI","p-value",
                             "Beta.Mult","LCI.Mult","UCI.Mult",
                             "K-fold RMSE","AIC","Posthoc_Power")
      }
    }
  } else {
    if(!is.null(ixterm)|predspline==T){
      if(any(logout)){
        names(Resultsmat)<-c(Outtitle,Predtitle,"N","Beta","LCI","UCI","p-value",
                             "PercBeta.Mult","PercLCI.Mult","PercUCI.Mult",
                             "K-fold RMSE","Variable","AIC")
      } else if(length(unique(Data[!is.na(Data[,outnames[i]]),outnames[i]]))==2){
        names(Resultsmat)<-c(Outtitle,Predtitle,"N","OR","ORLCI","ORUCI","p-value",
                             "OR.Mult","ORLCI.Mult","ORUCI.Mult",
                             "K-fold Accuracy","Variable","AIC")
      } else {
        names(Resultsmat)<-c(Outtitle,Predtitle,"N","Beta","LCI","UCI","p-value",
                             "Beta.Mult","LCI.Mult","UCI.Mult",
                             "K-fold RMSE","Variable","AIC")
      }
      
    } else {
      if(any(logout)){
        names(Resultsmat)<-c(Outtitle,Predtitle,"N","Beta","LCI","UCI","p-value",
                             "PercBeta.Mult","PercLCI.Mult","PercUCI.Mult",
                             "K-fold RMSE","AIC")
      } else if(length(unique(Data[!is.na(Data[,outnames[i]]),outnames[i]]))==2){
        names(Resultsmat)<-c(Outtitle,Predtitle,"N","OR","ORLCI","ORUCI","p-value",
                             "OR.Mult","ORLCI.Mult","ORUCI.Mult",
                             "K-fold Accuracy","AIC")
      } else {
        names(Resultsmat)<-c(Outtitle,Predtitle,"N","Beta","LCI","UCI","p-value",
                             "Beta.Mult","LCI.Mult","UCI.Mult",
                             "K-fold RMSE","AIC")
      }
    }
  }
  
  if(!is.null(ixterm)|predspline==T){
    if(post.power==T){
      Resultsmat<-Resultsmat[,c(1:2,12,3:11,13:14)]
    } else {
      Resultsmat<-Resultsmat[,c(1:2,12,3:11,13)]
    }
  }
  if(LOOCV==F){
    Resultsmat$'K-fold RMSE'<-NULL
    Resultsmat$'K-fold Accuracy'<-NULL
  }
  Resultsmat$Significant<-ifelse(Resultsmat$'p-value'<0.05,"Yes","")
  Resultsmat<-Resultsmat[which(apply(Resultsmat,1,function(x) !all(is.na(x)))),]
  if(extradiag){
    if(length(unique(Data[!is.na(Data[,outnames[i]]),outnames[i]]))==2){
      names(Extramat)[which(names(Extramat)=="Heterosk.p")]<-"QuantileSplinePval"
    }
    Resultsmat<-cbind(Resultsmat,Extramat)
  }
  Resmatout<-Resultsmat
  
  #Resmatout$Significant<-NULL #kable won't let me get rid of the Significant column, IDK why
  names(Resultsmat)[2]<-"Predictor"
  
  if(length(unique(Data[!is.na(Data[,outnames[i]]),outnames[i]]))==2){
    names(Resultsmat)[names(Resultsmat)%in%c("OR.Mult","ORLCI.Mult","ORUCI.Mult")]<-
      c("PercBeta.Mult","PercLCI.Mult","PercUCI.Mult")
    names(Resultsmat)[names(Resultsmat)%in%c("OR","ORLCI","ORUCI")]<-
      c("Beta","LCI","UCI")
  }
  
  if(!is.null(ixterm)|predspline==T){
    k1<-kable(Resmatout[,which(!grepl("Mult",names(Resmatout)))],"html") %>%
      kable_styling(font_size=14,full_width = F) %>%
      column_spec(1,bold=T) %>%
      collapse_rows(columns=c(1:2,9),valign="top")
  } else {
    k1<-kable(Resmatout[,which(!grepl("Mult",names(Resmatout)))],"html") %>%
      kable_styling(font_size=14,full_width = F) %>%
      column_spec(1,bold=T) %>%
      collapse_rows(columns=c(1:2,8),valign="top")
  }
  
  if(MyMult=="IQR"){
    Multname<-"IQR"
  } else {
    Multname<-Mult
  }
  
  if(is.null(facetcol)){
    facetcol<-length(outnames)
  }
  
  if(!is.null(ixterm)&class(Data[,ixterm])=="factor"){
    plotData<-Resultsmat[grep("|",Resultsmat$Variable,fixed=T),]
    plotData$Interaction_Level<-gsub(".*[|]","",plotData$Variable)
    # plotData$Interaction_Level<-factor(plotData$Interaction_Level,
    #                                       levels=levels(Data[,ixterm]))
    plotData$Predictor<-gsub("log10|_SGadj","",plotData$Predictor)
    plotData$Predictor<-factor(plotData$Predictor,levels=unique(plotData$Predictor))
    plotData[,Outtitle]<-factor(plotData[,Outtitle],levels=outnames)
    plotData$star<-ifelse(plotData$Significant=="Yes","*","")
    plotData.ix<-Resultsmat[grep("*",Resultsmat$Variable,fixed=T),]
    plotData.ix$Interaction_Level<-gsub(".*[*]","",plotData.ix$Variable)
    plotData.ix$Predictor<-gsub("log10|_SGadj","",plotData.ix$Predictor)
    plotData.ix$Predictor<-factor(plotData.ix$Predictor,levels=unique(plotData.ix$Predictor))
    plotData.ix$star<-ifelse(plotData.ix$Significant=="Yes","*","")
    plotData.ix[,Outtitle]<-factor(plotData.ix[,Outtitle],levels=outnames)
  } else if(!is.null(ixterm)) {
    plotData<-Resultsmat[-grep(":",Resultsmat$Variable,fixed=T),]
    plotData$Interaction_Level<-plotData$Variable
    plotData$star<-ifelse(plotData$Significant=="Yes","*","")
    plotData$Predictor<-gsub("log10|_SGadj","",plotData$Predictor)
    plotData$Predictor<-factor(plotData$Predictor,levels=unique(plotData$Predictor))
    plotData[,Outtitle]<-factor(plotData[,Outtitle],levels=outnames)
    plotData.ix<-Resultsmat[grep(":",Resultsmat$Variable,fixed=T),]
    plotData.ix$Interaction_Level<-plotData.ix$Variable
    plotData.ix$star<-ifelse(plotData.ix$Significant=="Yes","*","")
    plotData.ix$Predictor<-gsub("log10|_SGadj","",plotData.ix$Predictor)
    plotData.ix$Predictor<-factor(plotData.ix$Predictor,levels=unique(plotData.ix$Predictor))
    plotData.ix[,Outtitle]<-factor(plotData.ix[,Outtitle],levels=outnames)
  } else {
    plotData<-Resultsmat
    plotData$star<-ifelse(plotData$Significant=="Yes","*","")
    plotData$Predictor<-gsub("log10|_SGadj","",plotData$Predictor)
    plotData$Predictor<-factor(plotData$Predictor,levels=unique(plotData$Predictor))
    plotData[,Outtitle]<-factor(plotData[,Outtitle],levels=outnames)
  }
  
  
  if(predspline==F){
    if(any(logout)==T|length(unique(Data[!is.na(Data[,outnames[i]]),outnames[i]]))==2){
      if(is.null(ixterm)){
        if(plotOR==T){
          gg1<-ggplot(data=plotData,aes(x=Predictor,y=Beta,
                                        color=Significant))+
            geom_errorbar(aes(ymin=LCI,ymax=UCI),width=0.4,size=1)+
            geom_point()+geom_hline(aes(yintercept=1))+
            geom_text(aes(y=UCI+abs(UCI/20),
                          label=star),position=position_dodge(0.8),show.legend=F,size=6)+
            facet_wrap(formula(paste0(".~",Outtitle)),scales="free",ncol=facetcol)+
            ylab("OR")+xlab(paste0(Predtitle))+
            scale_color_manual(values=c("black","red"),guide="none")+
            theme(axis.text.x=element_text(angle=45,hjust=1))
        } else {
          if(!plotPercChange) {
            plotData$PercBeta.Mult<-log(plotData$Beta)
            plotData$PercLCI.Mult<-log(plotData$LCI)
            plotData$PercUCI.Mult<-log(plotData$UCI)
          }
          gg1<-ggplot(data=plotData,aes(x=Predictor,y=PercBeta.Mult,
                                        color=Significant))+
            geom_errorbar(aes(ymin=PercLCI.Mult,ymax=PercUCI.Mult),width=0.4,size=1)+
            geom_point()+geom_hline(aes(yintercept=0))+
            geom_text(aes(y=PercUCI.Mult+abs(PercUCI.Mult/20),
                          label=star),position=position_dodge(0.8),show.legend=F,size=6)+
            facet_wrap(formula(paste0(".~",Outtitle)),scales="free",ncol=facetcol)+
            ylab("")+xlab(paste0(Predtitle))+
            scale_color_manual(values=c("black","red"),guide="none")+
            theme(axis.text.x=element_text(angle=45,hjust=1))
        }
      } else {
        if(plotOR==T){
          gg1.x<-ggplot(data=plotData,aes(x=Predictor,y=Beta,
                                          color=Interaction_Level))+
            geom_errorbar(aes(ymin=LCI,ymax=UCI),width=0.4,size=1,
                          position=position_dodge(0.75))+
            geom_point(position=position_dodge(0.75))+geom_hline(aes(yintercept=1))+
            geom_text(aes(y=UCI+abs(UCI/20),label=star),
                      position=position_dodge(0.75),show.legend=F,size=6)+
            facet_wrap(formula(paste0(".~",Outtitle)),scales="free",ncol=facetcol)+
            ylab("OR")+xlab(paste0(Predtitle))+ggtitle("Marginal Coefficients")+
            scale_color_brewer(palette="Set1",name="Interaction Level")+
            theme(axis.text.x=element_text(angle=45,hjust=1),
                  legend.position = "bottom",
                  plot.title=element_text(hjust=0.5))
          gg2<-ggplot(data=plotData.ix,
                      aes(x=Predictor,y=Beta,color=Interaction_Level))+
            geom_errorbar(aes(ymin=LCI,ymax=UCI),width=0.4,size=1,
                          position=position_dodge(0.75))+
            geom_point(position=position_dodge(0.75))+geom_hline(aes(yintercept=1))+
            geom_text(aes(y=UCI+abs(UCI/20),label=star),
                      position=position_dodge(0.75),show.legend=F,size=6)+
            facet_wrap(formula(paste0(".~",Outtitle)),scales="free",ncol=facetcol)+
            ylab("OR")+xlab(paste0(Predtitle))+ggtitle("Interaction Coefficients")+
            scale_color_brewer(palette="Dark2",name="Interaction")+
            theme(axis.text.x=element_text(angle=45,hjust=1),
                  legend.position = "bottom",
                  plot.title=element_text(hjust=0.5))
          gg1<-cowplot::plot_grid(gg1.x,gg2,nrow=2)
        } else {
          if(!plotPercChange) {
            plotData$PercBeta.Mult<-log(plotData$Beta)
            plotData$PercLCI.Mult<-log(plotData$LCI)
            plotData$PercUCI.Mult<-log(plotData$UCI)
            plotData.ix$PercBeta.Mult<-log(plotData.ix$Beta)
            plotData.ix$PercLCI.Mult<-log(plotData.ix$LCI)
            plotData.ix$PercUCI.Mult<-log(plotData.ix$UCI)
          }
          gg1.x<-ggplot(data=plotData,aes(x=Predictor,y=PercBeta.Mult,
                                          color=Interaction_Level))+
            geom_errorbar(aes(ymin=PercLCI.Mult,ymax=PercUCI.Mult),width=0.4,size=1,
                          position=position_dodge(0.75))+
            geom_point(position=position_dodge(0.75))+geom_hline(aes(yintercept=0))+
            geom_text(aes(y=PercUCI.Mult+abs(PercUCI.Mult/20),label=star),
                      position=position_dodge(0.75),show.legend=F,size=6)+
            facet_wrap(formula(paste0(".~",Outtitle)),scales="free",ncol=facetcol)+
            ylab("")+xlab(paste0(Predtitle))+ggtitle("Marginal Coefficients")+
            scale_color_brewer(palette="Set1",name="Interaction Level")+
            theme(axis.text.x=element_text(angle=45,hjust=1),
                  legend.position = "bottom",
                  plot.title = element_text(hjust=0.5))
          gg2<-ggplot(data=plotData.ix,
                      aes(x=Predictor,y=PercBeta.Mult,color=Interaction_Level))+
            geom_errorbar(aes(ymin=PercLCI.Mult,ymax=PercUCI.Mult),width=0.4,size=1,
                          position=position_dodge(0.75))+
            geom_point(position=position_dodge(0.75))+geom_hline(aes(yintercept=0))+
            geom_text(aes(y=PercUCI.Mult+abs(PercUCI.Mult/20),label=star),
                      position=position_dodge(0.75),show.legend=F,size=6)+
            facet_wrap(formula(paste0(".~",Outtitle)),scales="free",ncol=facetcol)+
            ylab("")+xlab(paste0(Predtitle))+ggtitle("Interaction Coefficients")+
            scale_color_brewer(palette="Dark2",name="Interaction")+
            theme(axis.text.x=element_text(angle=45,hjust=1),
                  legend.position = "bottom",
                  plot.title=element_text(hjust=0.5))
          gg1<-cowplot::plot_grid(gg1.x,gg2,nrow=2)
        }
      }
      
    } else if(is.null(ixterm)){
      gg1<-ggplot(data=plotData,aes(x=Predictor,y=Beta.Mult,color=Significant))+
        geom_errorbar(aes(ymin=LCI.Mult,ymax=UCI.Mult),width=0.4,size=1)+
        geom_point()+geom_hline(aes(yintercept=0))+
        geom_text(aes(y=UCI.Mult+abs(UCI.Mult/20),
                      label=star),position=position_dodge(0.75),show.legend=F,size=6)+
        facet_wrap(formula(paste0(".~",Outtitle)),scales="free",ncol=facetcol)+
        ylab("")+xlab(paste0(Predtitle))+
        scale_color_manual(values=c("black","red"),guide="none")+
        theme(axis.text.x=element_text(angle=45,hjust=1))
    } else {
      gg1.x<-ggplot(data=plotData,aes(x=Predictor,y=Beta.Mult,color=Interaction_Level))+
        geom_errorbar(aes(ymin=LCI.Mult,ymax=UCI.Mult),width=0.4,size=1,
                      position=position_dodge(0.75))+
        geom_point(position=position_dodge(0.75))+geom_hline(aes(yintercept=0))+
        geom_text(aes(y=UCI.Mult+abs(UCI.Mult/20),
                      label=star),position=position_dodge(0.75),show.legend=F,size=6)+
        facet_wrap(formula(paste0(".~",Outtitle)),scales="free",ncol=facetcol)+
        ylab("")+xlab(paste0(Predtitle))+ggtitle("Marginal Coefficients")+
        scale_color_brewer(palette="Set1",name="Interaction Level")+
        theme(axis.text.x=element_text(angle=45,hjust=1),
              legend.position = "bottom",
              plot.title=element_text(hjust=0.5))
      gg2<-ggplot(data=plotData.ix,
                  aes(x=Predictor,y=Beta.Mult,color=Interaction_Level))+
        geom_errorbar(aes(ymin=LCI.Mult,ymax=UCI.Mult),width=0.4,size=1,
                      position=position_dodge(0.75))+
        geom_point(position=position_dodge(0.75))+geom_hline(aes(yintercept=0))+
        geom_text(aes(y=UCI.Mult+abs(UCI.Mult/20),
                      label=star),position=position_dodge(0.75),show.legend=F,size=6)+
        facet_wrap(formula(paste0(".~",Outtitle)),scales="free",ncol=facetcol)+
        ylab("")+xlab(paste0(Predtitle))+ggtitle("Interaction Coefficients")+
        scale_color_brewer(palette="Dark2",name="Interaction")+
        theme(axis.text.x=element_text(angle=45,hjust=1),
              legend.position = "bottom",
              plot.title=element_text(hjust=0.5))
      gg1<-cowplot::plot_grid(gg1.x,gg2,nrow=2)
    }
  } else {
    gg1<-list()
  }
  
  retlist<-list(k1,Resmatout,gg1,lmlist)
  names(retlist)<-c("Kable","Matrix","GGplot","LMlist")
  
  return(retlist)
}
