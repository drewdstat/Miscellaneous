InteractionCoefPlot<-function(model,data,pred,ixterm,multiplier=1,coeftransform=NULL,predname=NULL,ixname=NULL,
                              outname=NULL,title=NULL,fillcolor=NULL,autotitle=F,addpvallab=F,
                              labsize=4,lengthout=50,shadebysig=F,otherix=F,robust=T,logistic=F){
  if(is.null(outname)) outname<-names(model$model)[1]
  if(is.null(predname)) predname<-pred
  if(is.null(ixname)) ixname<-ixterm
  if(is.null(fillcolor)) fillcolor<-"grey50"
  predcoefdat<-data.frame(Outcome=outname,Predictor=predname,Interaction=ixname,
                          IxLevel=seq(min(data[,ixterm],na.rm=T),max(data[,ixterm],na.rm=T),length.out=lengthout))
  ixcoeflabel<-ifelse(paste0(pred,":",ixterm)%in%names(coef(model)),paste0(pred,":",ixterm),paste0(ixterm,":",pred))
  predcoefdat$Beta<-coef(model)[pred]+predcoefdat$IxLevel*coef(model)[ixcoeflabel]
  modcoef<-coef(model)
  if(robust){
    vc<-vcovHC(model,type="HC0")
  } else {
    vc<-vcov(model)
  }
  predcoefdat$SE<-sapply(predcoefdat$IxLevel,function(x) {
    msm::deltamethod(formula(paste0("~x",which(names(modcoef)==pred),"*",multiplier,"+x",
                                    which(names(modcoef)==ixcoeflabel),"*",x,"*",multiplier)),
                     modcoef,vc)
  })
  predcoefdat$LCI<-predcoefdat$Beta-(predcoefdat$SE*1.96)
  predcoefdat$UCI<-predcoefdat$Beta+(predcoefdat$SE*1.96)
  if(!is.null(coeftransform)){predcoefdat[,c("Beta","LCI","UCI")]<-
    lapply(predcoefdat[,c("Beta","LCI","UCI")],coeftransform)}
  if(logistic){
    predcoefdat[,c("Beta","LCI","UCI")]<-lapply(predcoefdat[,c("Beta","LCI","UCI")],function(x) exp(x))
    ylabaux<-" OR"
    yint<-1
  } else {
    ylabaux<-" Coefficient"
    yint<-0
  }
  
  rugdat<-data.frame(RugX=data[,ixterm])
  if(shadebysig){
    lowshadedat<-predcoefdat[predcoefdat$UCI<0&predcoefdat$LCI<0,]
    highshadedat<-predcoefdat[predcoefdat$LCI>0&predcoefdat$UCI>0,]
    g1<-ggplot(predcoefdat,aes(x=IxLevel,y=Beta))+geom_line()+theme_bw()+
      geom_ribbon(aes(ymin=LCI,ymax=UCI),alpha=0.3,fill=fillcolor)+
      geom_hline(yintercept = yint,linetype="dashed")+
      geom_ribbon(data=lowshadedat,aes(x=IxLevel,ymin=LCI,ymax=UCI),inherit.aes=F,alpha=0.3,fill=fillcolor)+
      geom_ribbon(data=highshadedat,aes(x=IxLevel,ymin=LCI,ymax=UCI),inherit.aes=F,alpha=0.3,fill=fillcolor)+
      geom_rug(data=rugdat,aes(x=RugX),inherit.aes=F)+
      ylab(paste0(predname,ylabaux))+xlab(ixname)
  } else {
    g1<-ggplot(predcoefdat,aes(x=IxLevel,y=Beta))+geom_line()+theme_bw()+
      geom_ribbon(aes(ymin=LCI,ymax=UCI),alpha=0.3,fill=fillcolor)+
      geom_hline(yintercept = yint,linetype="dashed")+
      geom_rug(data=rugdat,aes(x=RugX),inherit.aes=F)+
      ylab(paste0(predname,ylabaux))+xlab(ixname)
  }
  
  if(is.null(title)&autotitle){
    g1<-g1+ggtitle(paste0("Change in ",outname," Per ",multiplier,
                          " Unit Increase in\n",predname," Over Levels of ",ixname))+
      theme(plot.title=element_text(hjust=0.5))
  } else if(!is.null(title)){
    g1<-g1+ggtitle(title)+theme(plot.title=element_text(hjust=0.5))
  }
  if(addpvallab){
    if(robust){
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
      robustcoef<-robustse(model)
      templab<-paste0("Interaction p=",round(robustcoef[paste0(pred,":",ixterm),4],3))
      if(robustcoef[paste0(pred,":",ixterm),4]<0.05) templab<-paste0(templab,"*")
    } else {
      templab<-paste0("Interaction p=",round(summary(model)$coef[paste0(pred,":",ixterm),4],3))
      if(summary(model)$coef[paste0(pred,":",ixterm),4]<0.05) templab<-paste0(templab,"*")
    }
    
    g1<-g1+annotate(geom="text",x=mean(range(data[,ixterm])),y=max(predcoefdat$UCI),
                    label=templab,size=labsize,hjust=0.5,vjust=1)
    # +geom_text(data=labdat,aes(x=X,y=X,label=Label),inherit.aes=F,size=labsize,hjust=0.2,vjust=0.5)
  }
  retlist<-list(GGplot=g1,Matrix=predcoefdat)
  
  if(otherix){
    othercoefdat<-data.frame(Outcome=outname,Predictor=ixname,Interaction=predname,
                            IxLevel=seq(min(data[,pred],na.rm=T),max(data[,pred],na.rm=T),length.out=lengthout))
    othercoefdat$Beta<-coef(model)[ixterm]+othercoefdat$IxLevel*coef(model)[paste0(pred,":",ixterm)]
    modcoef<-coef(model)
    if(robust){
      vc<-vcovHC(model,type="HC0")
    } else {
      vc<-vcov(model)
    }
    othercoefdat$SE<-sapply(othercoefdat$IxLevel,function(x) {
      msm::deltamethod(formula(paste0("~x",which(names(modcoef)==ixterm),"*",multiplier,"+x",
                                      which(names(modcoef)==paste0(pred,":",ixterm)),"*",x,"*",multiplier)),
                       modcoef,vc)
    })
    othercoefdat$LCI<-othercoefdat$Beta-(othercoefdat$SE*1.96)
    othercoefdat$UCI<-othercoefdat$Beta+(othercoefdat$SE*1.96)
    
    if(!is.null(coeftransform)){othercoefdat[,c("Beta","LCI","UCI")]<-
      lapply(othercoefdat[,c("Beta","LCI","UCI")],coeftransform)}
    if(logistic){
      othercoefdat[,c("Beta","LCI","UCI")]<-lapply(othercoefdat[,c("Beta","LCI","UCI")],function(x) exp(x))
    } 
    if(logistic){othercoefdat[,c("Beta","LCI","UCI")]<-
      lapply(othercoefdat[,c("Beta","LCI","UCI")],function(x) exp(x))}
    otherrugdat<-data.frame(RugX=data[,pred])
    if(shadebysig){
      lowshadedat<-othercoefdat[othercoefdat$UCI<0&othercoefdat$LCI<0,]
      highshadedat<-othercoefdat[othercoefdat$LCI>0&othercoefdat$UCI>0,]
      g2<-ggplot(othercoefdat,aes(x=IxLevel,y=Beta))+geom_line()+theme_bw()+
        geom_ribbon(aes(ymin=LCI,ymax=UCI),alpha=0.3,fill=fillcolor)+
        geom_hline(yintercept = yint,linetype="dashed")+
        geom_ribbon(data=lowshadedat,aes(x=IxLevel,ymin=LCI,ymax=UCI),inherit.aes=F,alpha=0.3,fill=fillcolor)+
        geom_ribbon(data=highshadedat,aes(x=IxLevel,ymin=LCI,ymax=UCI),inherit.aes=F,alpha=0.3,fill=fillcolor)+
        geom_rug(data=otherrugdat,aes(x=RugX),inherit.aes=F)+
        ylab(paste0(ixname,ylabaux))+xlab(predname)
    } else {
      g2<-ggplot(othercoefdat,aes(x=IxLevel,y=Beta))+geom_line()+theme_bw()+
        geom_ribbon(aes(ymin=LCI,ymax=UCI),alpha=0.3,fill=fillcolor)+
        geom_hline(yintercept = yint,linetype="dashed")+
        geom_rug(data=otherrugdat,aes(x=RugX),inherit.aes=F)+
        ylab(paste0(ixname,ylabaux))+xlab(predname)
    }
    
    if(is.null(title)&autotitle){
      g2<-g2+ggtitle(paste0("Change in ",outname," Per ",multiplier,
                            " Unit Increase in\n",predname," Over Levels of ",ixname))+
        theme(plot.title=element_text(hjust=0.5))
    } else if(!is.null(title)){
      g2<-g2+ggtitle(title)+theme(plot.title=element_text(hjust=0.5))
    }
    if(addpvallab){
      if(robust){
        templab<-paste0("Interaction p=",round(robustcoef[paste0(pred,":",ixterm),4],3))
        if(robustcoef[paste0(pred,":",ixterm),4]<0.05) templab<-paste0(templab,"*")
      } else {
        templab<-paste0("Interaction p=",round(summary(model)$coef[paste0(pred,":",ixterm),4],3))
        if(summary(model)$coef[paste0(pred,":",ixterm),4]<0.05) templab<-paste0(templab,"*")
      }
      
      g2<-g2+annotate(geom="text",x=mean(range(data[,pred])),y=max(othercoefdat$UCI),
                      label=templab,size=labsize,hjust=0.5,vjust=1)
      #+geom_text(data=labdat,aes(x=X,y=X,label=Label),inherit.aes=F,size=labsize,vjust=0.2,hjust=0.5)
    }
    retlist<-list(MainGGplot=g1,OtherGGplot=g2,MainMatrix=predcoefdat,OtherMatrix=othercoefdat)
  }
  return(retlist)
}