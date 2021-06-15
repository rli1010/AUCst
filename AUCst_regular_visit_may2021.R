library(reshape2)

fp.func3=function(sk,t,theta){
  AUC.st.mat=c()
  for (s.curr in sk)
  { 
    tsvec=cbind(1, s.curr, t, s.curr^2, t^2, s.curr*t, t^3, s.curr^3, (t^2)*s.curr, t*(s.curr^2)) #specify basis
    teta=as.vector(as.matrix(tsvec)%*%theta)
    teta[teta>=100]=100
    teta[teta<=-100]=-100
    AUC.st=exp(teta)/(1+exp(teta))
    AUC.st.mat=rbind(AUC.st.mat, AUC.st)}
  return(AUC.st.mat)}   


logL=function(theta,sk,d.times,n1.h.mat,n2.h.mat,eta.sk.mat,sk.le.dh,Ewt.die){
  AUC.seq.mat=fp.func3(sk, d.times,theta)    
  logL.seq= Ewt.die*eta.sk.mat*sk.le.dh*(n1.h.mat*log(AUC.seq.mat)+n2.h.mat*log(1-AUC.seq.mat))
  logL=-sum(logL.seq)}


main1.sub.func=function(da, da.long, sk, par0=c(0.5,0,0,0,0,0,0,0,0,0), tseq.eval, resample=0, nsap=1,seed=12345){
  # Compute the AUC(s,t)  and related SE
  # da: the survival dataset
  # da.long: the long-format dataset with a longitudinal biomarker
  # sk: biomarker measurement time under the regular visit scenario
  # par0: initial value for theta, used in optimization
  # tseq.eval: time points for calculating the AUC
  # resample: indictacor of resampling. (resample=1: conduct perturbed-resampling, resample=0: do not contudct the perturbed-resampling)  
  # nsap: number of perturbation. To obtain SE, set resample to 1 and nsap to a large number (e.g., 250)
  
  s.num=length(sk)
  t.num=length(tseq.eval)
  resap.AUC.mat=c()
  n=dim(da)[1]  #sample size (do not allow missing data)
  
  long.basic=data.frame(cbind(rep(da$obs_id,rep(s.num,n)),rep(sk,n)))
  names(long.basic)=c("obs_id","vtime")
  marker0=merge(long.basic,da.long,by=c("obs_id","vtime"),all=T)
  marker0$eta=!is.na(marker0$Zt)*1
  marker0=merge(marker0,da,by="obs_id")

  set.seed(seed) #seed is used for resampling
  for( sap in 1:nsap){ 
    
    if (resample==1){ EXPwt=rexp(n, 1)}  #perturbation weight
    if (resample==0){ EXPwt=rep(1, n) }
    
    d.times=sort(unique(da$Y[da$delta==1]))
    md=length(d.times)
    Ymat=matrix(rep(da$Y,md),n)
    dt.mat=matrix(rep(d.times,rep(n,md)),n)
    dt.mat.tb=dt.mat+0.00001; 
    n.h=apply((Ymat>dt.mat.tb),2,sum)
    death.mat=(Ymat==dt.mat)*da$delta #currently do not allow tied failure times;
    marker0$eta[marker0$eta==0]=0
    RS.data=marker0
    
    # Calculate the unweighted n1 and n2, and perturbed n1 and n2
    RS.seq.list=RS.mat.list=RS.d.list=RS.dmat.list=n1.h.list=n2.h.list=list()
    Wn1.h.list=Wn2.h.list=list() 
    n1n2.wt=matrix(rep(EXPwt,md),n)
    for (i in 1:s.num)
    {
      RS.seq.list[[i]]=RS.data[RS.data$vtime==sk[i],]
      RS.mat.list[[i]]<-matrix(rep(RS.seq.list[[i]]$Zt,md),n)
      RS.d.list[[i]]=apply(RS.mat.list[[i]]*death.mat,2,sum,na.rm=TRUE)
      RS.d.list[[i]][RS.d.list[[i]]==0]<-NA 
      RS.dmat.list[[i]]=matrix(rep(RS.d.list[[i]],rep(n,md)),n)
      n1.h.list[[i]]=apply((RS.dmat.list[[i]]>RS.mat.list[[i]])*(Ymat>dt.mat.tb),2,sum,na.rm=TRUE)
      n2.h.list[[i]]=apply((RS.dmat.list[[i]]<=RS.mat.list[[i]])*(Ymat>dt.mat.tb),2,sum,na.rm=TRUE)
      
      Wn1.h.list[[i]]=apply((RS.dmat.list[[i]]>RS.mat.list[[i]])*(Ymat>dt.mat.tb)*n1n2.wt,2,sum,na.rm=TRUE)
      Wn2.h.list[[i]]=apply((RS.dmat.list[[i]]<=RS.mat.list[[i]])*(Ymat>dt.mat.tb)*n1n2.wt,2,sum,na.rm=TRUE)
    }
    
    n1.h.mat=n2.h.mat=c()
    Wn1.h.mat=Wn2.h.mat=c()
    eta.sk.mat=c()
    for (i in 1:s.num)
    {
      n1.h.mat=rbind(n1.h.mat,n1.h.list[[i]])
      n2.h.mat=rbind(n2.h.mat,n2.h.list[[i]])
      Wn1.h.mat=rbind(Wn1.h.mat, Wn1.h.list[[i]])
      Wn2.h.mat=rbind(Wn2.h.mat, Wn2.h.list[[i]])
      eta.sk.mat=rbind(eta.sk.mat, RS.d.list[[i]])
    }  
    eta.sk.mat=(eta.sk.mat/eta.sk.mat)
    eta.sk.mat[is.na(eta.sk.mat)]=0
    sk.mat=matrix(rep(sk,each=md), ncol=md, byrow=TRUE)
    d.times.mat=matrix(rep(d.times,rep(s.num,md)),s.num)
    sk.le.dh=(sk.mat<=d.times.mat)
    
    # Extract the weights(Ewt.die) for subjects who experienced the event of interest
    da0<-cbind(da$Y,da$delta,EXPwt) 
    da0.sort=da0[order(da0[,1], decreasing = F),]
    Ewt.dseq=da0.sort[da0.sort[,2]==1,3] 
    Ewt.die=matrix(rep(Ewt.dseq,rep(s.num,md)),s.num) 
    
    #Computed the perturbed estimates of theta and AUC(s,t)
    resap.fit=optim(par=par0, fn=logL, d.times=d.times, n1.h.mat=Wn1.h.mat, n2.h.mat=Wn2.h.mat, sk=sk, 
                     eta.sk.mat=eta.sk.mat, sk.le.dh=sk.le.dh, control = list(maxit = 10000,reltol=1e-8),Ewt.die=Ewt.die)
    
    s.mat=matrix(rep(sk,each=t.num), ncol=t.num, byrow=TRUE)
    t.mat=matrix(rep(tseq.eval,rep(s.num,t.num)),s.num)
    slet.mat=1*(s.mat<=t.mat)
    slet.mat[slet.mat==0]=NA
    resap.AUC0=fp.func3(sk,tseq.eval,resap.fit$par)*slet.mat 
    resap.AUC.vec= c(as.vector(t(resap.AUC0)), resap.fit$convergence)
    resap.AUC.mat=rbind(resap.AUC.mat, resap.AUC.vec)
    print(sap)
  }
  
  #Remove non-converged or abnormal perturbed AUC(s,t) estimates 
  #Compute the standard deviation of these resampled AUCs
  resap.AUC.mat=as.data.frame(resap.AUC.mat)
  names(resap.AUC.mat)[ncol(resap.AUC.mat)]<-"cov1"
  resap.AUC.mat$minAUC=apply(resap.AUC.mat[,1:(t.num*s.num)], 1, FUN=min, na.rm=T)
  resap.AUC.mat$maxAUC=apply(resap.AUC.mat[,1:(t.num*s.num)], 1, FUN=max, na.rm=T)
  resap.AUC.mat$cov2=ifelse(resap.AUC.mat$minAUC<0.00001,1,0)
  resap.AUC.mat$cov3=ifelse(resap.AUC.mat$maxAUC>0.99999,1,0)
  resap.AUC.mat.cov=as.matrix(resap.AUC.mat[resap.AUC.mat$cov1==0&resap.AUC.mat$cov2==0&resap.AUC.mat$cov3==0,])  
  RESAP.sd=apply(resap.AUC.mat.cov[,1:(t.num*s.num)], 2, FUN=sd) 
  RESAP.sd.mat=t(matrix(RESAP.sd, t.num, s.num))

  #Compute the estimated theta and AUCs without perturbed-resampling
  fit.eta=optim(par=par0, fn=logL, d.times=d.times, n1.h.mat=n1.h.mat, n2.h.mat=n2.h.mat, sk=sk, 
                eta.sk.mat=eta.sk.mat, sk.le.dh=sk.le.dh, control = list(maxit = 10000,reltol=1e-8),Ewt.die=1)
  
  theta.est=fit.eta$par
  conv.est=fit.eta$convergence
  AUC.est=fp.func3(sk,tseq.eval,theta.est)*slet.mat
  return(list(conv=conv.est, AUC=AUC.est, theta=theta.est, RS= RS.seq.list, SE=RESAP.sd.mat))
}

####################################
#########Analysis Start Here########
####################################
# da.short includes subject id obs_id,  observed event time Y, and status delta  
# da.long  includes subject id obs_id,  measurement time vtime, and biomaker value Zt

da.sim.short <- read.delim(url("https://raw.githubusercontent.com/rli1010/AUCst/main/Reg_data_da.sim_short.txt"),sep=" ")
da.sim.long <- read.delim(url("https://raw.githubusercontent.com/rli1010/AUCst/main/Reg_data_da.sim_long.txt"),sep=" ")

sort(unique(da.sim.long$vtime))  #regular visit times in the data 


# Estimate the AUC(s,t); to get SE set resample=1 and nsap to a large number (e.g., 250)
obj.i=main1.sub.func(da=da.sim.short, 
                     da.long=da.sim.long,
                     sk=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2), 
                     tseq.eval=c(0.4,1.11,1.50,1.90), resample =0, nsap=3)

# AUC.result: the estimated AUC(s,t), with s corresponding to row and t corresponding to column
# Only has value for s<=t

AUC.result=obj.i$AUC          
AUC.result

