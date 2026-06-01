################################################################################
#
# This script computes the numerical experiments presented in:
# Singh, S. K., Cruz, F. R. B. & Quinino, R. C. (2026) Bayesian Inference
# in Erlang Single Server Queueing Model Based on Queue Length, American 
# Journal of Management Sciences (in press).
#
# Programmed by:
#
# Frederico R. B. Cruz, Roberto C. Quinino
# Universidade Federal de Minas Gerais
# E-mail: {fcruz,roberto}@est.ufmg.br
# (c) 2026 Cruz & Quinino
# v.2026.05.26
#
################################################################################

rm(list=ls())
library(this.path) # use relative path
setwd(this.path::here())

################################################################################
# Export to Excel
################################################################################
ExExcel <- function(data){
  write.table(format(data, decimal.mark = '.'),
              'clipboard', sep='\t', na = '*')}

################################################################################
# remove all variables, except functions
################################################################################
rem_var <- function() {
  all_objects <- ls(envir = .GlobalEnv)
  functions <- all_objects[sapply(all_objects,
            function(x) is.function(get(x, envir = .GlobalEnv)))]
  rm(list = setdiff(all_objects, functions), envir = .GlobalEnv)
}

################################################################################
# discrete event simulation for a M/Er/1 queue
################################################################################
rMEr1<-function(n,rho,r) { # simulate
  ##############################################################################
  # initial set up
  ##############################################################################
  users=0
  ssize=0
  nusers=numeric(n)
  x=numeric(n)
  fx=numeric(n)
  lambda=1.0
  warmup<-100/lambda
  mu=lambda/rho
  event.type<-character(1000)
  event.time<-numeric(1000)
  #  cat("rMEr1: processing first events\n")
  # generate first arrival dt~exp(lambda)
  time.event<-0.0+rexp(n=1,rate=lambda)
  # enlist first event
  event.size<-1
  event.time[event.size]<-time.event
  event.type[event.size]<-"arr"
  # generate final event (end)
  time.event<-Inf
  # enlist final event
  event.size<-2
  event.time[event.size]<-time.event
  event.type[event.size]<-"end"
  #  cat("event list\n")
  #  for (i in 1:event.size) cat(event.type[i]," ",event.time[i],"\n")
  while (event.type[1]!="end") {
    ################################################################################
    if (event.type[1]=="arr") {
      #      cat("rMEr1: processing arrival\n")
      ################################################################################
      users=users+1
      if (users==1) { # server is free
        # generate departure dt~rgamma(n,shape,rate)
        time.event<-event.time[1]+rgamma(n=1,shape=r,rate=r*mu)
        # find position in list and enlist event
        event.size<-event.size+1
        i<-event.size
        while ((i>1)&&(event.time[i-1]>time.event)) {
          event.type[i]<-event.type[i-1]
          event.time[i]<-event.time[i-1]
          i<-i-1}
        event.time[i]<-time.event
        event.type[i]<-"dep"}
      # generate next arrival dt~exp(lambda)
      time.event<-event.time[1]+rexp(n=1,rate=lambda)
      # find position in list and enlist event
      event.size<-event.size+1
      i<-event.size
      while ((i>1)&&(event.time[i-1]>time.event)) {
        event.type[i]<-event.type[i-1]
        event.time[i]<-event.time[i-1]
        i<-i-1}
      event.time[i]<-time.event
      event.type[i]<-"arr"
      # remove current event (arrival)
      event.size<-event.size-1
      for (i in 1:event.size) {
        event.type[i]<-event.type[i+1]
        event.time[i]<-event.time[i+1]}
      ################################################################################
      
    } else if (event.type[1]=="dep") {
      #      cat("rMEr1: processing departure\n")
      ################################################################################
      # register number of users in system at departure
      users=users-1
      if (event.time[1]>warmup) {
        ssize=ssize+1
        nusers[ssize]=users
        x[ssize]=event.time[1]
        fx[ssize]=users}
      if (users!=0) { # queue is not empty
        # generate departure dt~rgamma(n,shape,rate)
        time.event<-event.time[1]+rgamma(n=1,shape=r,rate=r*mu)
        # find position in list and enlist event
        event.size<-event.size+1
        i<-event.size
        while ((i>1)&&(event.time[i-1]>time.event)) {
          event.type[i]<-event.type[i-1]
          event.time[i]<-event.time[i-1]
          i<-i-1}
        event.time[i]<-time.event
        event.type[i]<-"dep"}
      # remove current event (departure)
      event.size<-event.size-1
      for (i in 1:event.size) {
        event.type[i]<-event.type[i+1]
        event.time[i]<-event.time[i+1]}
      if (ssize>=n) {
        # end of simulation
        event.type[1]<-"end"}
      ################################################################################
    } else {
      cat("rMEr1: error, unknown event\n")
      ################################################################################
      return(0)}
    #    cat("event list\n")
    #    for (i in 1:event.size) {
    #      cat(event.type[i]," ",event.time[i],"\n")}
  }
  return(nusers)}

################################################################################
# expected number in the system, L, as function of rho
################################################################################
LMEr1<-function(rho,r) {
  L<-LqMEr1(rho,r)+rho
  return(L)}

################################################################################
# average queue length, Lq, as function of rho
################################################################################
LqMEr1<-function(rho,r) {
  Lq<-((r+1)/(2*r))*(rho^2)/(1-rho)
  return(Lq)}

################################################################################
# traffic intensity, rho, as function of L (fLx=LMEr1) or Lq (fLx=LqMEr1)
################################################################################
rhoLxMEr1<-function(Lx,r,fLx){
  epsilon<-1E-06
  gr<-(sqrt(5.0)+1)/2
  rhoL<-0
  rhoU<-1
  frhoL<-abs(Lx-fLx(rhoL,r))
  frhoU<-abs(Lx-fLx(rhoU,r))
  rho1<-rhoU-(rhoU-rhoL)/gr
  rho2<-rhoL+(rhoU-rhoL)/gr
  frho1<-abs(Lx-fLx(rho1,r))
  frho2<-abs(Lx-fLx(rho2,r))
  #  cat("Search=[", rhoL, ",", rhoU, "]->")
  #  cat("[", frhoL, ",", frhoU, "]\n")
  #  cat("[", rho1, ",", rho2, "]=")
  #  cat("[", frho1, ",", frho2, "]\n")
  repeat {
    if (frho1<frho2) {
      rhoU<-rho2
      frhoU<-frho2
      rho2<-rho1
      frho2<-frho1
      rho1<-rhoU-(rhoU-rhoL)/gr
      frho1<-abs(Lx-fLx(rho1,r))
      #       cat("new rhoU=", rhoU, "\n")
    } else {
      rhoL<-rho1
      frhoL<-frho1
      rho1<-rho2
      frho1<-frho2
      rho2<-rhoL+(rhoU-rhoL)/gr
      frho2<-abs(Lx-fLx(rho2,r))
      #       cat("new rhoL=", rhoL, "\n")
    }
    #   cat("[", rho1, ",", rho2, "]=[", frho1, ",", frho2, "]\n")
    if (frhoU<epsilon) break
  }
  return(rhoU)}

################################################################################
# traffic intensity, rho, as function of Lq
################################################################################
rhoLqMEr1<-function(Lq,r){
  return(rhoLxMEr1(Lq,r,LqMEr1))}

################################################################################
# traffic intensity, rho, as function of L
################################################################################
rhoLMEr1<-function(L,r){
  return(rhoLxMEr1(L,r,LMEr1))}

################################################################################
# maximum likelihood estimate for rho (state-of-the-art)
################################################################################
MEr1RhoMLESoA<-function(samp,r) {
  Eps.MLE<-1E-06
  n<-length(samp)
  y<-sum(samp)
  if (y<n)
    return(y/n)
  else
    return(1-Eps.MLE)}

################################################################################
# maximum likelihood estimate for rho
################################################################################
MEr1RhoMLEf<-function(samp,r) {
  Eps.MLE<-1E-06
#  log-likelihood function
  loglike.f<-function(rho,samp,r) {
    n<-length(samp)
    infty<-max(samp)
    nij<-matrix(0,infty+1,infty+1)
    phi1<-0.0
    phi2<-0.0
    loglikef=0.0
    i<-2
    while (i<=n) {
     # cat("i=",i,"\n")
      nij[samp[i-1]+1,samp[i]+1]<-nij[samp[i-1]+1,samp[i]+1]+1
      i<-i+1}
    j=0
    while (j<=infty) {
     # cat("j=",j,"\n")
      phi1<-phi1+j*nij[0+1,j+1]
      phi2<-phi2+(r+j)*nij[0+1,j+1]
      loglikef=loglikef+nij[0+1,j+1]*log(choose(r+j-1,j))
      j=j+1}
    i=1
    while (i<=infty) {
      j<-i-1
      while (j<=infty) {
#        cat("(i,j)=(",i,",",j,")\n")
        phi1<-phi1+(j-i+1)*nij[i+1,j+1]
        phi2<-phi2+(r+j-i+1)*nij[i+1,j+1]
        loglikef=loglikef+nij[i+1,j+1]*log(choose(r+j-i,j-i+1))
        j<-j+1}
      i<-i+1}
    loglikef=loglikef+phi1*log(rho/r)-phi2*log(1+rho/r)
   # cat("phi1=",phi1,"\t","phi2=",phi2,"\t",loglikef,"\n")
    return(loglikef)}
#  likelihood function
  like.f<-function(rho,samp,r) {
    n<-length(samp)
    infty<-max(samp)
    nij<-matrix(0,infty+1,infty+1)
    phi1<-0.0
    phi2<-0.0
    likef=1.0
    i<-2
    while (i<=n) {
      #      cat("i=",i,"\n")
      nij[samp[i-1]+1,samp[i]+1]<-nij[samp[i-1]+1,samp[i]+1]+1
      i<-i+1}
    j=0
    while (j<=infty) {
      #      cat("j=",j,"\n")
      phi1<-phi1+j*nij[0+1,j+1]
      phi2<-phi2+(r+j)*nij[0+1,j+1]
      # likef=likef*(choose(r+j-1,j))^nij[0+1,j+1]
      j=j+1}
    i=1
    while (i<=infty) {
      j<-i-1
      while (j<=infty) {
        #        cat("(i,j)=(",i,",",j,")\n")
        phi1<-phi1+(j-i+1)*nij[i+1,j+1]
        phi2<-phi2+(r+j-i+1)*nij[i+1,j+1]
        # likef=likef*(choose(r+j-i,j-i+1))^nij[i+1,j+1]
        j<-j+1}
      i<-i+1}
    likef=likef*(rho/r)^phi1*(1+rho/r)^(-phi2)
   # cat("phi1=",phi1,"\t","phi2=",phi2,"\t",likef,"\n")
    return(likef)}
#   numerical maximization
    # rho<-seq(0,1,0.01)
    # plot(rho,like.f(rho,samp,r),type="l")
    # res<-optimize(loglike.f,c(Eps.MLE,1-Eps.MLE),samp,r,maximum=TRUE,tol=Eps.MLE)$maximum
#   analytical maximization
#   find phi1 and phi2
  n<-length(samp)
  infty<-max(samp)
  nij<-matrix(0,infty+1,infty+1)
  phi1<-0.0
  phi2<-0.0
  i<-2
  while (i<=n) {
    #    cat("i=",i,"\n")
    nij[samp[i-1]+1,samp[i]+1]<-nij[samp[i-1]+1,samp[i]+1]+1
    i<-i+1}
  j=0
  while (j<=infty) {
    #    cat("j=",j,"\n")
    phi1<-phi1+j*nij[0+1,j+1]
    phi2<-phi2+(r+j)*nij[0+1,j+1]
    j=j+1}
  i=1
  while (i<=infty) {
    j<-i-1
    while (j<=infty) {
      #      cat("(i,j)=(",i,",",j,")\n")
      phi1<-phi1+(j-i+1)*nij[i+1,j+1]
      phi2<-phi2+(r+j-i+1)*nij[i+1,j+1]
      j<-j+1}
    i<-i+1}
  #  cat("phi1=",phi1,"\t","phi2=",phi2,"\n")
  #   find maximum
    if (r*phi1<(phi2-phi1))
      return(r*phi1/(phi2-phi1))
    else
      return(1-Eps.MLE)}

################################################################################
# maximum likelihood estimate for Lq
################################################################################
MEr1LqMLEf<-function(samp,r,rho0=0.95) {
  rhoMLE<-min(MEr1RhoMLEf(samp,r),rho0)
  return(LqMEr1(rhoMLE,r))}

################################################################################
# Gauss hypergeometric prior
################################################################################
dGaussHypf<-function(x,a,b,c,z) {
  #  cat("dGaussHypf(",x,",",a,",",b,",",c,",",z,"):\n")
  GaussHypItg<-function(u,a,b,c,z) {
    U<-u^(b-1)*(1-u)^(c-b-1)*(1-u*z)^(-a)
    #    cat("x^(b-1):\n",x^(b-1),"\n")
    #    cat("(1-x)^(c-b-1):\n",(1-x)^(c-b-1),"\n")
    #    cat("(1-x*z)^(-a):\n",(1-x*z)^(-a),"\n")
    #    for (i in 1:length(x)) {
    #      cat("U:\n",x[i]^(b-1),"*",(1-x[i])^(c-b-1),"*",(1-x[i]*z)^(-a),"=",U[i],"\n")
    #      cat("logU:\n",log(x[i]^(b-1)),"+",log((1-x[i])^(c-b-1)),"+",log(1-x[i]*z)^(-a),"=",exp((b-1)*log(x[i])+(c-b-1)*log(1-x[i])-a*log(1-x[i]*z)),"\n")}
    #    U<-exp((b-1)*log(x)+(c-b-1)*log(1-x)-a*log(1-x*z))
    return(U)}
  d<-GaussHypItg(x,a,b,c,z)
  GaussHyp1F2<-integrate(GaussHypItg,0,1,a,b,c,z)[[1]]/beta(b,c-b)
  dGaussHyp<-d/beta(b,c-b)/GaussHyp1F2
  #  cat("dGaussHypf(",x,")=",d,"/",beta(b,c-b),"/",GaussHyp1F2,"=",dGaussHyp,"\n")
  return(dGaussHyp)}

################################################################################
# Gauss hypergeometric posterior
################################################################################
dGaussHypPostf<-function(x,samp,r,a,b,c,z) {
  GaussHypItg<-function(u,a,b,c,z) {
    U<-u^(b-1)*(1-u)^(c-b-1)*(1-u*z)^(-a)
    #    cat("x^(b-1):\n",x^(b-1),"\n")
    #    cat("(1-x)^(c-b-1):\n",(1-x)^(c-b-1),"\n")
    #    cat("(1-x*z)^(-a):\n",(1-x*z)^(-a),"\n")
    #    for (i in 1:length(x)) {
    #      cat("U:\n",x[i]^(b-1),"*",(1-x[i])^(c-b-1),"*",(1-x[i]*z)^(-a),"=",U[i],"\n")
    #      cat("logU:\n",log(x[i]^(b-1)),"+",log((1-x[i])^(c-b-1)),"+",log(1-x[i]*z)^(-a),"=",exp((b-1)*log(x[i])+(c-b-1)*log(1-x[i])-a*log(1-x[i]*z)),"\n")}
    #    U<-exp((b-1)*log(x)+(c-b-1)*log(1-x)-a*log(1-x*z))
    return(U)}
  n<-length(samp)
  infty<-max(samp)
  nij<-matrix(0,infty+1,infty+1)
  phi1<-0.0
  phi2<-0.0
  i<-2
  while (i<=n) {
    #    cat("i=",i,"\n")
    nij[samp[i-1]+1,samp[i]+1]<-nij[samp[i-1]+1,samp[i]+1]+1
    i<-i+1}
  j=0
  while (j<=infty) {
    #    cat("j=",j,"\n")
    phi1<-phi1+j*nij[0+1,j+1]
    phi2<-phi2+(r+j)*nij[0+1,j+1]
    j=j+1}
  i=1
  while (i<=infty) {
    j<-i-1
    while (j<=infty) {
      #      cat("(i,j)=(",i,",",j,")\n")
      phi1<-phi1+(j-i+1)*nij[i+1,j+1]
      phi2<-phi2+(r+j-i+1)*nij[i+1,j+1]
      j<-j+1}
    i<-i+1}
  # cat("phi1=",phi1,"\t","phi2=",phi2,"\n")
  d<-GaussHypItg(x,phi2+a,phi1+b,phi1+c,z)
  GaussHyp1F2<-integrate(GaussHypItg,0,1,phi2+a,phi1+b,phi1+c,z)[[1]]/beta(phi1+b,c-b)
  dGaussHypPostf<-d/beta(phi1+b,c-b)/GaussHyp1F2
  #  cat("dGaussHypPostf(",x,")=",d,"/",beta(phi1+b,c-b),"/",GaussHyp1F2,"=",dGaussHyp,"\n")
  return(dGaussHypPostf)}

################################################################################
# MEr1RhoSELFf implements Bayesian estimation under SELF for rho
################################################################################
MEr1RhoSELFf<-function(samp,r,a,b,c,z) {
  GaussHypItg<-function(u,a,b,c,z) {
    U<-u^(b-1)*(1-u)^(c-b-1)*(1-u*z)^(-a)
    return(U)}
  n<-length(samp)
  infty<-max(samp)
  nij<-matrix(0,infty+1,infty+1)
  phi1<-0.0
  phi2<-0.0
  i<-2
  while (i<=n) {
    #    cat("i=",i,"\n")
    nij[samp[i-1]+1,samp[i]+1]<-nij[samp[i-1]+1,samp[i]+1]+1
    i<-i+1}
  j=0
  while (j<=infty) {
    #    cat("j=",j,"\n")
    phi1<-phi1+j*nij[0+1,j+1]
    phi2<-phi2+(r+j)*nij[0+1,j+1]
    j=j+1}
  i=1
  while (i<=infty) {
    j<-i-1
    while (j<=infty) {
      #      cat("(i,j)=(",i,",",j,")\n")
      phi1<-phi1+(j-i+1)*nij[i+1,j+1]
      phi2<-phi2+(r+j-i+1)*nij[i+1,j+1]
      j<-j+1}
    i<-i+1}
  #  cat("phi1=",phi1,"\t","phi2=",phi2,"\n")
  GaussHyp1F2a<-integrate(GaussHypItg,0,1,phi2+a,phi1+b+1,phi1+c+1,z)[[1]]/beta(phi1+b+1,c-b)
  GaussHyp1F2b<-integrate(GaussHypItg,0,1,phi2+a,phi1+b  ,phi1+c  ,z)[[1]]/beta(phi1+b  ,c-b)
  MEr1RhoSELF<-(phi1+b)*GaussHyp1F2a/(phi1+c)/GaussHyp1F2b
  #  cat("MEr1RhoSELF=",phi1+b,"/",GaussHyp1F2a,"/",(phi1+c),"/",GaussHyp1F2b,"=",MEr1RhoSELFf,"\n")
  return(MEr1RhoSELF)}

################################################################################
# MEr1LqSELFf implements Bayesian estimation under SELF for Lq
################################################################################
MEr1LqSELFf<-function(samp,r,a,b,c,z) {
  GaussHypItg<-function(u,a,b,c,z) {
    U<-u^(b-1)*(1-u)^(c-b-1)*(1-u*z)^(-a)
    return(U)}
  n<-length(samp)
  infty<-max(samp)
  nij<-matrix(0,infty+1,infty+1)
  phi1<-0.0
  phi2<-0.0
  i<-2
  while (i<=n) {
    #    cat("i=",i,"\n")
    nij[samp[i-1]+1,samp[i]+1]<-nij[samp[i-1]+1,samp[i]+1]+1
    i<-i+1}
  j=0
  while (j<=infty) {
    #    cat("j=",j,"\n")
    phi1<-phi1+j*nij[0+1,j+1]
    phi2<-phi2+(r+j)*nij[0+1,j+1]
    j=j+1}
  i=1
  while (i<=infty) {
    j<-i-1
    while (j<=infty) {
      #      cat("(i,j)=(",i,",",j,")\n")
      phi1<-phi1+(j-i+1)*nij[i+1,j+1]
      phi2<-phi2+(r+j-i+1)*nij[i+1,j+1]
      j<-j+1}
    i<-i+1}
  #  cat("phi1=",phi1,"\t","phi2=",phi2,"\n")
  GaussHyp1F2a<-integrate(GaussHypItg,0,1,phi2+a,phi1+b+2,phi1+c+1,z)[[1]]/beta(phi1+b+2,c-b-1)
  GaussHyp1F2b<-integrate(GaussHypItg,0,1,phi2+a,phi1+b  ,phi1+c  ,z)[[1]]/beta(phi1+b  ,c-b  )
  MEr1LqSELF<-(r+1)*(phi1+b+1)*(phi1+b)*GaussHyp1F2a/(2*r)/(phi1+c)/(c-b-1)/GaussHyp1F2b
  # cat("MEr1LqSELF=", MEr1LqSELF,"\n")
  return(MEr1LqSELF)}

################################################################################
# MEr1RhoPLFf implements Bayesian estimation under PLF for rho in M/Er/1 queues
################################################################################
MEr1RhoPLFf<-function(samp,r,a,b,c,z) {
  GaussHypItg<-function(u,a,b,c,z) {
    U<-u^(b-1)*(1-u)^(c-b-1)*(1-u*z)^(-a)
    return(U)}
  n<-length(samp)
  infty<-max(samp)
  nij<-matrix(0,infty+1,infty+1)
  phi1<-0.0
  phi2<-0.0
  i<-2
  while (i<=n) {
    #    cat("i=",i,"\n")
    nij[samp[i-1]+1,samp[i]+1]<-nij[samp[i-1]+1,samp[i]+1]+1
    i<-i+1}
  j=0
  while (j<=infty) {
    #    cat("j=",j,"\n")
    phi1<-phi1+j*nij[0+1,j+1]
    phi2<-phi2+(r+j)*nij[0+1,j+1]
    j=j+1}
  i=1
  while (i<=infty) {
    j<-i-1
    while (j<=infty) {
      #      cat("(i,j)=(",i,",",j,")\n")
      phi1<-phi1+(j-i+1)*nij[i+1,j+1]
      phi2<-phi2+(r+j-i+1)*nij[i+1,j+1]
      j<-j+1}
    i<-i+1}
  #  cat("phi1=",phi1,"\t","phi2=",phi2,"\n")
  GaussHyp1F2a<-integrate(GaussHypItg,0,1,phi2+a,phi1+b+2,phi1+c+2,z)[[1]]/beta(phi1+b+2,c-b)
  GaussHyp1F2b<-integrate(GaussHypItg,0,1,phi2+a,phi1+b  ,phi1+c  ,z)[[1]]/beta(phi1+b  ,c-b)
  MEr1RhoPLF<-sqrt((phi1+b+1)*(phi1+b)*GaussHyp1F2a/(phi1+c+1)/(phi1+c)/GaussHyp1F2b)
#  cat("MEr1RhoPLF=",phi1+b+1,"/",phi1+b,"/",GaussHyp1F2a,"/",phi1+c+1,"/",phi1+c,"/",GaussHyp1F2b,"=",MEr1RhoPLFf,"\n")
  return(MEr1RhoPLF)}

################################################################################
# MEr1LqPLFf implements Bayesian estimation under PLF for Lq in M/Er/1 queues
################################################################################
MEr1LqPLFf<-function(samp,r,a,b,c,z) {
  GaussHypItg<-function(u,a,b,c,z) {
    U<-u^(b-1)*(1-u)^(c-b-1)*(1-u*z)^(-a)
    return(U)}
  n<-length(samp)
  infty<-max(samp)
  nij<-matrix(0,infty+1,infty+1)
  phi1<-0.0
  phi2<-0.0
  i<-2
  while (i<=n) {
    #    cat("i=",i,"\n")
    nij[samp[i-1]+1,samp[i]+1]<-nij[samp[i-1]+1,samp[i]+1]+1
    i<-i+1}
  j=0
  while (j<=infty) {
    #    cat("j=",j,"\n")
    phi1<-phi1+j*nij[0+1,j+1]
    phi2<-phi2+(r+j)*nij[0+1,j+1]
    j=j+1}
  i=1
  while (i<=infty) {
    j<-i-1
    while (j<=infty) {
      #      cat("(i,j)=(",i,",",j,")\n")
      phi1<-phi1+(j-i+1)*nij[i+1,j+1]
      phi2<-phi2+(r+j-i+1)*nij[i+1,j+1]
      j<-j+1}
    i<-i+1}
  #  cat("phi1=",phi1,"\t","phi2=",phi2,"\n")
  GaussHyp1F2a<-integrate(GaussHypItg,0,1,phi2+a,phi1+b+4,phi1+c+2,z)$value/beta(phi1+b+4,c-b-2)
  GaussHyp1F2b<-integrate(GaussHypItg,0,1,phi2+a,phi1+b  ,phi1+c  ,z)$value/beta(phi1+b  ,c-b  )
  MEr1LqPLF<-(r+1)/(2*r)*sqrt((phi1+b+3)*(phi1+b+2)*(phi1+b+1)*(phi1+b)*GaussHyp1F2a/(phi1+c+1)/(phi1+c)/(c-b-1)/(c-b-2)/GaussHyp1F2b)
  # cat("MEr1LqPLF=", MEr1LqPLF,"\n")
  return(MEr1LqPLF)}

################################################################################
# predictive distributions
################################################################################
MEr1Predf<-function(x,samp,r,a,b,c,z) {
  GaussHypItg<-function(u,a,b,c,z) {
    U<-u^(b-1)*(1-u)^(c-b-1)*(1-u*z)^(-a)
    return(U)}
  n<-length(samp)
  infty<-max(samp)
  nij<-matrix(0,infty+1,infty+1)
  phi1<-0.0
  phi2<-0.0
  i<-2
  while (i<=n) {
    #    cat("i=",i,"\n")
    nij[samp[i-1]+1,samp[i]+1]<-nij[samp[i-1]+1,samp[i]+1]+1
    i<-i+1}
  j=0
  while (j<=infty) {
    #    cat("j=",j,"\n")
    phi1<-phi1+j*nij[0+1,j+1]
    phi2<-phi2+(r+j)*nij[0+1,j+1]
    j=j+1}
  i=1
  while (i<=infty) {
    j<-i-1
    while (j<=infty) {
      #      cat("(i,j)=(",i,",",j,")\n")
      phi1<-phi1+(j-i+1)*nij[i+1,j+1]
      phi2<-phi2+(r+j-i+1)*nij[i+1,j+1]
      j<-j+1}
    i<-i+1}
  # cat("phi1=",phi1,"\t","phi2=",phi2,"\n")
  MEr1Predf<-numeric(length(x))
  for (i in 1:length(x)){
    GaussHyp1F2a<-integrate(GaussHypItg,0,1,phi2+a+x[i]+r,phi1+b+x[i],phi1+c+x[i],z)[[1]]/beta(phi1+b+x[i],c-b)
    GaussHyp1F2b<-integrate(GaussHypItg,0,1,phi2+a,phi1+b,phi1+c,z)[[1]]/beta(phi1+b,c-b)
    MEr1Predf[i]<-choose(x[i]+r-1,r-1)*r^{-x[i]}*beta(phi1+b+x[i],c-b)*GaussHyp1F2a/beta(phi1+b,c-b)/GaussHyp1F2b
    # cat(MEr1Predf[i],"\n")
  }
  return(MEr1Predf)}

################################################################################
# Bayesian credible region for rho in M/Er/1 queues
################################################################################
MEr1CR<-function(samp,r,a,b,c,z){
  epsilon<-1e-04
  sl<-0.05  # significance level
  accum <- 0.0
  x<-seq(0,1,by=epsilon)
  dens<-dGaussHypPostf(x,samp,r,a,b,c,z)
  # find inferior limit
  i<-1
  while (accum<(sl/2) && i<length(dens)){
    i<-i+1
    accum<-accum+epsilon*(dens[i-1]+dens[i])/2}
  linf<-x[i-1]
  # find superior limit
  while (accum<(1-sl/2) && i<length(dens)){
    i<-i+1
    accum<-accum+epsilon*(dens[i-1]+dens[i])/2}
  lsup<-x[i]
  return(list(linf=linf,lsup=lsup))}

################################################################################
# Monte Carlo simulation for rho MLE state-of-the-art
################################################################################
MEr1MonteCarloRhoSoA<-function(size,rho,r){
  set.seed(2026)
  mcRep<-1000
  est<-numeric(mcRep)
  for(i in 1:mcRep){
    # cat("MEr1MonteCarloSoA:i=",i,"\n")
    smp<-rnbinom(size,r,1/(1+rho/r))
    est[i]<-MEr1RhoMLESoA(smp,r)
    # cat("MEr1MonteCarloSoA: estimate=",est[i],"\n")
  }
  mest<-mean(est); vest<-var(est)
  return(c(mest,vest,mest-rho,sqrt(vest+(mest-rho)^2)))}

################################################################################
# Monte Carlo simulation for rho
################################################################################
MEr1MonteCarloRho<-function(size,rho,r,fEst,...){
  set.seed(2026)
  mcRep<-1000
  est<-numeric(mcRep)
  for(i in 1:mcRep){
    # cat("MEr1MonteCarlo:i=",i,"\n")
    smp<-rMEr1(size,rho,r)
    est[i]<-fEst(smp,r,...)
    # cat("MEr1MonteCarlo: estimate=",est[i],"\n")
  }
  mest<-mean(est); vest<-var(est)
  return(c(mest,vest,mest-rho,sqrt(vest+(mest-rho)^2)))}

################################################################################
# Monte Carlo simulation for Lq
################################################################################
MEr1MonteCarloLq<-function(size,Lq,r,fEst,...){
  set.seed(2026)
  mcRep<-1000
  est<-numeric(mcRep)
  for(i in 1:mcRep){
    # cat("MEr1MonteCarlo:i=",i,"\n")
    rho<-rhoLqMEr1(Lq,r)
    smp<-rMEr1(size,rho,r)
    est[i]<-fEst(smp,r,...)
    # cat("MEr1MonteCarlo: estimate=",est[i],"\n")
  }
  mest<-mean(est); vest<-var(est)
  return(c(mest,vest,mest-Lq,sqrt(vest+(mest-Lq)^2)))}

################################################################################
# Monte Carlo simulation for credible regions
################################################################################
MEr1MonteCarloCR<-function(size,rho,r,a,b,c){
  set.seed(2026)
  mcRep<-1000
  crcover<-0
  crlength<-numeric(length(size))
  est<-numeric(mcRep)
  for(i in 1:mcRep){
    # cat("MEr1MonteCarloCR:i=",i,"\n")
    smp<-rMEr1(size,rho,r)
    cr<-MEr1CR(smp,r,a,b,c,-1/r)
    # cat("MEr1MonteCarloCR: ci=",cr$linf,cr$lsup,"\n")
    if ((cr$linf<=rho)&&(rho<=cr$lsup))
      crcover=crcover+1
    end
    crlength[i]=cr$lsup-cr$linf
  }
  crmean<-mean(crlength)
  crcover<-crcover/mcRep
  return(c(crmean=crmean,crcover=crcover))}

################################################################################
# Monte Carlo table for rho MLE state-of-the-art
################################################################################
MEr1MonteCarloRhoSoATab<-function(sizes,rhos,r){
  tab<-matrix(nrow=length(rhos),ncol=1+4*length(sizes))
  for (i in 1:length(rhos)){
    tab[i,1]<-rhos[i]
    for (j in 1:length(sizes)){
      est<-MEr1MonteCarloRhoSoA(sizes[j],rhos[i],r)
      tab[i,(4*(j-1)+2):(4*(j-1)+5)]<-est
    }
  }
  return(tab)}

################################################################################
# Monte Carlo table for rho
################################################################################
MEr1MonteCarloRhoTab<-function(sizes,rhos,r,fEst,...){
  tab<-matrix(nrow=length(rhos),ncol=1+4*length(sizes))
  for (i in 1:length(rhos)){
    tab[i,1]<-rhos[i]
    for (j in 1:length(sizes)){
      est<-MEr1MonteCarloRho(sizes[j],rhos[i],r,fEst,...)
      tab[i,(4*(j-1)+2):(4*(j-1)+5)]<-est
    }
  }
  return(tab)}

################################################################################
# Monte Carlo table for Lq
################################################################################
MEr1MonteCarloLqTab<-function(sizes,Lqs,r,fEst,...){
  tab<-matrix(nrow=length(Lqs),ncol=1+4*length(sizes))
  for (i in 1:length(Lqs)){
    tab[i,1]<-Lqs[i]
    for (j in 1:length(sizes)){
      est<-MEr1MonteCarloLq(sizes[j],Lqs[i],r,fEst,...)
      tab[i,(4*(j-1)+2):(4*(j-1)+5)]<-est
    }
  }
  return(tab)}

################################################################################
# Monte Carlo table for confidence region
################################################################################
MEr1MonteCarloCRTab<-function(sizes,rhos,r,a,b,c){
  tab<-matrix(nrow=length(rhos),ncol=1+2*length(sizes))
  for (i in 1:length(rhos)){
    tab[i,1]<-rhos[i]
    for (j in 1:length(sizes)){
      est<-MEr1MonteCarloCR(sizes[j],rhos[i],r,a,b,c)
      tab[i,(2*(j-1)+2):(2*(j-1)+3)]<-est
    }
  }
  return(tab)}

################################################################################
# experiments
################################################################################

rem_var()
################################################################################
# plotting prior Gaussian hypergeometric distributions
################################################################################
{setEPS()
postscript(paste("FiGaussHypR.eps",sep=""),width=10.5*1.0,height=8*1.0)
# dev.new(width=10.5*0.5,height=8*0.5)
par(mfrow=c(1,1))
# rem_var()
a<-1
b<-2
c<-6
r<-2
gammap<-a
alphap<-b
betap <-c-b
dGaussHypf(0.5,gammap,alphap,alphap+betap,-1/r)
integrate(dGaussHypf,0,1,gammap,alphap,alphap+betap,-1/r)[[1]]
rho<-seq(0.01,0.99,0.01)
plot(rho,dGaussHypf(rho,gammap,alphap,alphap+betap,-1/r),type="l",lty=1,col=1,
     lwd=2,xlab=expression(rho),ylab=expression(pi[1](rho)))
a<-0
b<-1
c<-2
r<-2
gammap<-a
alphap<-b
betap <-c-b
integrate(dGaussHypf,0,1,gammap,alphap,alphap+betap,-1/r)[[1]]
rho<-seq(0.01,0.99,0.01)
lines(rho,dGaussHypf(rho,gammap,alphap,alphap+betap,-1/r),type="l",lty=2,col=2,
      lwd=2)
a<-1
b<-4
c<-6
r<-2
gammap<-a
alphap<-b
betap <-c-b
integrate(dGaussHypf,0,1,gammap,alphap,alphap+betap,-1/r)[[1]]
lines(rho,dGaussHypf(rho,gammap,alphap,alphap+betap,-1/r),type="l",lty=3,col=3,
      lwd=2)
legend("bottom",lty=c(1,2,3),col=c(1,2,3),lwd=2, legend=c(
  expression(paste("GH(",gamma,"=1;",alpha,"=2;",alpha+beta,"=6;-1/r)")),
  expression(paste("GH(",gamma,"=0;",alpha,"=1;",alpha+beta,"=2;-1/r)")),
  expression(paste("GH(",gamma,"=1;",alpha,"=4;",alpha+beta,"=6;-1/r)"))))
graphics.off()
}

################################################################################
# plotting prior beta distributions
################################################################################
{setEPS()
postscript(paste("FiBetaR.eps",sep=""),width=10.5*1.0,height=8*1.0)
# dev.new(width=10.5*0.5,height=8*0.5)
par(mfrow=c(1,1))
# rem_var()
gammap<-0
alphap<-2
betap <-4
a<-gammap
b<-alphap
c<-alphap+betap
r<-0
dGaussHypf(0.5,a,b,c,-1/r)
integrate(dGaussHypf,0,1,a,b,c,-1/r)[[1]]
rho<-seq(0.01,0.99,0.01)
plot(rho,dGaussHypf(rho,a,b,c,-1/r),type="l",lty=1,col=1,
     lwd=2,xlab=expression(rho),ylab=expression(pi[2](rho)))
gammap<-0
alphap<-1
betap <-1
a<-gammap
b<-alphap
c<-alphap+betap
r<-0
integrate(dGaussHypf,0,1,a,b,c,-1/r)[[1]]
lines(rho,dGaussHypf(rho,a,b,c,-1/r),type="l",lty=2,col=2,
      lwd=2)
gammap<-0
alphap<-4
betap <-2
a<-gammap
b<-alphap
c<-alphap+betap
r<-0
integrate(dGaussHypf,0,1,a,b,c,-1/r)[[1]]
lines(rho,dGaussHypf(rho,a,b,c,-1/r),type="l",lty=3,col=3,
      lwd=2)
legend("bottom",lty=c(1,2,3),col=c(1,2,3),lwd=2, legend=c(
  expression(paste("Beta(",alpha,"=2;",beta,"=4)")),
  expression(paste("Beta(",alpha,"=1;",beta,"=1)")),
  expression(paste("Beta(",alpha,"=4;",beta,"=2)"))))
graphics.off()
}

################################################################################
# Monte Carlo tables for rho r = 1
################################################################################
{# rem_var()
  # Monte Carlo table for rho under MLE state-of-the-art r = 1
  r<-1
  rhos<-c(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
  sizes<-c(10,20,50,100,200)
  MEr1RhoMLESoA1<-MEr1MonteCarloRhoSoATab(sizes,rhos,r)
  
  # Monte Carlo table for rho under MLE r = 1
  r<-1
  rhos<-c(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
  sizes<-c(10,20,50,100,200)
  MEr1RhoMLE1<-MEr1MonteCarloRhoTab(sizes,rhos,r,MEr1RhoMLEf)
  
  # Monte Carlo table for rho under HG-SELF r = 1
  # rem_var()
  r<-1
  rhos<-c(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
  sizes<-c(10,20,50,100,200)
  a<-1
  b<-2
  c<-6
  gammap<-a
  alphap<-b
  betap <-c-b
  MEr1RhoHGSELF1<-MEr1MonteCarloRhoTab(sizes,rhos,r,MEr1RhoSELFf,gammap,alphap,alphap+betap,-1/r)
  
  # Monte Carlo table for rho under HG-PLF r = 1
  # rem_var()
  r<-1
  rhos<-c(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
  sizes<-c(10,20,50,100,200)
  a<-1
  b<-2
  c<-6
  gammap<-a
  alphap<-b
  betap <-c-b
  MEr1RhoHGPLF1<-MEr1MonteCarloRhoTab(sizes,rhos,r,MEr1RhoPLFf,gammap,alphap,alphap+betap,-1/r)
  
  # Monte Carlo table for rho under B-SELF r = 1
  # rem_var()
  r<-1
  rhos<-c(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
  sizes<-c(10,20,50,100,200)
  gammap<-0
  alphap<-2
  betap <-4
  a<-gammap
  b<-alphap
  c<-alphap+betap
  MEr1RhoBSELF1<-MEr1MonteCarloRhoTab(sizes,rhos,r,MEr1RhoSELFf,a,b,c,-1/r)
  
  # Monte Carlo table for rho under B-PLF r = 1
  # rem_var()
  r<-1
  rhos<-c(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
  sizes<-c(10,20,50,100,200)
  gammap<-0
  alphap<-2
  betap <-4
  a<-gammap
  b<-alphap
  c<-alphap+betap
  MEr1RhoBPLF1<-MEr1MonteCarloRhoTab(sizes,rhos,r,MEr1RhoPLFf,a,b,c,-1/r)
}

# save results
MEr1RhoRes1<-rbind(MEr1RhoMLESoA1,MEr1RhoMLE1,MEr1RhoHGSELF1,
                   MEr1RhoHGPLF1,MEr1RhoBSELF1,MEr1RhoBPLF1)
save(MEr1RhoRes1,file='MEr1RhoRes1.rdata')

# load results
# load(file='MEr1RhoRes1.rdata')

# export results to Excel
ExExcel(MEr1RhoRes1)

################################################################################
# Monte Carlo tables for rho r = 4
################################################################################
{# rem_var()
  # Monte Carlo for rho under MLE state-of-the-art r = 4
  r<-4
  rhos<-c(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
  sizes<-c(10,20,50,100,200)
  MEr1RhoMLESoA4<-MEr1MonteCarloRhoSoATab(sizes,rhos,r)
  
  # Monte Carlo table for rho under MLE r = 4
  r<-4
  rhos<-c(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
  sizes<-c(10,20,50,100,200)
  MEr1RhoMLE4<-MEr1MonteCarloRhoTab(sizes,rhos,r,MEr1RhoMLEf)
  
  # Monte Carlo table for rho under HG-SELF r = 4
  # rem_var()
  r<-4
  rhos<-c(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
  sizes<-c(10,20,50,100,200)
  a<-1
  b<-2
  c<-6
  gammap<-a
  alphap<-b
  betap <-c-b
  MEr1RhoHGSELF4<-MEr1MonteCarloRhoTab(sizes,rhos,r,MEr1RhoSELFf,gammap,alphap,alphap+betap,-1/r)
  
  # Monte Carlo table for rho under HG-PLF r = 4
  # rem_var()
  r<-4
  rhos<-c(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
  sizes<-c(10,20,50,100,200)
  a<-1
  b<-2
  c<-6
  gammap<-a
  alphap<-b
  betap <-c-b
  MEr1RhoHGPLF4<-MEr1MonteCarloRhoTab(sizes,rhos,r,MEr1RhoPLFf,gammap,alphap,alphap+betap,-1/r)
  
  # Monte Carlo table for rho under B-SELF r = 4
  # rem_var()
  r<-4
  rhos<-c(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
  sizes<-c(10,20,50,100,200)
  gammap<-0
  alphap<-2
  betap <-4
  a<-gammap
  b<-alphap
  c<-alphap+betap
  MEr1RhoBSELF4<-MEr1MonteCarloRhoTab(sizes,rhos,r,MEr1RhoSELFf,a,b,c,-1/r)
  
  # Monte Carlo table for rho under B-PLF r = 4
  # rem_var()
  r<-4
  rhos<-c(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
  sizes<-c(10,20,50,100,200)
  gammap<-0
  alphap<-2
  betap <-4
  a<-gammap
  b<-alphap
  c<-alphap+betap
  MEr1RhoBPLF4<-MEr1MonteCarloRhoTab(sizes,rhos,r,MEr1RhoPLFf,a,b,c,-1/r)
}

# save results
MEr1RhoRes4<-rbind(MEr1RhoMLESoA4,MEr1RhoMLE4,MEr1RhoHGSELF4,
                   MEr1RhoHGPLF4,MEr1RhoBSELF4,MEr1RhoBPLF4)
save(MEr1RhoRes1,file='MEr1RhoRes4.rdata')

# load results
# load(file='MEr1RhoRes4.rdata')

# export results to Excel
ExExcel(MEr1RhoRes4)

################################################################################
# Monte Carlo tables for Lq r = 1
################################################################################
{# rem_var()
  # Monte Carlo table for Lq under MLE r = 1
  r<-1
  Lqs<-c(0.50,1.0,2.0,4.0)
  sizes<-c(10,20,50,100,200)
  MEr1LqMLE1<-MEr1MonteCarloLqTab(sizes,Lqs,r,MEr1LqMLEf,)
  
  # Monte Carlo table for Lq under HG-SELF r = 1
  # rem_var()
  r<-1
  Lqs<-c(0.50,1.0,2.0,4.0)
  sizes<-c(10,20,50,100,200)
  a<-1
  b<-2
  c<-6
  gammap<-a
  alphap<-b
  betap <-c-b
  MEr1LqHGSELF1<-MEr1MonteCarloLqTab(sizes,Lqs,r,MEr1LqSELFf,gammap,alphap,alphap+betap,-1/r)
  
  # Monte Carlo table for Lq under HG-PLF r = 1
  # rem_var()
  r<-1
  Lqs<-c(0.50,1.0,2.0,4.0)
  sizes<-c(10,20,50,100,200)
  rhos<-sapply(Lqs,rhoLqMEr1,r)
  a<-1
  b<-2
  c<-6
  gammap<-a
  alphap<-b
  betap <-c-b
  MEr1LqHGPLF1<-MEr1MonteCarloLqTab(sizes,Lqs,r,MEr1LqPLFf,gammap,alphap,alphap+betap,-1/r)
  
  # Monte Carlo table for Lq under B-SELF r = 1
  # rem_var()
  r<-1
  Lqs<-c(0.50,1.0,2.0,4.0)
  sizes<-c(10,20,50,100,200)
  rhos<-sapply(Lqs,rhoLqMEr1,r)
  gammap<-0
  alphap<-2
  betap <-4
  a<-gammap
  b<-alphap
  c<-alphap+betap
  MEr1LqBSELF1<-MEr1MonteCarloLqTab(sizes,Lqs,r,MEr1LqSELFf,a,b,c,-1/r)
  
  # Monte Carlo table for Lq under B-PLF r = 1
  # rem_var()
  r<-1
  Lqs<-c(0.50,1.0,2.0,4.0)
  sizes<-c(10,20,50,100,200)
  rhos<-sapply(Lqs,rhoLqMEr1,r)
  gammap<-0
  alphap<-2
  betap <-4
  a<-gammap
  b<-alphap
  c<-alphap+betap
  MEr1LqBPLF1<-MEr1MonteCarloLqTab(sizes,Lqs,r,MEr1LqPLFf,a,b,c,-1/r)
}

# save results
MEr1LqRes1<-rbind(MEr1LqMLE1,MEr1LqHGSELF1,MEr1LqHGPLF1,
                  MEr1LqBSELF1,MEr1LqBPLF1)
save(MEr1RhoRes1,file='MEr1LqRes1.rdata')

# load results
# load(file='MEr1LqRes1.rdata')

# export results to Excel
ExExcel(MEr1LqRes1)

################################################################################
# Monte Carlo tables for Lq r = 4
################################################################################
{# rem_var()
  # Monte Carlo table for Lq under MLE r = 4
  r<-4
  Lqs<-c(0.50,1.0,2.0,4.0)
  sizes<-c(10,20,50,100,200)
  MEr1LqMLE4<-MEr1MonteCarloLqTab(sizes,Lqs,r,MEr1LqMLEf,)
  
  # Monte Carlo table for Lq under HG-SELF r = 4
  # rem_var()
  r<-4
  Lqs<-c(0.50,1.0,2.0,4.0)
  sizes<-c(10,20,50,100,200)
  a<-1
  b<-2
  c<-6
  gammap<-a
  alphap<-b
  betap <-c-b
  MEr1LqHGSELF4<-MEr1MonteCarloLqTab(sizes,Lqs,r,MEr1LqSELFf,gammap,alphap,alphap+betap,-1/r)
  
  # Monte Carlo table for Lq under HG-PLF r = 4
  # rem_var()
  r<-4
  Lqs<-c(0.50,1.0,2.0,4.0)
  sizes<-c(10,20,50,100,200)
  rhos<-sapply(Lqs,rhoLqMEr1,r)
  a<-1
  b<-2
  c<-6
  gammap<-a
  alphap<-b
  betap <-c-b
  MEr1LqHGPLF4<-MEr1MonteCarloLqTab(sizes,Lqs,r,MEr1LqPLFf,gammap,alphap,alphap+betap,-1/r)
  
  # Monte Carlo table for Lq under B-SELF r = 4
  # rem_var()
  r<-4
  Lqs<-c(0.50,1.0,2.0,4.0)
  sizes<-c(10,20,50,100,200)
  gammap<-0
  alphap<-2
  betap <-4
  a<-gammap
  b<-alphap
  c<-alphap+betap
  MEr1LqBSELF4<-MEr1MonteCarloLqTab(sizes,Lqs,r,MEr1LqSELFf,a,b,c,-1/r)
  
  # Monte Carlo table for Lq under B-PLF r = 4
  # rem_var()
  r<-4
  Lqs<-c(0.50,1.0,2.0,4.0)
  sizes<-c(10,20,50,100,200)
  gammap<-0
  alphap<-2
  betap <-4
  a<-gammap
  b<-alphap
  c<-alphap+betap
  MEr1LqBPLF4<-MEr1MonteCarloLqTab(sizes,Lqs,r,MEr1LqPLFf,a,b,c,-1/r)
}

# save results
MEr1LqRes4<-rbind(MEr1LqMLE4,MEr1LqHGSELF4,MEr1LqHGPLF4,
                  MEr1LqBSELF4,MEr1LqBPLF4)
save(MEr1RhoRes4,file='MEr1LqRes4.rdata')

# load results
# load(file='MEr1LqRes4.rdata')

# export results to Excel
ExExcel(MEr1LqRes4)

################################################################################
# computing predictive distribution for a prior Gaussian hypergeometric r = 1
################################################################################
{# rem_var()
  a<-1
  b<-2
  c<-6
  r<-1
  gammap<-a
  alphap<-b
  betap <-c-b
  x<-seq(0,10,1)
  set.seed(2026)
  smer1<-rMEr1(n=20,rho=0.50,r)
  PAxGH1<-MEr1Predf(x,smer1,r,gammap,alphap,alphap+betap,-1/r)
  #
  set.seed(2026)
  smer1<-rMEr1(n=50,rho=0.50,r)
  PAxGH1<-cbind(PAxGH1,MEr1Predf(x,smer1,r,gammap,alphap,alphap+betap,-1/r))
  #
  set.seed(2026)
  smer1<-rMEr1(n=100,rho=0.50,r)
  PAxGH1<-cbind(PAxGH1,MEr1Predf(x,smer1,r,gammap,alphap,alphap+betap,-1/r))
  #
  set.seed(2026)
  smer1<-rMEr1(n=200,rho=0.50,r)
  PAxGH1<-cbind(PAxGH1,MEr1Predf(x,smer1,r,gammap,alphap,alphap+betap,-1/r))
}

################################################################################
# computing predictive distribution for a prior Gaussian hypergeometric r = 4
################################################################################
{# rem_var()
  a<-1
  b<-2
  c<-6
  r<-4
  gammap<-a
  alphap<-b
  betap <-c-b
  x<-seq(0,10,1)
  set.seed(2026)
  smer1<-rMEr1(n=20,rho=0.50,r)
  PAxGH4<-MEr1Predf(x,smer1,r,gammap,alphap,alphap+betap,-1/r)
  #
  set.seed(2026)
  smer1<-rMEr1(n=50,rho=0.50,r)
  PAxGH4<-cbind(PAxGH4,MEr1Predf(x,smer1,r,gammap,alphap,alphap+betap,-1/r))
  #
  set.seed(2026)
  smer1<-rMEr1(n=100,rho=0.50,r)
  PAxGH4<-cbind(PAxGH4,MEr1Predf(x,smer1,r,gammap,alphap,alphap+betap,-1/r))
  #
  set.seed(2026)
  smer1<-rMEr1(n=200,rho=0.50,r)
  PAxGH4<-cbind(PAxGH4,MEr1Predf(x,smer1,r,gammap,alphap,alphap+betap,-1/r))
}

################################################################################
# computing predictive distribution for a prior beta r = 1
################################################################################
{# rem_var()
  gammap<-0
  alphap<-2
  betap <-4
  a<-gammap
  b<-alphap
  c<-alphap+betap
  r<-1
  x<-seq(0,10,1)
  set.seed(2026)
  smer1<-rMEr1(n=20,rho=0.50,r)
  PAxB1<-MEr1Predf(x,smer1,r,gammap,alphap,alphap+betap,-1/r)
  #
  set.seed(2026)
  smer1<-rMEr1(n=50,rho=0.50,r)
  PAxB1<-cbind(PAxB1,MEr1Predf(x,smer1,r,gammap,alphap,alphap+betap,-1/r))
  #
  set.seed(2026)
  smer1<-rMEr1(n=100,rho=0.50,r)
  PAxB1<-cbind(PAxB1,MEr1Predf(x,smer1,r,gammap,alphap,alphap+betap,-1/r))
  #
  set.seed(2026)
  smer1<-rMEr1(n=200,rho=0.50,r)
  PAxB1<-cbind(PAxB1,MEr1Predf(x,smer1,r,gammap,alphap,alphap+betap,-1/r))
}

################################################################################
# computing predictive distribution for a prior beta r = 4
################################################################################
{# rem_var()
  gammap<-0
  alphap<-2
  betap <-4
  a<-gammap
  b<-alphap
  c<-alphap+betap
  r<-4
  x<-seq(0,10,1)
  set.seed(2026)
  smer1<-rMEr1(n=20,rho=0.50,r)
  PAxB4<-MEr1Predf(x,smer1,r,gammap,alphap,alphap+betap,-1/r)
  #
  set.seed(2026)
  smer1<-rMEr1(n=50,rho=0.50,r)
  PAxB4<-cbind(PAxB4,MEr1Predf(x,smer1,r,gammap,alphap,alphap+betap,-1/r))
  #
  set.seed(2026)
  smer1<-rMEr1(n=100,rho=0.50,r)
  PAxB4<-cbind(PAxB4,MEr1Predf(x,smer1,r,gammap,alphap,alphap+betap,-1/r))
  #
  set.seed(2026)
  smer1<-rMEr1(n=200,rho=0.50,r)
  PAxB4<-cbind(PAxB4,MEr1Predf(x,smer1,r,gammap,alphap,alphap+betap,-1/r))
}

# export results to Excel
ExExcel(rbind(cbind(PAxGH1,PAxGH4),cbind(PAxB1,PAxB4)))

################################################################################
# Monte Carlo tables for credible region
################################################################################
{# rem_var()
  # Monte Carlo table for CR for rho under HG-SELF r = 1
  # rem_var()
  r<-1
  rhos<-c(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
  sizes<-c(10,20,50,100,200)
  a<-1
  b<-2
  c<-6
  gammap<-a
  alphap<-b
  betap <-c-b
  MEr1CRHGSELF1<-MEr1MonteCarloCRTab(sizes,rhos,r,gammap,alphap,alphap+betap)
  
  # Monte Carlo table for CR for rho under B-SELF r = 1
  # rem_var()
  r<-1
  rhos<-c(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
  sizes<-c(10,20,50,100,200)
  gammap<-0
  alphap<-2
  betap <-4
  a<-gammap
  b<-alphap
  c<-alphap+betap
  MEr1CRBSELF1<-MEr1MonteCarloCRTab(sizes,rhos,r,a,b,c)
}

# save results
MEr1CRRes1<-rbind(MEr1CRHGSELF1,MEr1CRBSELF1)
save(MEr1CRRes1,file='MEr1CRRes1.rdata')

# load results
# load(file='MEr1CRRes1')

# export results to Excel
ExExcel(MEr1CRRes1)

{# rem_var()
  # Monte Carlo table for CR for rho under HG-SELF r = 4
  # rem_var()
  r<-4
  rhos<-c(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
  sizes<-c(10,20,50,100,200)
  a<-1
  b<-2
  c<-6
  gammap<-a
  alphap<-b
  betap <-c-b
  MEr1CRHGSELF4<-MEr1MonteCarloCRTab(sizes,rhos,r,gammap,alphap,alphap+betap)
  
  # Monte Carlo table for CR for rho under B-SELF r = 4
  # rem_var()
  r<-4
  rhos<-c(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
  sizes<-c(10,20,50,100,200)
  gammap<-0
  alphap<-2
  betap <-4
  a<-gammap
  b<-alphap
  c<-alphap+betap
  MEr1CRBSELF4<-MEr1MonteCarloCRTab(sizes,rhos,r,a,b,c)
}

# save results
MEr1CRRes4<-rbind(MEr1CRHGSELF4,MEr1CRBSELF4)
save(MEr1CRRes4,file='MEr1CRRes4.rdata')

# load results
# load(file='MEr1CRRes4')

# export results to Excel
ExExcel(MEr1CRRes4)

stop("Simulation concluded")

################################################################################
# tests
################################################################################

# testing discrete event simulation for a M/Er/1 queue
# rem_var()
par(mfrow=c(1,2))
res1<-rMEr1(n=100,rho=0.5,r=1)
table(res1)
res2<-rgeom(n=100,p=(1-0.5))
table(res2)
plot(table(res1),type="h")
plot(table(res2),type="h")
(rMEr1(n=10,rho=0.50,r=4))
(rMEr1(n=10,rho=0.99,r=4))

# testing function for traffic intensity, rho, for a given Lq or L
# rem_var()
r<-4
rho<-0.9
for (rho in seq(0.001,0.999,0.05)) {
  Lq<-LqMEr1(rho,r)
  cat("Lq(rho=", rho, ")=", Lq, "\n")
  rhoEst<-rhoLqMEr1(Lq,r)
  cat("rho(Lq=", Lq, ")=", rhoEst, "\n")
  L<-LMEr1(rho,r)
  cat("L(rho=", rho, ")=", L, "\n")
  rhoEst<-rhoLMEr1(L,r)
  cat("rho(L=", L, ")=", rhoEst, "\n")}

# testing posterior distribution
# dev.new(width=10.5*0.5,height=8*0.5)
par(mfrow=c(1,1))
r<-1
# rem_var()
a<-1
b<-2
c<-6
gammap<-a
alphap<-b
betap <-c-b
rho<-0.50
set.seed(2026)
smer1<-rMEr1(n=200,rho,r)
dGaussHypPostf(0.5,smer1,r,gammap,alphap,alphap+betap,-1/r)
integrate(dGaussHypPostf,0,1,smer1,r,gammap,alphap,alphap+betap,-1/r)[[1]]
rho<-seq(0,1,0.01)
plot(rho,dGaussHypPostf(rho,smer1,r,gammap,alphap,alphap+betap,-1/r),
     type="l",lty=1,col=2,lwd=2,xlab=expression(rho),ylab=expression(pi[1](rho)))

# testing MLE for rho (state of the art)
# rem_var()
r<-4
rho<-0.50
set.seed(2026)
smer12<-rnbinom(n=200,r,1/(1+rho/r))
(MEr1RhoMLESoA(smer12,r))

# testing MLE for rho
# rem_var()
r<-4
rho<-0.50
set.seed(2026)
smer1<-rMEr1(n=200,rho,r)
(MEr1RhoMLEf(smer1,r))

# testing MLE for Lq
# rem_var()
r<-4
Lq<-1.0
rho<-rhoLqMEr1(Lq,r)
set.seed(2026)
smer1<-rMEr1(n=200,rho,r)
(MEr1LqMLEf(smer1,r))

# testing Bayes estimator under squared error loss function (SELF) for rho 
# rem_var()
r<-4
rho<-0.50
a<-1
b<-2
c<-6
gammap<-a
alphap<-b
betap <-c-b
set.seed(2026)
smer1<-rMEr1(n=200,rho,r)
(MEr1RhoSELF<-MEr1RhoSELFf(smer1,r,gammap,alphap,alphap+betap,-1/r))
# numerical checking
mer1ErhoItg<-function(x,samp,r,a,b,c,z) return(x*dGaussHypPostf(x,samp,r,a,b,c,z))
(MEr1RhoSELFChk<-integrate(mer1ErhoItg,0,1,smer1,r,gammap,alphap,alphap+betap,-1/r)$value)
cat("Error=", MEr1RhoSELF-MEr1RhoSELFChk,"\n")

# testing Bayes estimates under squared error loss function (SELF) for Lq 
# rem_var()
r<-4
a<-1
b<-2
c<-6
gammap<-a
alphap<-b
betap <-c-b
Lq=1.0
rho<-rhoLqMEr1(Lq,r)
set.seed(2026)
smer1<-rMEr1(n=200,rho,r)
(MEr1LqSELF<-MEr1LqSELFf(smer1,r,gammap,alphap,alphap+betap,-1/r))
# numerical checking
mer1ErhoItg<-function(x,samp,r,a,b,c,z)
  return(((r+1)*x^2/(2*r)/(1-x))*dGaussHypPostf(x,samp,r,a,b,c,z))
(MEr1LqSELFChk<-integrate(mer1ErhoItg,0,1,smer1,r,gammap,alphap,alphap+betap,-1/r)$value)
cat("Error=", MEr1LqSELF-MEr1LqSELFChk,"\n")

# testing Bayes estimator under precautionary loss function (PLF) for rho
r<-4
rho<-0.50
a<-1
b<-2
c<-6
gammap<-a
alphap<-b
betap <-c-b
set.seed(2026)
smer1<-rMEr1(n=200,rho,r)
(MEr1RhoPLF<-MEr1RhoPLFf(smer1,r,gammap,alphap,alphap+betap,-1/r))
# numerical checking
mer1ErhoItg<-function(x,samp,r,a,b,c,z)
              return(x^2*dGaussHypPostf(x,samp,r,a,b,c,z))
(MEr1RhoPLFChk<-sqrt(integrate(mer1ErhoItg,0,1,smer1,r,gammap,alphap,alphap+betap,-1/r)$value))
cat("Error=", MEr1RhoPLF-MEr1RhoPLFChk,"\n")

# testing Bayes estimates for Lq under precautionary loss function (PLF)
# rem_var()
r<-4
a<-1
b<-2
c<-6
gammap<-a
alphap<-b
betap <-c-b
Lq=1.0
rho<-rhoLqMEr1(Lq,r)
set.seed(2026)
smer1<-rMEr1(n=200,rho,r)
(MEr1LqPLF<-MEr1LqPLFf(smer1,r,gammap,alphap,alphap+betap,-1/r))
# numerical checking
mer1ErhoItg<-function(x,samp,r,a,b,c,z)
  return(((r+1)*x^2/(2*r)/(1-x))^2*dGaussHypPostf(x,samp,r,a,b,c,z))
(MEr1LqPLFChk<-sqrt(integrate(mer1ErhoItg,0,1,smer1,r,gammap,alphap,alphap+betap,-1/r)$value))
cat("Error=", MEr1LqPLF-MEr1LqPLFChk,"\n")

# testing Monte Carlo functions
# MLE state-of-the-art
# rem_var()
r<-1
rho<-0.50
size<-100
MEr1MonteCarloRhoSoA(size,rho,r)

# MLE
# rem_var()
r<-1
rho<-0.50
size<-100
MEr1MonteCarloRho(size,rho,r,MEr1RhoMLEf)

# Gauss hypergeometric prior
r<-1
rho<-0.50
size<-100
a<-1
b<-2
c<-6
gammap<-a
alphap<-b
betap <-c-b
MEr1MonteCarloRho(size,rho,r,MEr1RhoSELFf,gammap,alphap,alphap+betap,-1/r)
MEr1MonteCarloRho(size,rho,r,MEr1RhoPLFf,gammap,alphap,alphap+betap,-1/r)

# beta prior
r<-1
rho<-0.50
size<-100
gammap<-0
alphap<-2
betap <-4
a<-gammap
b<-alphap
c<-alphap+betap
MEr1MonteCarloRho(size,rho,r,MEr1RhoSELFf,a,b,c,-1/r)
MEr1MonteCarloRho(size,rho,r,MEr1RhoPLFf,a,b,c,-1/r)

# testing Monte Carlo table for rho MLE state-of-the-art
# rem_var()
r<-1
rhos<-c(0.50)
sizes<-c(100)
MEr1MonteCarloRhoSoATab(sizes,rhos,r)


# testing Monte Carlo table for rho MLE
r<-1
rhos<-c(0.50)
sizes<-c(100)
MEr1MonteCarloRhoTab(sizes,rhos,r,MEr1RhoMLEf)

# testing Monte Carlo table for Lq MLE
r<-1
Lqs<-c(0.50)
sizes<-c(100)
MEr1MonteCarloLqTab(sizes,rhos,r,MEr1LqMLEf)

# testing Monte Carlo table for creadible region for Gauss hypergeometric
r<-1
rho<-0.50
size<-100
a<-1
b<-2
c<-6
gammap<-a
alphap<-b
betap <-c-b
MEr1MonteCarloCR(size,rho,r,gammap,alphap,alphap+betap)

################################################################################
# THE END
################################################################################
