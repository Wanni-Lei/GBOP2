##library(Rcpp)
##library(RcppArmadillo)
##library(Rcpp11)
#library(MASS)
##library(RcppEigen)
##library(progress)
# sourceCpp("Calculation.cpp")
# sourceCpp(file="Calculation.cpp",cacheDir="cache")



maxresp <- function(cutoff,dprior,ntr,phi,contrast)
{ ##dprior: Dirichlet prior
  ##contrast: design vector: b with elements of 0 and 1
  ##ntr=n=nobs.seq[i]=nobs[i]:  sample sizes
  ##ntr total number of patients treated at current,0:ntr means number of repsonses from 0 to ntr.
  ncut=NULL
  #if(ncol(contrast)==1) contrast = t(contrast)
  contrast = matrix(contrast,nrow=length(phi))
  for (i in 1:length(phi))
  ##0:ntr+sum(dprior[which(contrast[i,]==1)]): sum(b_k(a_k+x_k)) b_k=(1,0..,1), sum(x_k)=number of patients belongs to kth category
  ## number of reponses from 0 to ntr=n to find out the maximum reponses to stop the trial
  ## (whihc returns order)  -1=number of patients
  {ncut <- c(ncut, max(which(1-pbeta(phi[i],0:ntr+sum(dprior[which(contrast[i,]==1)]),ntr:0+sum(dprior)-sum(dprior[which(contrast[i,]==1)]))<cutoff))-1)}
  return (ncut)
}


##Binary Efficacy, Ordinal Efficacy, Multiple Efficacy
getboundary <- function(dprior,contrast,nobs,b,pow,phi)
{ ##nobs: interim data   n: the number of patients, see interims:
  stopbound = NULL
  nmax=max(nobs)
  for (n in nobs)
  { ## Pr(b*theta<=phi|D)>C = 1-b*(n/nmax)^pow
    stopbound = cbind(stopbound,maxresp(cutoff=b*(n/nmax)^pow,dprior=dprior,contrast=contrast,ntr=n,phi=phi))
  }
  boundary <- rbind(nobs,stopbound)
  boundary[which(boundary>999)]=999
  boundary[which(boundary< -999)]=-999
  return(boundary)
}

##bound of Binary Toxicity
getboundary.bitox <- function(dprior,contrast=as.matrix(1),nobs,b,pow,phi)
{
  stopbound = NULL
  nmax=max(nobs)
  for (n in nobs)
  {
    stopbound = cbind(stopbound,maxresp(cutoff=b*(n/nmax)^(pow/3),dprior=dprior,contrast=contrast,ntr=n,phi=phi))
  }
  boundary <- rbind(nobs,nobs-stopbound)
  boundary[which(boundary>999)]=999
  boundary[which(boundary< -999)]=-999
  ###patients who do not experience toxitity
  return(boundary)
}

##bound of Efficacy & Toxicity
#dimnames(bound.1) <- list(c("# patients treated","Stop if # CR <=","Stop if # CR/PR <="),NULL)
getboundary.tox <- function(dprior,contrast,nobs,nobsTX,b,pow,phi)
{
  stopbound1 = NULL
  stopbound2 = NULL
  nm = sort(union(nobs,nobsTX))
  nmax=max(nm)
  for (n in nm)
  {
    stopbound1 = cbind(stopbound1,maxresp(cutoff=b*(n/nmax)^pow,dprior=dprior,contrast=as.matrix(contrast[1,]),ntr=n,phi=phi[1]))
    if(!(n %in% nobs)) stopbound1[ncol(stopbound1)] = -Inf
    stopbound2 = cbind(stopbound2,maxresp(cutoff=b*(n/nmax)^(pow/3),dprior=dprior,contrast=as.matrix(contrast[2,]),ntr=n,phi=phi[2]))
    if(!(n %in% nobsTX)) stopbound2[ncol(stopbound2)] = -Inf

    # stopbound1 = cbind(stopbound1,maxresp(cutoff=b*(n/nmax)^pow,dprior=dprior,contrast=contrast,ntr=n,phi=phi))

  }
  # boundary <- rbind(nm,stopbound1[1,],nm-stopbound1[2,])

  boundary <- rbind(nm,stopbound1,nm-stopbound2)
  boundary[which(boundary>999)]=999
  boundary[which(boundary< -999)]=-999
  ##first row efficacy
  ## second row toxicity
  return(boundary)
}

getboundary.TTE <- function(gprior, phi, coeff, pow, n.interim, n.column=11)
{
  nmax = max(n.interim)
  boundary.full = c()
  #boundary.full = vector("list",length=length(n.interim))
  for (k in 1:length(n.interim))
  {
    boundary.temp = c()
    for (ss in 0:n.interim[k])
    {

      st = 0
      t.fu = 0
      while(st==0)
      {
        #if(pinvgamma(phi/log(2), shape = gprior[1]+ss, rate = gprior[2]+t.fu) <= 1 - coeff * (n.interim[k]/nmax)^pow) st=1

        if(pgamma(log(2)/phi, shape = gprior[1]+ss, rate = gprior[2]+t.fu,lower.tail = FALSE) <= 1 - coeff * (n.interim[k]/nmax)^pow) st=1
        else t.fu=t.fu+1
      }
      if(t.fu>0) {
        t.seq = seq(t.fu-1,t.fu,by=0.03)
        #t.bound = t.seq[max(which(pinvgamma(phi/log(2), shape = gprior[1]+ss, rate = gprior[2]+t.seq) > 1 - coeff * (n.interim[k]/nmax)^pow))]

        t.bound = t.seq[max(which(pgamma(log(2)/phi, shape = gprior[1]+ss, rate = gprior[2]+t.seq,lower.tail = FALSE) > 1 - coeff * (n.interim[k]/nmax)^pow))]
      }
      else t.bound=0

      boundary.temp = cbind(boundary.temp,c(n.interim[k],ss,t.bound))
    }
    boundary.full = cbind(boundary.full,boundary.temp)
  }
  return(boundary.full)
}


resp_list_r=function(seed, nsim,nobs,ptrue){

set.seed(seed)

resplist = list()
resplist_1=matrix(nrow = length(nobs), ncol = nsim)
resplist_2=matrix(nrow = length(nobs), ncol = nsim)
resplist_3=matrix(nrow = length(nobs), ncol = nsim)
resplist_4=matrix(nrow = length(nobs), ncol = nsim)
resplist[[1]] <- rmultinom(nsim,size=nobs[1],prob=ptrue)
resplist_1[1,] <- resplist[[1]][1,]
resplist_2[1,] <- resplist[[1]][2,]
resplist_3[1,] <- resplist[[1]][3,]
resplist_4[1,] <- resplist[[1]][4,]
for (j in 2:length(nobs))
{
  resplist[[j]] <- rmultinom(nsim,size=nobs[j]-nobs[j-1],prob=ptrue)
  resplist_1[j,]=resplist[[j]][1,]
  resplist_2[j,]=resplist[[j]][2,]
  resplist_3[j,]=resplist[[j]][3,]
  resplist_4[j,]=resplist[[j]][4,]
}
resplist_r = list()
resplist_r[[1]]=resplist_1
resplist_r[[2]]=resplist_2
resplist_r[[3]]=resplist_3
resplist_r[[4]]=resplist_4
return(resplist_r)
}

