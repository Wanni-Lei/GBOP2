

#library(data.table)
library(tidyr)
library(doParallel)
library(foreach)
library(parallel)
registerDoParallel(10)

#
den <-function(x,y,n,p){
   # p is a vector (p00, p01, p10, p11)
   #p11 is patients having response and toxicity
   #p01 = Pt - p11
   #p10 = Pr - p11
   Pr =p[3]+p[4] # total prob of response
   Pt = p[2] +p[4]
   Pt_1 = p[4]/Pr #condition toxicity from response patients p(0|1)
   Pt_0 = p[2]/(1-Pr) #condition toxicity from non response patients p(0|0)

   # This is done by iterating over all possible ways to distribute the y toxicities 
   # among the x response patients and the remaining n-x non-response patients.
   dx = dbinom(x,n,Pr)
   sumb=0
   for(j in 0:min(x,y)){
     # print(j)
      sumb = sumb+ dbinom(j, x,Pt_1)*dbinom((y-j), (n-x), Pt_0)
   }
   den = dx*sumb
   return(den)
}

#
# Calculate and fill the npts(new cohort) density matrix
#
cell_density <-function(interm, scenario){
   nintm = length(interm)
   den_list=list(nintm)
   new_patient = c(0, interm)
   for(n in seq(nintm)){
      npt = new_patient[n+1] - new_patient[n]
      #print(npt)
      den_list[[n]] =matrix(NA, (npt+1),(npt+1))
      for(i in 0:npt){
         for(j in 0:npt){
            den_list[[n]][i+1,j+1] =den(i,j,npt,scenario)
         }
      }
   }
   return(den_list)
}
#
#This function used to calculate the density of each cell of k-1 full matrix
#
D_den <-function(r,  # number of efficacy
                 t,  # number of toxicty
                 RD, # cumulative of events to make go decision
                 TD, # cumulative of toxicity events to make go decision
                 s,  # stage of interm
                 Dk_1, # last D matrix
                 cr,   # efficacy boundary [s-1 :s]
                 den_list,  # function to calculate the d
                 env = parent.frame()# boundary of last two for index
){
   sump=0
   rc= c(RD[RD$s==r,1:2]+1)
   if(cr >0 ){ rc[[1]] = rc[[1]] -cr-1}
   tc=c(TD[TD$s==t,1:2]+1)
   for(r in 1:length(rc[[1]])){
      for(t in 1:length(tc[[1]])){
         a1=rc[[1]][r] # previous round
         a2=rc[[2]][r] # this round
         b1=tc[[1]][t] # previous round
         b2=tc[[2]][t] # this round 
         sump=sump+Dk_1[a1,b1]*den_list[[s]][a2,b2]
      }
   }
   return(sump)
}


#
#This function used to fill in the density of each cell of k-1 full matrix
#
D_matrix <-function(s,r_bound,t_bound,npt,r_interm,t_interm,Dk_1,den_list){

  # s=2
  # r_bound=bond.eff.temp
  # t_bound = bond.tox.temp
  # r_interm = interm.temp
  # t_interm = interm.temp
  # Dk_1 = D_mat[[s-1]]
  # den_list

   rr = rev(expand.grid(0:npt[s], (r_bound[s-1]+1):r_interm[s-1])) # previous round number of events (make go decision), current round number of events
   rr$s=rowSums(rr)
   RR = rr[rr$s >r_bound[s],] # RR is the number of cumulative efficacy events by s condition on going from s-1; still make go decision


   tt = rev(expand.grid(0:npt[s],0:(t_bound[s-1]-1))) # previous round number of events (make go decision), current round number of events
   tt$s =rowSums(tt)
   TT = tt[tt$s<t_bound[s],] # TT is the number of cumulative toxicity events by s  condition on going from s-1; still make go decision


   D_mat = matrix(NA, nrow=(r_interm[s]-r_bound[s]), ncol=t_bound[s]) # probability matrix for making go decision at this round
   # k_com = as.matrix(CJ((r_bound[s]+1):r_interm[s],0:(t_bound[s]-1))) 
   k_com =as.matrix(expand.grid(0:(t_bound[s]-1),(r_bound[s]+1):r_interm[s])) # event counts for making go decision
   k_com = k_com[,c(2,1)]  # combination of efficacy and toxicity events (all possible X 2)
   print(nrow(k_com))
   cr=r_bound[(s-1):s]
   for(i in 1:nrow(k_com)){
      r=k_com[i,1] # efficacy events
      t=k_com[i,2] # toxicity events
      #D_mat[[s]][r-11,t+1] = den_k3(r=r, t=t)
      D_mat[r-cr[2],t+1] = D_den(r=r, t=t,RD=RR, TD=TT,s=s,Dk_1=Dk_1,cr=cr[1],den_list) # r-cr[2]: efficacy count for this round; # t+1: toxic count for this round
   }

   return(D_mat)

}

BOP2_TE<- function(interm.eff, interm.tox, boundary.eff, boundary.tox, den_list){


   interm.temp = unique(sort(c(interm.eff, interm.tox)))
   bond.eff.temp <- numeric(length(interm.temp))  # Create a new numeric vector a1 with the same length as a
   bond.tox.temp <- numeric(length(interm.temp))
   for (i in 1:length(interm.temp)) {
      if (interm.temp[i] %in% interm.eff) {
         bond.eff.temp[i] <- boundary.eff[interm.eff == interm.temp[i]]  # Fill a1 with the corresponding value from b
      } else {
         bond.eff.temp[i] <- -1  # Fill a1 with -1 if the element is not in b
      }
   }

   j =1
   for (i in 1:length(interm.temp)) {
      if (interm.temp[i] %in% interm.tox) {
         bond.tox.temp[i] <- boundary.tox[interm.tox == interm.temp[i]]  # Fill a1 with the corresponding value from b
      } else {
         bond.tox.temp[i] <- boundary.tox[j]
      }
      j = ifelse(j >= length(interm.tox), length(interm.tox), j+1)
   }
   nintm = length(interm.temp)
   new_patient = c(0, interm.temp)
   npt=rep(NA, nintm)
   for(j in seq(nintm)){
      npt[j] = new_patient[j+1] - new_patient[j]
   }
   #den_list=cell_density(interm.temp,scenario)
   D_mat=list(nintm) # save the probability of D() as a list for each interm from npt patients
   D_mat[[1]] = matrix(den_list[[1]][(bond.eff.temp[1]+2):(npt[1]+1),1:bond.tox.temp[1]],ncol=bond.tox.temp[1])
   
   for(s in 2:nintm){
      D_mat[[s]] = D_matrix(s,bond.eff.temp,bond.tox.temp,npt,interm.temp,interm.temp,D_mat[[s-1]],
                            den_list=den_list)
   }
   return(sum(D_mat[[nintm]]))
}





##Calculate, probability, average sample size and early stop rate
BOP2_TE.OC<- function(interm.eff, interm.tox, boundary.eff, boundary.tox, scenario){


   interm.temp = unique(sort(c(interm.eff, interm.tox)))
   bond.eff.temp <- numeric(length(interm.temp))  # Create a new numeric vector a1 with the same length as a
   bond.tox.temp <- numeric(length(interm.temp))
   for (i in 1:length(interm.temp)) {
      if (interm.temp[i] %in% interm.eff) {
         bond.eff.temp[i] <- boundary.eff[interm.eff == interm.temp[i]]  # Fill a1 with the corresponding value from b
      } else {
         bond.eff.temp[i] <- -1  # Fill a1 with -1 if the element is not in b
      }
   }

   j =1
   for (i in 1:length(interm.temp)) {
      if (interm.temp[i] %in% interm.tox) {
         bond.tox.temp[i] <- boundary.tox[interm.tox == interm.temp[i]]  # Fill a1 with the corresponding value from b
      } else {
         bond.tox.temp[i] <- boundary.tox[j]
      }
      j = ifelse(j >= length(interm.tox), length(interm.tox), j+1)
   }
   nintm = length(interm.temp)
   new_patient = c(0, interm.temp)
   npt=rep(NA, nintm)
   for(j in seq(nintm)){
      npt[j] = new_patient[j+1] - new_patient[j]
   }
   den_list=cell_density(interm.temp,scenario)
   D_mat=list(nintm) # save the probability of D() as a list for each interm from npt patients
   N=0
   D_mat[[1]] = matrix(den_list[[1]][(bond.eff.temp[1]+2):(npt[1]+1),1:bond.tox.temp[1]],ncol=bond.tox.temp[1])
   N= npt[1]
   
   for(s in 2:nintm){
     D_mat[[s]] = D_matrix(s,bond.eff.temp,bond.tox.temp,npt,interm.temp,interm.temp,D_mat[[s-1]],den_list)
     N = N + (npt[s]*sum(D_mat[[s-1]]))
   }
   
   

   ESS=round(N,2)

   PET=round((1-sum(D_mat[[nintm-1]])),4) # prob. of early stopping
   PCP=round(sum(D_mat[[nintm]]),4)


   res= matrix(NA, nrow=1, ncol=3)
   colnames(res)= c("PCP", 'PET', "ESS" )
   res[1,]=c(PCP=PCP, PET=PET,ESS=ESS)
   return(res)
}


get_boundarycpp <-function(interm.eff, interm.tox, lambda_e, lambda_t,
                           gamma,
                           prior= prior, r0=0.3, t0=0.4){
  
  contrast = rbind(c(0,0,1,1), c(0,1,0,1))
  interm = unique(sort(c(interm.eff, interm.tox)))
  N = max(interm)
  boundary = NULL
  for(nn in interm){
    p= matrix(NA, nrow=2, ncol=nn)
    for(n in seq(nn)){
      alpha = n +  contrast %*% prior
      beta = nn + sum(prior) - alpha
      # print(alpha)
      # print(beta)
      # if (n > nn){
      #   print(interm.eff)
      #   print(interm)
      #   warning(interm)
      #   warning(paste(lambda_e, lambda_t, gamma, 
      #                 n,nn,"n is out of the bound."))
      #   warning(prior)
      # }
      if (n > nn){
        n = nn
      }
      
      if(length(alpha) == 2 && length(beta) == 2){  # Check to prevent out-of-bounds error
        p[1,n] = 1 - pbeta(r0, alpha[1], beta[1])
        p[2,n] = pbeta(t0, alpha[2], beta[2])
      } else {
        warning("Alpha or beta has unexpected dimensions.")
      }
    }
    
    cr = lambda_e*(nn/N)^gamma
    ct = lambda_t*(nn/N)^(gamma/3)
    
    if(min(p[1,]) < cr){
      eff = max(which(p[1,]< cr))
    }
    else{eff = 0}
    tox =min(which(p[2,]<ct))
    b= rbind(nn,eff, tox)
    boundary = cbind(boundary, b)
  }
  
  bond.eff = boundary[c(1,2),][,boundary[1, ] %in% interm.eff][2,]
  bond.tox = boundary[c(1,3),][,boundary[1, ] %in% interm.tox][2,]
  
  boundary1 = list(
    boundary.eff = bond.eff,
    boundary.tox = bond.tox
  )
  return(boundary1)
  # return the list of boundary eff and tox, the length of eff and tox may be different
}

get_boundarycpp_2lambda <-function(interm.eff, interm.tox, lambda_e, lambda_t,
                           lambda_e2, lambda_t2, gamma,
                           prior= prior, r0=0.3, t0=0.4){

   contrast = rbind(c(0,0,1,1), c(0,1,0,1))
   interm = unique(sort(c(interm.eff, interm.tox)))
   N = max(interm)
   boundary = NULL
   for(nn in interm){
      p= matrix(NA, nrow=2, ncol=nn)
      for(n in seq(nn)){
         alpha = n +  contrast %*% prior
         beta = nn + sum(prior) - alpha
         # print(alpha)
         # print(beta)
         # if (n > nn){
         #   print(interm.eff)
         #   print(interm)
         #   warning(interm)
         #   warning(paste(lambda_e, lambda_t, gamma, 
         #                 n,nn,"n is out of the bound."))
         #   warning(prior)
         # }
         if (n > nn){
           n = nn # to prevent error from C++
         }
         
         if(length(alpha) == 2 && length(beta) == 2){  # Check to prevent out-of-bounds error
           p[1,n] = 1 - pbeta(r0, alpha[1], beta[1])
           p[2,n] = pbeta(t0, alpha[2], beta[2])
         } else {
           warning("Alpha or beta has unexpected dimensions.")
         }
      }
      
      if (nn == interm[length(interm)]){
        cr = lambda_e2
        ct = lambda_t2
      } else {
        cr = lambda_e*(nn/N)^gamma
        ct = lambda_t*(nn/N)^(gamma/3)
      }
    
      if(min(p[1,]) < cr){
         eff = max(which(p[1,]< cr))
      }
      else{eff = 0}
      
      tox =min(which(p[2,]<ct))
      b= rbind(nn,eff, tox)
      boundary = cbind(boundary, b)
   }

   bond.eff = boundary[c(1,2),][,boundary[1, ] %in% interm.eff][2,]
   bond.tox = boundary[c(1,3),][,boundary[1, ] %in% interm.tox][2,]

   boundary1 = list(
      boundary.eff = bond.eff,
      boundary.tox = bond.tox
   )
   return(boundary1)
   # return the list of boundary eff and tox, the length of eff and tox may be different
}

get_boundarycpp_efftox <-function(interm.eff, interm.tox, interm.efftox, 
                                  lambda_e, lambda_t, lambda_et, gamma_e, gamma_t, gamma_et,
                           prior= prior, criteria){
  
  contrast = rbind(c(0,0,1,0), c(0,1,0,0), c(0,0,0,1))
  interm = unique(sort(c(interm.eff, interm.tox, interm.efftox)))
  N = max(interm)
  boundary = NULL
  for(nn in interm){
    p= matrix(NA, nrow=3, ncol=nn)
    for(n in seq(nn)){
      alpha = n +  contrast %*% prior
      beta = nn + sum(prior) - alpha
      p[1,n] =1-pbeta(criteria[3], alpha[1], beta[1])
      p[2,n] =pbeta(criteria[2], alpha[2], beta[2])
      p[3,n] =pbeta(criteria[4], alpha[3], beta[3])
    }
    
    cr = lambda_e*(nn/N)^gamma_e
    ct = lambda_t*(nn/N)^(gamma_t)
    crt = lambda_et*(nn/N)^(gamma_et)
    
    if(min(p[1,]) < cr){
      eff = max(which(p[1,]< cr))
    }
    else{eff = 0}
    tox =min(which(p[2,]<ct))
    efftox =min(which(p[3,]<crt))
    
    b= rbind(nn,eff, tox, efftox)
    boundary = cbind(boundary, b)
  }
  
  bond.eff = boundary[c(1,2),][,boundary[1, ] %in% interm.eff][2,]
  bond.tox = boundary[c(1,3),][,boundary[1, ] %in% interm.tox][2,]
  bond.efftox = boundary[c(1,4),][,boundary[1, ] %in% interm.efftox][2,]
  
  boundary1 = list(
    boundary.eff = bond.eff,
    boundary.tox = bond.tox,
    boundary.efftox = bond.efftox
  )
  return(boundary1)
  # return the list of boundary eff and tox, the length of eff and tox may be different
}


# tt=get_boundarycpp(interm.eff=c(20,40), interm.tox=c(10,20,40), lambda_e=0.9, lambda_t=0.85, gamma=0.93,
#                            prior= c(0.42 ,0.28, 0.18, 0.12), r0=0.3, t0=0.4)
#
# tt2=get_boundarycpp(interm.eff=c(20,40), interm.tox=c(10,20,40), lambda_e=0.9, lambda_t=0.85, gamma=0.93,
#                    prior= c(0.3 ,0.1, 0.5, 0.1), r0=0.3, t0=0.4)


get_boundary_table <-function(interm.eff, interm.tox, lambda_e, lambda_t, gamma,
                              prior= prior, r0=0.3, t0=0.4){

   boundary=get_boundarycpp(interm.eff= interm.eff,interm.tox =interm.tox,
                            lambda_e=lambda_e, lambda_t=lambda_t, gamma=gamma,
                            prior=prior , r0=r0, t0=t0)
   boundary.eff = boundary[[1]]
   boundary.tox = boundary[[2]]
   interm.temp = unique(sort(c(interm.eff, interm.tox)))
   bond.eff.temp <- numeric(length(interm.temp))  # Create a new numeric vector a1 with the same length as a
   bond.tox.temp <- numeric(length(interm.temp))

   for (i in 1:length(interm.temp)) {
      if (interm.temp[i] %in% interm.eff) {
         bond.eff.temp[i] <- boundary.eff[interm.eff == interm.temp[i]]  # Fill a1 with the corresponding value from b
      } else {
         bond.eff.temp[i] <- NA  # Fill a1 with -1 if the element is not in b
      }
   }

   j =1
   for (i in 1:length(interm.temp)) {
      if (interm.temp[i] %in% interm.tox) {
         bond.tox.temp[i] <- boundary.tox[interm.tox == interm.temp[i]]  # Fill a1 with the corresponding value from b
      } else {
         bond.tox.temp[i] <- NA
      }
      j = ifelse(j >= length(interm.tox), length(interm.tox), j+1)
   }


   boundary = cbind(interm.temp,bond.eff.temp,bond.tox.temp )
   boundary[is.na(boundary)] = "NA"
   colnames(boundary) = c("# patients treated","Stop if # response <="," OR # toxicity >=")
   return(boundary)
}





grid_search_cpp <-function(interm.eff,interm.tox, scenario,lambda_r,lambda_t,
                           gamma_space,
                           r1, t1,
                           r0, t0,
                           typeI01,
                           typeI10,
                           typeI00,
                           lpow=0.2,
                           prior= scenario[1,],
                           true_opt=TRUE,
                           cppfunction
){
   res = matrix(NA, nrow=length(lambda_r)*length(lambda_t)*length(gamma_space)*2, ncol=9)
   colnames(res)= c("lambda_le","lambda_lt", "gamma", "a01", 'a10', 'a11','a00','pts_e', 'pts_t')
   i=s=0

   boundary_last = list(c(NA,1),c(NA,1))
   boundary_hist = matrix(NA, nrow=length(lambda_r)*length(lambda_t)*length(gamma_space),
                          ncol=(length(interm.eff)) + (length(interm.tox)))
   for(le in lambda_r){
      for(g in gamma_space){
         for(lt in lambda_t){
            i = i+1
            boundary=get_boundarycpp(interm.eff= interm.eff,interm.tox =interm.tox,
                                     lambda_e=le, lambda_t=lt, gamma=g,
                                     prior=prior , r0=r0, t0=t0)


            if(s==0 && identical(boundary, boundary_last)){
               res[i,] = c(le,lt,g, a01, a10, a11,a00,temp_pe$pts, temp_pt$pts)
               next
            }
            s=0
            boundary_last = boundary

            #scenario = c(PE.a, PE.n PT.n, PT.a)
            temp_pe = Exacterror(nobs=interm.eff,ncut=boundary[[1]],
                                             pnull=r0,palter=r1)

            temp_pt = Exacterror(nobs=interm.tox,ncut=interm.tox-boundary[[2]],
                                             pnull=1-t0,palter=1-t1)

            a00 = temp_pe$t1err * temp_pt$t1err
            a01 = temp_pe$t1err * temp_pt$power
            a10 = temp_pe$power * temp_pt$t1err
            a11 = temp_pe$power * temp_pt$power

            if(a11 <lpow){
               next
            }

            if(a10 >(1.5*typeI10) || a01 > (1.5*typeI01) || a00 > (1.5*typeI00)){
               s =1
               break
            }
            res[i,] = c(le,lt,g, a01, a10, a11, a00, temp_pe$pts, temp_pt$pts)
         }
      }
   }
   # Optimization with highest power
   #best = li[order(li[,5],decreasing=TRUE),]
   res =  na.omit(res)


   res1 = subset(res, res[,4] <typeI01 & res[,5] < typeI10)
   best = round(res1[order(res1[,6],decreasing=TRUE), ],4)

   ## saved for the boundary that are slight over the type I error
   best_all = round(res[order(res[,6],decreasing=TRUE), ],4)


   best_all[,4] = abs(best_all[,4] - typeI01)
   best_all[,5] = abs(best_all[,5] - typeI10)
   best_all[,7] = abs(best_all[,7] - typeI10)

   m1=best_all[best_all[,4] == min(best_all[,4]),]
   m1 = rbind(m1,m1) # to adoive issue when only one row return from above
   m1=m1[m1[,5] == min(m1[,5]),]

   m2=best_all[best_all[,5] == min(best_all[,5]),]
   m2= rbind(m2,m2)  # to adoive issue when only one row return from above
   m2=m2[m2[,4] == min(m2[,4]),]


   m3 = best_all[best_all[,7] == min(best_all[,7]),]




   if((m2[1,4]+m2[1,5]) >(m1[1,4] + m1[1,5])){
      optim_backup = m1[1,]
   }else{
      optim_backup = m2[1,]
   }

   if(true_opt==FALSE & optim_backup[6] > best[1,6] ){
      optimal = optim_backup
   }else{
      optimal = best[1,]
   }

   return(optimal)

}




grid_search_cpp_minN <-function(interm.eff,interm.tox, scenario,lambda_r,lambda_t,
                           gamma_space,
                           r1, t1,
                           r0, t0,
                           typeI01,
                           typeI10,
                           typeI00,
                           lpow=0.2,
                           power_e,
                           power_t,
                           prior= scenario[1,],
                           true_opt=TRUE,
                           cppfunction
){
  res = matrix(NA, nrow=length(lambda_r)*length(lambda_t)*length(gamma_space), ncol=12)
  colnames(res)= c("lambda_le","lambda_lt", "gamma", "a01", 'a10', 'a11','a00','power_e','power_t','pts_e', 'pts_t', 'pts')
  i=s=0
  
  boundary_last = list(c(NA,1),c(NA,1))
  boundary_hist = matrix(NA, nrow=length(lambda_r)*length(lambda_t)*length(gamma_space),
                         ncol=(length(interm.eff)) + (length(interm.tox)))
  for(le in lambda_r){
    for(g in gamma_space){
      for(lt in lambda_t){
        i = i+1
        boundary=get_boundarycpp(interm.eff= interm.eff,interm.tox =interm.tox,
                                 lambda_e=le, lambda_t=lt, gamma=g,
                                 prior=prior , r0=r0, t0=t0)
        
        
        if(s==0 && identical(boundary, boundary_last)){
          res[i,] = c(le,lt,g, a01, a10, a11,a00,temp_pe$power, temp_pt$power, temp_pe$pts, temp_pt$pts, temp_pe$pts+temp_pt$pts)
          next
        }
        s=0
        boundary_last = boundary
        
        #scenario = c(PE.a, PE.n PT.n, PT.a)
        temp_pe = Exacterror(nobs=interm.eff,ncut=boundary[[1]],
                             pnull=r0,palter=r1)
        
        temp_pt = Exacterror(nobs=interm.tox,ncut=interm.tox-boundary[[2]],
                             pnull=1-t0,palter=1-t1)
        
        a00 = temp_pe$t1err * temp_pt$t1err
        a01 = temp_pe$t1err * temp_pt$power
        a10 = temp_pe$power * temp_pt$t1err
        a11 = temp_pe$power * temp_pt$power
        
        if(a11 <lpow){
          next
        }
        
        if(a10 >(1.5*typeI10) || a01 > (1.5*typeI01) || a00 > (1.5*typeI00)){
          s =1
          break
        }
        res[i,] = c(le,lt,g, a01, a10, a11, a00, temp_pe$power, temp_pt$power, temp_pe$pts, temp_pt$pts, temp_pe$pts+temp_pt$pts)
      }
    }
  }
  # Optimization with highest power
  #best = li[order(li[,5],decreasing=TRUE),]
  res =  na.omit(res)
  
  res1 = subset(res, res[,4] <= typeI01 & res[,5] <= typeI10 & res[,8] >= power_e & res[,9] >= power_t)
  
  best = round(res1[order(res1[,12],decreasing=FALSE), ],4)
  optimal = best[1,]
  
  return(optimal)
  
}




# scen = hypotheses_ind(0.3,0.6,0.4,0.2)
#
# cut.seq = rev(c(seq(0.5,0.8,0.025),seq(0.81,0.96,0.01)))
# power.seq = log(seq(1,0.5,by=-0.025))/log(0.5)
#
# out1=grid_search_cpp(interm.eff=c(10,30,50)
#                            ,interm.tox=c(10,30,50),
#                            scenario=scen,
#                            lambda_r=cut.seq,
#                            lambda_t=cut.seq,
#                            gamma_space=power.seq,
#                            r1=0.6, t1=0.2,
#                            r0=0.3, t0=0.4,
#                            typeI01=0.10,
#                            typeI10=0.10,
#                            typeI00=0.025,
#                            lpow=0.2,
#                            true_opt=TRUE,
#                            cppfunction=Calc
# )
#




###### Based on boundary searching ####
# find out the unique boundary from  parameter space


bond.unique <-function(interm.eff,interm.tox, scenario,lambda_r,lambda_t,
                       gamma_space, r0, t0, prior = scenario[1,]

){
   boundary_all = matrix(NA, nrow=length(lambda_r)*length(lambda_t)*length(gamma_space),
                          ncol=(length(interm.eff)) + (length(interm.tox)))
   para_all = matrix(NA, nrow=length(lambda_r)*length(lambda_t)*length(gamma_space),
                     ncol=3)
   i=0
   for(le in lambda_r){
      for(g in gamma_space){
         for(lt in lambda_t){
            i = i+1
            boundary=get_boundarycpp(interm.eff= interm.eff,interm.tox =interm.tox,
                                     lambda_e=le, lambda_t=lt, gamma=g,
                                     prior=prior , r0=r0, t0=t0)
            boundary_all[i,] =c(boundary[[1]],boundary[[2]])
            para_all[i,] = c(le, lt, g)

         }
      }
   }
   return(list(boundary_search = unique(boundary_all),
               parameter_search = para_all,
               boundary_all =boundary_all)
          )
}

bond.unique_2lambda <-function(interm.eff,interm.tox, scenario,lambda_r,lambda_t,lambda_r2,lambda_t2,
                       gamma_space, r0, t0, prior = scenario[1,]
                       
){
  boundary_all = matrix(NA, nrow=length(lambda_r)*length(lambda_t)*length(lambda_r2)*length(lambda_t2)*length(gamma_space),
                        ncol=(length(interm.eff)) + (length(interm.tox)))
  para_all = matrix(NA, nrow=length(lambda_r)*length(lambda_t)*length(lambda_r2)*length(lambda_t2)*length(gamma_space),
                    ncol=5)
  i=0
  for(le in lambda_r){
    for(le2 in lambda_r2){
      for(g in gamma_space){
        for(lt in lambda_t){
          for(lt2 in lambda_t2){
            i = i+1
            boundary=get_boundarycpp_2lambda(interm.eff= interm.eff,interm.tox =interm.tox,
                                     lambda_e=le, lambda_t=lt, lambda_e2=le2, lambda_t2=lt2, gamma=g,
                                     prior=prior , r0=r0, t0=t0)
            boundary_all[i,] =c(boundary[[1]],boundary[[2]])
            para_all[i,] = c(le, lt, le2, lt2, g)
          }
        }
      }
    }
    
  }
  return(list(boundary_search = unique(boundary_all),
              parameter_search = para_all,
              boundary_all =boundary_all)
  )
}

grid_search_corr<-function(boundary_search,interm.eff,interm.tox, scenario,
                           typeI01, typeI10,typeI00, lpow){


   ll = dim(boundary_search)[1]
   eff= boundary_search[ , 1:length(interm.eff)]
   tox = boundary_search[ , (length(interm.eff)+1):(length(interm.eff)+length(interm.tox))]


   pow = rep(0, ll)
   a10_e = rep(0, ll)
   a01_t = rep(0, ll)
   a00_et = rep(0, ll)
   
   N00_et = rep(0, ll)
   N01_t = rep(0, ll)
   N10_e = rep(0, ll)
   N11_power = rep(0, ll)

   for(i in 1:ll){
      
      temp = BOP2_TE.OC(interm.eff, interm.tox, boundary.eff=eff[i,], boundary.tox=tox[i,], scenario=scenario[2,])
      N01 = temp[3]
      a01= temp[1]

      if(a01 > typeI01){
         next
      }

      temp= BOP2_TE.OC(interm.eff, interm.tox, boundary.eff=eff[i,], boundary.tox=tox[i,], scenario=scenario[3,])
      N10 = temp[3]
      a10= temp[1]
      if(a10 >typeI10){
         next
      }

      temp = BOP2_TE.OC(interm.eff, interm.tox, boundary.eff=eff[i,], boundary.tox=tox[i,], scenario=scenario[1,])
      N00 = temp[3]
      a00= temp[1]
      if(a00 >typeI00){
         next
      }

      temp = BOP2_TE.OC(interm.eff, interm.tox, boundary.eff=eff[i,], boundary.tox=tox[i,], scenario=scenario[4,])
      N11 = temp[3]
      a11= temp[1]
      pow[i] = a11
      a10_e[i] = a10
      a01_t[i] = a01
      a00_et[i] = a00
      N00_et[i] = N00
      N10_e[i] = N10
      N01_t[i] = N01
      N11_power[i] = N11
      
   }
   if(max(pow) >0){
      bond2 = boundary_search[pow== max(pow),]
      a10 = a10_e[pow== max(pow)]
      a01 = a01_t[pow==max(pow)]
      a00 = a00_et[pow==max(pow)]
      return(c(a01, a10, max(pow),a00,bond2))
   }
}


grid_search_corr_minN<-function(boundary_search,interm.eff,interm.tox, scenario,
                           typeI01, typeI10,typeI00, power11,lpow){
  
  
  ll = dim(boundary_search)[1]
  eff= boundary_search[ , 1:length(interm.eff)]
  tox = boundary_search[ , (length(interm.eff)+1):(length(interm.eff)+length(interm.tox))]
  
  
  pow = rep(0, ll)
  a10_e = rep(0, ll)
  a01_t = rep(0, ll)
  a00_et = rep(0, ll)
  
  N00_et = rep(0, ll)
  N01_t = rep(0, ll)
  N10_e = rep(0, ll)
  N11_power = rep(0, ll)
  
  for(i in 1:ll){
    
    temp = BOP2_TE.OC(interm.eff, interm.tox, boundary.eff=eff[i,], boundary.tox=tox[i,], scenario=scenario[2,])
    N01 = temp[3]
    a01= temp[1]
    
    if(a01 > typeI01){
      next
    }
    
    temp= BOP2_TE.OC(interm.eff, interm.tox, boundary.eff=eff[i,], boundary.tox=tox[i,], scenario=scenario[3,])
    N10 = temp[3]
    a10= temp[1]
    if(a10 >typeI10){
      next
    }
    
    temp = BOP2_TE.OC(interm.eff, interm.tox, boundary.eff=eff[i,], boundary.tox=tox[i,], scenario=scenario[1,])
    N00 = temp[3]
    a00= temp[1]
    if(a00 >typeI00){
      next
    }
    
    temp = BOP2_TE.OC(interm.eff, interm.tox, boundary.eff=eff[i,], boundary.tox=tox[i,], scenario=scenario[4,])
    N11 = temp[3]
    a11= temp[1]
    if(a11 < power11){
      next
    }
    
    pow[i] = a11
    a10_e[i] = a10
    a01_t[i] = a01
    a00_et[i] = a00
    N00_et[i] = N00
    N10_e[i] = N10
    N01_t[i] = N01
    N11_power[i] = N11
    
    
  }
  
  ind = (N11_power != 0)
  # print(boundary_search[ind,])
  res = data.frame(matrix(boundary_search[ind,], ncol=4))
  res$a00 = a00_et[ind]
  res$a10 = a10_e[ind]
  res$a01 = a01_t[ind]
  res$power = pow[ind]
  res$N00 = N00_et[ind]
  res$N10 = N10_e[ind]
  res$N01 = N01_t[ind]
  res$N11 = N11_power[ind]
  
  return(res)
  
}



optimization_corr <- function(interm.eff, interm.tox, scenario, search.temp,
                       r0=0.3, t0=0.4, typeI01=0.10, typeI10=0.10, typeI00=0.025, lpow=0.5,
                       func=c("BOP2_TE", "grid_search_corr", "get_boundarycpp", "bond.unique",
                              "D_matrix","BOP2_TE.OC",
                              "D_den","cell_density","den"),
                       ncore=10,
                       Shiny=TRUE){

   t1 = Sys.time()


   # # search.temp = bond.unique(interm.eff=interm.eff, interm.tox=interm.tox, scenario=scenario,
   # #                        lambda_r=lambda_E_space,
   # #                        lambda_t=lambda_T_space,
   # #                        gamma_space=gamma_space,
   # #                        r0=r0, t0 =t0
   # )

   boundary_search = search.temp[[1]]
   parameter_search = search.temp[[2]]
   boundary_all = search.temp[[3]]
   print(dim(boundary_search)[1])
   nn = floor(dim(boundary_search)[1]/ncore)
   print(nn)


   opt <-foreach(
      nc=1:ncore,
      .export = func,
      .combine = 'rbind'
   ) %dopar%{
      a= nn*(nc-1)+1
      b=nn*nc

      if(nc ==ncore){
         b=dim(boundary_search)[1]
      }
      re1 = grid_search_corr(boundary_search=boundary_search[a:b,], interm.eff=interm.eff, interm.tox=interm.tox,
                             scenario= scenario,
                             typeI01= typeI01, typeI10 =typeI10,
                             typeI00 = typeI00, lpow=lpow)

      return(re1)
   }

    

   opt2 = opt[order(opt[,3], decreasing = T),]
   boundary_temp = opt2[1, 5:(dim(opt2)[2])]
   print(boundary_temp)


   row_matches <- apply(boundary_all, 1, function(row) all(row == boundary_temp))
   parameter_temp = parameter_search[row_matches,]

   print(parameter_temp[1:5,])

   t2 = Sys.time()
   print(t2-t1)

   if(Shiny==TRUE){

      return(c(lambda_E = parameter_temp[1,1], lambda_T = parameter_temp[1,2],
               gamma= parameter_temp[1,3],a01 = opt2[1,1], a10 =opt2[1,2], power=opt2[1,3], a00=opt2[1,4]))
   }
   else{
      rest = list(parameter =c(lambda_E = parameter_temp[1,1], lambda_T = parameter_temp[1,2],
                               gamma= parameter_temp[1,3]),
                  prob =c(a01 = opt2[1,1], a10 =opt2[1,2], power=opt2[1,3],a00=opt2[1,4]),
                  boundary =boundary_temp
      )
      return(rest)

   }

}



optimization_corr_minN <- function(interm.eff, interm.tox, scenario, search.temp,
                              r0=0.3, t0=0.4, typeI01=0.10, typeI10=0.10, typeI00=0.025, power11=0.8,lpow=0.5,
                              func=c("BOP2_TE", "grid_search_corr_minN", "get_boundarycpp", "bond.unique",
                                     "D_matrix","BOP2_TE.OC",
                                     "D_den","cell_density","den"),
                              ncore=10,
                              Shiny=TRUE){
  
  t1 = Sys.time()
  
  
  # # search.temp = bond.unique(interm.eff=interm.eff, interm.tox=interm.tox, scenario=scenario,
  # #                        lambda_r=lambda_E_space,
  # #                        lambda_t=lambda_T_space,
  # #                        gamma_space=gamma_space,
  # #                        r0=r0, t0 =t0
  # )
  
  boundary_search = search.temp[[1]]
  parameter_search = search.temp[[2]]
  boundary_all = search.temp[[3]]
  print(dim(boundary_search)[1])
  nn = floor(dim(boundary_search)[1]/ncore)
  print(nn)
  
  
  opt <-foreach(
    nc=1:ncore,
    .export = func,
    .combine = 'rbind'
  ) %dopar%{
    a= nn*(nc-1)+1
    b=nn*nc
    
    if(nc ==ncore){
      b=dim(boundary_search)[1]
    }
    re1 = grid_search_corr_minN(boundary_search=boundary_search[a:b,], interm.eff=interm.eff, interm.tox=interm.tox,
                           scenario= scenario,
                           typeI01= typeI01, typeI10 =typeI10,
                           typeI00 = typeI00, power11=power11,lpow=lpow)
    
    return(re1)
  }
  
  length_boundary = ncol(boundary_search)
  length_t1e_power = 4
  
  opt$Nfinal = (rowSums(opt[,c("N00", "N10", "N01")])/3 + opt[,"N11"])/2
  # opt$Nfinal = rowSums(opt[,(length_boundary+length_t1e_power+1):ncol(opt)])
  
  opt2 = opt[order(opt$Nfinal, -opt$power),] # order Nfinal in ascending and power in descending
  print(opt2)
  boundary_temp = opt2[1,1:length_boundary]
  print(boundary_temp)
  
  
  row_matches <- apply(boundary_all, 1, function(row) all(row == boundary_temp))
  parameter_temp = parameter_search[row_matches,]
  
  print(parameter_temp[1:5,])
  
  t2 = Sys.time()
  print(t2-t1)
  
  rest = list(parameter =c(lambda_E = parameter_temp[1,1], lambda_T = parameter_temp[1,2],
                           gamma= parameter_temp[1,3]),
              prob =c(a01 = opt2[1,"a01"], a10 =opt2[1,"a10"], a00=opt2[1,"a00"],power=opt2[1,"power"]),
              boundary =boundary_temp,
              SampleSize = opt2[1,(ncol(opt)-4):ncol(opt)]
  )
  return(rest)
  
}



scenario <-function(Pr, Pt, phi=1){
   ## Calculate Prt
   a = 1 + (Pr+Pt)*(phi -1)
   b = -4*phi*(phi-1)*Pr*Pt

   if(phi==1)
   {Prt = Pr*Pt}
   else{
      Prt = (a - sqrt(a^2 + b))/(2*(phi-1))
   }

   p11 = Prt
   p01 = Pt- Prt
   p10 = Pr - Prt
   p00 = 1 - Pr - Pt + Prt
   pi = c(p00, p01, p10, p11)
   return(list(P = pi, oddsratio = phi))
}



hypotheses_ind <- function(Ne, Ae, Nt, At){
   h00 = scenario(Ne, Nt, phi=1)$P
   h01 = scenario(Ne, At, phi=1)$P
   h10 = scenario(Ae, Nt, phi=1)$P
   h11 = scenario(Ae, At, phi=1)$P
   h4 = rbind(h00, h01, h10, h11)
   return(h4)
}

hypotheses_corr <-function(Ne, Ae, Nt, At, PA_ET, PN_ET){

   # joint probabilities of four possible outcomes 
   # (no efficacy and no toxicity, no efficacy and toxicity, efficacy and no toxicity, both efficacy and toxicity)
   h00 = c(1-Ne-Nt+PN_ET, Nt-PN_ET, Ne-PN_ET, PN_ET)

   phi = (PN_ET*(1-Ne-Nt+PN_ET))/((Ne-PN_ET)*(Nt-PN_ET))
   h01 = scenario(Ne, At, phi=phi)$P
   h10 = scenario(Ae, Nt, phi=phi)$P
   ### All the null hypotheses share the same correlation phi

   h11 = c(1-At-Ae+PA_ET, At-PA_ET, Ae-PA_ET, PA_ET)
   h4 = rbind(h00, h01, h10, h11)
   return(h4)
}















#
# parameter.tune <-function(interm.eff, interm.tox,
#                           r0, t0,
#                           boundary,
#                           prior,
#                           lambda_E_space,
#                           lambda_T_space,
#                           gamma_space){
#
#    all.space = as.matrix(expand.grid(lambda_E_space,lambda_T_space,
#                                      gamma_space))
#    param = NULL
#
#    for(i in 1:dim(all.space)[1]){
#
#       bond.temp = get_boundarycpp(interm.eff= interm.eff,interm.tox =interm.tox,
#                                   r0=r0, t0=t0,
#                                   lambda_e = all.space[i,1],
#                                   lambda_t= all.space[i,2],
#                                   gamma= all.space[i,3],
#                                   prior=prior
#       )
#
#       if(identical(bond.temp[[1]], boundary[[1]]) &
#          identical(bond.temp[[2]], boundary[[2]])){
#          param = rbind(param, all.space[i,])
#
#       }
#
#    }
#
#    return(param)
#
#
# }

# optim_bond_add <- function(boundary.add, interm.eff, interm.tox, scenario,lambda_E_space,lambda_T_space, gamma_space,
#                        r0=0.2, t0=0.4, typeI01=0.10, typeI10=0.10, lpow=0.5,
#                        func=c("BOP2_TE", "grid_search_bond", "get_boundarycpp", "bond.unique",
#                               "parameter.tune","D_matrix",
#                               "D_den","cell_density","den"),
#                        ncore=10){
#
#    t1 = Sys.time()
#
#    all.temp = boundary.add
#
#    print(dim(all.temp)[1])
#    nn = floor(dim(all.temp)[1]/ncore)
#    print(nn)
#
#    # for( nc in 1:10){
#    #    a= nn*(nc-1)+1
#    #    b=nn*nc
#    #    if(b >dim(all.temp)[1]){b=dim(all.temp)[1]}
#    #    print(c(a, b))
#    #
#    # }
#
#
#    opt <-foreach(
#       nc=1:ncore,
#       .export = func,
#       .combine = 'rbind'
#    ) %dopar%{
#       a= nn*(nc-1)+1
#       b=nn*nc
#
#       if(nc ==ncore){
#          b=dim(all.temp)[1]
#       }
#       re1 = grid_search_bond(boundary_all=all.temp[a:b,], interm.eff=interm.eff, interm.tox=interm.tox,
#                              scenario= scenario,
#                              typeI01= typeI01, typeI10 =typeI10, lpow=0.1)
#
#       return(re1)
#    }
#
#    opt2 = opt[order(opt[,3], decreasing = T),]
#
#
#    boundary_temp = list(
#       opt2[1, 4:(3+length(interm.eff))],
#       opt2[1,(length(interm.eff)+4): (dim(opt2)[2])]
#    )
#    print(boundary_temp)
#    para = parameter.tune(interm.eff=interm.eff, interm.tox=interm.tox,
#                          r0=r0, t0=t0,
#                          boundary = boundary_temp,
#                          prior = scenario[1,],
#                          lambda_E_space = lambda_E_space,
#                          lambda_T_space = lambda_T_space,
#                          gamma_space = gamma_space
#    )
#    #opt1 =as.data.frame(opt)
#    t2 = Sys.time()
#    print(t2-t1)
#
#    rest = list(parameter =c(lambda_E = para[1,1], lambda_T = para[1,2], gamma= para[1,3]),
#                Prob =c(a01 = opt2[1,1], a10 =opt2[1,2], power=opt2[1,3]),
#                Boundary =boundary_temp
#    )
#    return(rest)
# }
#

# grid_search_general<-function(interm.eff,interm.tox, scenario,
#                               lambda_r,lambda_t, gamma_space,
#                               r0, t0, typeI01, typeI10, lpow){
#
#    res = matrix(NA, nrow=length(lambda_r)*length(lambda_t)*length(gamma_space), ncol=6)
#    colnames(res)= c("lambda_le","lambda_lt", "gamma", "a01", 'a10', 'a11')
#    i=s=0
#
#    boundary_last = list(c(NA,1),c(NA,1))
#
#    for(g in gamma_space){
#       for(le in lambda_r ){
#          # Find the le first
#          lt= 0.98
#          boundary0=get_boundarycpp(interm.eff= interm.eff,interm.tox =interm.tox,
#                                    lambda_e=le, lambda_t=lt, gamma=g,
#                                    prior=scenario[1,] , r0=r0, t0=t0)
#
#          a10 = BOP2_TE(interm.eff= interm.eff, interm.tox = interm.tox,
#                        boundary.eff=boundary0[[1]], boundary.tox = boundary0[[2]],
#                        scenario=scenario[3,])
#          if(a10 >typeI10){
#             break
#          }
#
#          else{  # then find the lt
#             for(lt in lambda_t){
#                i = i+1
#                boundary=get_boundarycpp(interm.eff= interm.eff,interm.tox =interm.tox,
#                                         lambda_e=le, lambda_t=lt, gamma=g,
#                                         prior=scenario[1,] , r0=r0, t0=t0)
#
#                if(s==0 && identical(boundary, boundary_last)){
#                   res[i,] = c(le,lt,g, a01, a10, a11)
#                   next
#                }
#
#                s=0
#                boundary_last = boundary
#
#
#                a11 = BOP2_TE(interm.eff= interm.eff, interm.tox = interm.tox,
#                              boundary.eff=boundary[[1]], boundary.tox = boundary[[2]],
#                              scenario=scenario[4,])
#                if(a11 <lpow){
#                   s=1
#                   next
#                }
#
#                a01 = BOP2_TE(interm.eff= interm.eff, interm.tox = interm.tox,
#                              boundary.eff=boundary[[1]], boundary.tox = boundary[[2]],
#                              scenario=scenario[2,])
#                if(a01 >typeI01){
#                   s=1
#                   break
#                }
#
#                a10 = BOP2_TE(interm.eff= interm.eff, interm.tox = interm.tox,
#                              boundary.eff=boundary[[1]], boundary.tox = boundary[[2]],
#                              scenario=scenario[3,])
#                if(a10 >typeI10){
#                   s=1
#                   break
#                }
#
#                res[i,] = c(le,lt,g, a01, a10, a11)
#                print(c(le, lt, g))
#
#             }
#          }
#       }
#
#    }
#    res =  na.omit(res)
#    best = round(res [order(res[,6],decreasing=TRUE), ],4)
#    return(best)
# }
#
#
#
# optim_general <- function(interm.eff, interm.tox, scenario,lambda_E_space,lambda_T_space, gamma_space,
#                           r0=0.2, t0=0.4, typeI01=0.10, typeI10=0.10, lpow=0.5,
#                           func=c("BOP2_TE", "grid_search_general", "get_boundarycpp", "D_matrix",
#                                  "D_den","cell_density","den"),
#                           ncore=10){
#
#    t1 = Sys.time()
#    nn = ceiling(length(lambda_T_space)/ncore)
#
#    opt <-foreach(
#       nc=1:ncore,
#       .export = func,
#       .combine = 'rbind'
#    ) %dopar%{
#       a= nn*(nc-1)+1
#       b=nn*nc
#       if(b >length(lambda_T_space)){b=length(lambda_T_space)}
#       re1= grid_search_general(interm.eff=interm.eff,
#                                interm.tox = interm.tox,
#                                scenario=scenario,
#                                lambda_r=lambda_E_space,
#                                lambda_t=lambda_T_space[a:b],
#                                gamma_space=gamma_space,
#                                r0=r0, t0=t0, typeI01=typeI01,
#                                typeI10=typeI10, lpow=lpow)
#       return(re1)
#    }
#    t2 = Sys.time()
#    print(t2-t1)
#    opt2 = opt[order(opt[,6], decreasing = T),]
#    #opt1 =as.data.frame(opt)
#    return(opt2)
# }





















# out <-function(sc, bt,interm, eff_bond_less=TRUE){
#
#   oc= matrix(NA, nrow=3, ncol=4)
#   colnames(oc)= c('a00', 'a01', 'a10','a11')
#   rownames(oc) = c("Go(%)", 'PET(%)', "ASS" )
#   param = rep(NA,3)
#
#   res = bt[bt[,"a11"] == bt[1,"a11"], , drop=F]
#   if(dim(res)[2]==5){
#     res=res[order(res[,"lambda"]), ,drop=F]
#     param[1:2] = res[1,1]
#     param[3] =res[1,2]}
#   else{
#     res=res[order(res[,"lambda_le"]), ]
#     param[1:2] = res[1,1:2]
#     param[3] =res[1,3]
#   }
#   bond =get_boundary(interm=interm, r0=(sc[1,3]+sc[1,4]), t0=(sc[1,2]+sc[1,4]),lambda_e = param[1],
#
#                      lambda_t=param[2], gamma=param[3], prior=sc[1,])
#
#   for(i in 1:4){
#     oc[,i] = t(ASS_PET(interm,bond, sc[i,]))
#   }
#   bond= rbind(r=bond[[1]],t=bond[[2]])
#
#   out=list(boundary=bond, OC=oc, parameter=param)
#   return(out)
# }




