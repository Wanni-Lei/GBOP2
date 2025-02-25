

# ------------------------
DoStoppingBoundaries = function( input, globalPara, Calc, nobs.seq, nobs.seqTX, isSBcalculated, session )
{
  #browser();
  if(input$e1n>0.0001 & input$e1a>0.0001 & input$e1n<1 & input$e1a<1 &
     input$e2n>0.0001 & input$e2a>0.0001 & input$e2n<1 & input$e2a<1 &
     input$e3n>0.0001 & input$e3a>0.0001 & input$e3n<1 & input$e3a<1 &
     !is.na(input$e1n) & !is.na(input$e1a) & !is.na(input$e2n) & !is.na(input$e2a) & !is.na(input$e3n) & !is.na(input$e3a)){

    if(input$e1n<=input$e1a & input$e2n>=input$e2a) {

      if(input$e1n+input$e2n-input$e3n<1 & input$e1a+input$e2a-input$e3a<1 &
         input$e3n<min(input$e1n,input$e2n) & input$e3a<min(input$e1a,input$e2a)){

        if((input$ETprior1>=0 & input$ETprior1<=1 & input$ETprior2>=0 & input$ETprior2<=1 & input$ETprior3>=0 & input$ETprior3<=1 &
            !is.na(input$ETprior1) & !is.na(input$ETprior2) & !is.na(input$ETprior3) & !is.na(input$ETpriorSS) & input$ETpriorSS>0)| input$priorspec){

          if(input$ETprior1+input$ETprior2-input$ETprior3 <= 1  &  input$ETprior3 <= min(input$ETprior1,input$ETprior2)){

            # cut.seq = sort(c(seq(0.5,0.95,by=0.025),seq(0.96,0.99,0.01),seq(0.9-input$err1,1-input$err1,0.005)))
            cut.seq = rev(c(seq(0.5,0.8,0.025),seq(0.81,0.96,0.01)))
            power.seq = log(seq(1,0.5,by=-0.025))/log(0.5)
            p.n = c(input$e3n,input$e1n-input$e3n,input$e2n-input$e3n,1-input$e1n-input$e2n+input$e3n)
            ##input$e3n: Pr(Eff & Tox), input$e1n: Pr(Eff), input$e2n: Pr(Tox)
            ##input$e1n-input$e3n: Pr(Eff & no Tox)
            ##input$e2n-input$e3n: Pr(no Eff & Tox)
            ##1-input$e1n-input$e2n+input$e3n: Pr(no Eff & no Tox)
            p.a = c(input$e3a,input$e1a-input$e3a,input$e2a-input$e3a,1-input$e1a-input$e2a+input$e3a)
            # if(input$priorspec){
            #    prior = p.n
            # }else{
            #    prior = c(input$ETprior3*input$ETpriorSS,(input$ETprior1-input$ETprior3)*input$ETpriorSS,
            #              (input$ETprior2-input$ETprior3)*input$ETpriorSS,(1-input$ETprior1-input$ETprior2+input$ETprior3)*input$ETpriorSS)
            # }

            oc.mat.efftox = c()
            cut.start = 1

            # oc.mat.efftox = Calc$GridSearchToxEffRcpp(seed=input$seed, contrast = matrix(c(1,1,0,0,0,1,0,1),nrow=2,byrow=TRUE),
            #                                           nobs = nobs.seq, nobsTX = nobs.seqTX, dprior = p.n,
            #                                           b = cut.seq, pow = power.seq, pn=p.n, pa = p.a,
            #                                           phi=c(input$e1n,1-input$e2n), cutstart=cut.start,
            #                                           nsim = input$numOfSimForTiralSetting,
            #                                           err1 = input$err1, maxresp,resp_list_r)

            ########################################################
            r0 = input$e1n
            t0 = input$e2n
            r1 = input$e1a
            t1 = input$e2a

            interm.eff = nobs.seq #
            interm.tox = nobs.seqTX #
            print(unique(sort(c(interm.eff, interm.tox)))) #
            scen = hypotheses_ind(r0, r1, t0, t1) ##independent
            ##########################################################
            if(input$priorspec){
              prior = scen[1,]
            }else{
              # pi = c(p00, p01, p10, p11)
              prior= c(1-input$ETprior1-input$ETprior2+input$ETprior3,
                       input$ETprior2-input$ETprior3,
                       input$ETprior1-input$ETprior3,
                       input$ETprior3)*input$ETpriorSS

            }


            oc.mat.efftox.temp= grid_search_cpp(interm.eff=interm.eff,
                                                interm.tox=interm.tox,
                                                scenario=scen,
                                                r0 = r0, r1= r1,
                                                t0=t0, t1= t1,
                                                lambda_r=cut.seq,
                                                lambda_t= cut.seq,
                                                gamma_space=power.seq,
                                                typeI01=input$err_eff, typeI10=input$err_tox,
                                                typeI00=input$err_all,
                                                true_opt = input$truet1e,
                                                prior = scen[1,],
                                                cppfunction=GridSearchBiRcpp) # return a vector

            if(input$e3n == (input$e1n*input$e2n)& input$e3a == (input$e1a*input$e2a) ){

              oc.mat.efftox =oc.mat.efftox.temp
              print("Independent")

            }else{
              print("Not independent")
              optim.temp=oc.mat.efftox.temp[1:3]
              lambda_E_space = seq(optim.temp[1]-0.2, min(1,optim.temp[1]+0.05 ),by=0.01)
              lambda_T_space = seq(optim.temp[2]-0.2, min(1,optim.temp[2]+0.05 ),by=0.01)
              gamma_space = seq(optim.temp[3]-0.2, min(1,optim.temp[3]+0.1 ),by=0.02)


              scen = hypotheses_corr(r0,r1,t0,t1,PA_ET=input$e3a, PN_ET = input$e3n)
              print(scen)
              # joint probabilities of four possible outcomes
              # (no efficacy and no toxicity, no efficacy and toxicity, efficacy and no toxicity, both efficacy and toxicity)

              if(input$priorspec){
                prior = scen[1,]
              }else{
                # pi = c(p00, p01, p10, p11)
                prior= c(1-input$ETprior1-input$ETprior2+input$ETprior3,
                         input$ETprior2-input$ETprior3,
                         input$ETprior1-input$ETprior3,
                         input$ETprior3)*input$ETpriorSS

              }
              search_space = bond.unique(interm.eff=interm.eff, interm.tox=interm.tox, scenario=scen,
                                         lambda_r = lambda_E_space,
                                         lambda_t = lambda_T_space,
                                         gamma_space= gamma_space,
                                         r0=r0,t0=t0,
                                         prior = scen[1,])

              oc.mat.efftox=optimization_corr(interm.eff=interm.eff,
                                              interm.tox=interm.tox,
                                              search.temp = search_space,
                                              scenario=scen,
                                              r0=r0, t0=t0,
                                              typeI01=input$err_eff,
                                              typeI10=input$err_tox,
                                              typeI00=input$err_all,
                                              lpow=0.3,
                                              ncore=10)

            }

            if(is.null(oc.mat.efftox) | length(oc.mat.efftox)==0) {
              globalPara$error="The power is too small (<5%) under the current setting. Please increase the type I error or the sample size."
              globalPara$calculated=FALSE
              globalPara$errFlag=TRUE
              globalPara$powertext=""
              globalPara$boundary=c()
              globalPara$prior=c()
            } else {
              optim=oc.mat.efftox
              #print(oc.mat.efftox[1:10,])
              # oc.mat.efftoxS = oc.mat.efftox[oc.mat.efftox[,3]<=input$err1,]
              # oc.mat.efftoxL = oc.mat.efftox[oc.mat.efftox[,3]>input$err1,]
              #
              # #print(oc.mat.efftoxS)
              # #print(oc.mat.efftoxL)
              #
              # if(is.null(oc.mat.efftoxS) | length(oc.mat.efftoxS)==0) {
              #   globalPara$error="The power is too small (<5%) under the current setting. Please increase the type I error or the sample size."
              #   globalPara$calculated=FALSE
              #   globalPara$errFlag=TRUE
              #   globalPara$powertext=""
              #   globalPara$boundary=c()
              #   globalPara$prior=c()
              # }
              #
              # else{
              #
              #   if(is.null(oc.mat.efftoxL) | length(oc.mat.efftoxL)==0){
              #     if(length(oc.mat.efftoxS)==4) optim <- oc.mat.efftoxS
              #     else optim=oc.mat.efftoxS[oc.mat.efftoxS[,4]==max(oc.mat.efftoxS[,4]),]
              #
              #   }else{
              #     if(length(oc.mat.efftoxS)==4) optimS <- oc.mat.efftoxS
              #     else optimS <- oc.mat.efftoxS[oc.mat.efftoxS[,4]==max(oc.mat.efftoxS[,4]),]
              #     if(length(oc.mat.efftoxL)==4) optimL <- oc.mat.efftoxL
              #     else optimL <- oc.mat.efftoxL[oc.mat.efftoxL[,3]==min(oc.mat.efftoxL[,3]),]
              #
              #     if(is.matrix(optimS)) optimS=optimS[1,]
              #     if(is.matrix(optimL)) optimL=optimL[1,]
              #
              #     #print(optimS)
              #     #print(optimL)
              #
              #     if(optimS[4]>=optimL[4]) {optim=optimS}
              #     else if(abs(optimS[3]-input$err1)>=abs(optimL[3]-input$err1) & (!input$truet1e)) {optim=optimL}
              #     else {optim=optimS}
              #   }
              #
              #   if(is.matrix(optim)) optim=optim[1,]


              globalPara$lambda_E <- optim[1]
              globalPara$lambda_T <- optim[2]
              globalPara$gamma <- optim[3]
              globalPara$t1e <- optim[6]
              # globalPara$power= optim[6]
              globalPara$a01= optim[4]
              globalPara$a10= optim[5]
              globalPara$a00= optim[7]

              cat("Prior is " ,prior)
              boundary2=get_boundarycpp(interm.eff=interm.eff,
                                        interm.tox=interm.tox,
                                        prior=prior,
                                        lambda_e = optim[1], lambda_t =optim[2],
                                        gamma = optim[3], r0=r0, t0=t0)


              globalPara$boundary.eff = boundary2[[1]]
              globalPara$boundary.tox = boundary2[[2]]

              bound_temp=get_boundary_table(interm.eff=interm.eff,
                                            interm.tox=interm.tox,
                                            prior=prior,
                                            lambda_e = optim[1], lambda_t =optim[2],
                                            gamma = optim[3], r0=r0, t0=t0)

              globalPara$power = BOP2_TE.OC(interm.eff=nobs.seq,
                                            interm.tox=nobs.seqTX,
                                            boundary.eff=globalPara$boundary.eff,
                                            boundary.tox=globalPara$boundary.tox,
                                            scenario=scen[4,])[1]
              # globalPara$boundary =cbind(unique(sort(c(interm.eff, interm.tox))),as.character(bound_temp[[1]]),
              #                            as.character(bound_temp[[2]]))

              globalPara$pts_e = optim[8]
              globalPara$pts_t = optim[9]

              globalPara$boundary =bound_temp
              # dimnames(globalPara$boundary) <- list(NULL,c("# patients treated","Stop if # response <="," OR # toxicity >="))
              globalPara$calculated=TRUE

              # temp2 = Calc$Getoc_Rcpp_Tox_Eff(input$seed, input$numOfSimForTiralSetting, matrix(c(1,1,0,0,0,1,0,1),
              #                                                                                   nrow=2,byrow=TRUE), nobs.seq,nobs.seqTX, optim[1], optim[2], prior, p.a, c(input$e1n,1-input$e2n), maxresp,resp_list_r);
              # globalPara$power=temp2[[3]]
              globalPara$powertext <- paste("The power of this trial is: ",round(globalPara$power,4),sep="")
              isSBcalculated$ET <- TRUE
              globalPara$prior=prior

            }
          }
          else{globalPara$error="Error: The value of Prob(Eff & Tox) in Prior Specification is invalid (either too large or too small).
                         Given the values of Prob(Eff) and Prob(Tox), Prob(Eff & Tox) must satisfy the constraint Prob(Eff) + Prob(Tox) -1 <= Prob(Eff & Tox) <= min{Prob(Eff), Prob(Tox)}.
                         The left-hand side of the constraint is to ensure that Prob(no Eff & no Tox)>=0."
          globalPara$calculated=FALSE
          globalPara$errFlag=TRUE
          globalPara$powertext=""
          globalPara$prior=c()
          globalPara$boundary=c()
          }
        }
        else {globalPara$error="Error: The value in Prior Specification is missing or invalid (either too large or too small)."
        globalPara$calculated=FALSE
        globalPara$errFlag=TRUE
        globalPara$powertext=""
        globalPara$prior=c()
        globalPara$boundary=c()}
      }
      else {globalPara$error="Error: The value of Pr(Eff & Tox) is invalid (either too large or too small).
                   Given the values of Pr(Eff) and Pr(Tox), Pr(Eff & Tox) must satisfy the constraint Pr(Eff) + Pr(Tox) -1 < Pr(Eff & Tox) < min{Pr(Eff), Pr(Tox)}.
                   The left-hand side of the constraint is to ensure that Pr(no Eff & no Tox)>0."
      globalPara$calculated=FALSE
      globalPara$errFlag=TRUE
      globalPara$powertext=""
      globalPara$prior=c()
      globalPara$boundary=c()}
    } else {globalPara$error="Error: The value of Pr(Eff)/Pr(Tox) in alternative hypothesis must be greater/smaller than the null."
    globalPara$calculated=FALSE
    globalPara$errFlag=TRUE
    globalPara$powertext=""
    globalPara$prior=c()
    globalPara$boundary=c()}
  }
  else {globalPara$error="Error: The value of null/alternative hypothesis is missing, too small (negative) or greater than 1"
  globalPara$calculated=FALSE
  globalPara$errFlag=TRUE
  globalPara$powertext=""
  globalPara$prior=c()
  globalPara$boundary=c()}

  return(globalPara)

}


DoStoppingBoundaries_minN = function( input, globalPara, Calc, nobs.seq, nobs.seqTX, isSBcalculated, session )
{
  #browser();
  if(input$e1n>0.0001 & input$e1a>0.0001 & input$e1n<1 & input$e1a<1 &
     input$e2n>0.0001 & input$e2a>0.0001 & input$e2n<1 & input$e2a<1 &
     input$e3n>0.0001 & input$e3a>0.0001 & input$e3n<1 & input$e3a<1 &
     !is.na(input$e1n) & !is.na(input$e1a) & !is.na(input$e2n) & !is.na(input$e2a) & !is.na(input$e3n) & !is.na(input$e3a)){

    if(input$e1n<=input$e1a & input$e2n>=input$e2a) {

      if(input$e1n+input$e2n-input$e3n<1 & input$e1a+input$e2a-input$e3a<1 &
         input$e3n<min(input$e1n,input$e2n) & input$e3a<min(input$e1a,input$e2a)){

        if((input$ETprior1>=0 & input$ETprior1<=1 & input$ETprior2>=0 & input$ETprior2<=1 & input$ETprior3>=0 & input$ETprior3<=1 &
            !is.na(input$ETprior1) & !is.na(input$ETprior2) & !is.na(input$ETprior3) & !is.na(input$ETpriorSS) & input$ETpriorSS>0)| input$priorspec){

          if(input$ETprior1+input$ETprior2-input$ETprior3 <= 1  &  input$ETprior3 <= min(input$ETprior1,input$ETprior2)){

            # cut.seq = sort(c(seq(0.5,0.95,by=0.025),seq(0.96,0.99,0.01),seq(0.9-input$err1,1-input$err1,0.005)))
            cut.seq = rev(c(seq(0.5,0.8,0.025),seq(0.81,0.96,0.01)))
            power.seq = log(seq(1,0.5,by=-0.025))/log(0.5)
            p.n = c(input$e3n,input$e1n-input$e3n,input$e2n-input$e3n,1-input$e1n-input$e2n+input$e3n)
            ##input$e3n: Pr(Eff & Tox), input$e1n: Pr(Eff), input$e2n: Pr(Tox)
            ##input$e1n-input$e3n: Pr(Eff & no Tox)
            ##input$e2n-input$e3n: Pr(no Eff & Tox)
            ##1-input$e1n-input$e2n+input$e3n: Pr(no Eff & no Tox)
            p.a = c(input$e3a,input$e1a-input$e3a,input$e2a-input$e3a,1-input$e1a-input$e2a+input$e3a)
            # if(input$priorspec){
            #    prior = p.n
            # }else{
            #    prior = c(input$ETprior3*input$ETpriorSS,(input$ETprior1-input$ETprior3)*input$ETpriorSS,
            #              (input$ETprior2-input$ETprior3)*input$ETpriorSS,(1-input$ETprior1-input$ETprior2+input$ETprior3)*input$ETpriorSS)
            # }

            oc.mat.efftox = c()

            # oc.mat.efftox = Calc$GridSearchToxEffRcpp(seed=input$seed, contrast = matrix(c(1,1,0,0,0,1,0,1),nrow=2,byrow=TRUE),
            #                                           nobs = nobs.seq, nobsTX = nobs.seqTX, dprior = p.n,
            #                                           b = cut.seq, pow = power.seq, pn=p.n, pa = p.a,
            #                                           phi=c(input$e1n,1-input$e2n), cutstart=cut.start,
            #                                           nsim = input$numOfSimForTiralSetting,
            #                                           err1 = input$err1, maxresp,resp_list_r)

            r0 = input$e1n
            t0 = input$e2n
            r1 = input$e1a
            t1 = input$e2a

            interm.eff = nobs.seq
            interm.tox = nobs.seqTX
            print(unique(sort(c(interm.eff, interm.tox))))
            scen = hypotheses_ind(r0, r1, t0, t1)
            #
            if(input$priorspec){
              prior = scen[1,]
            }else{
              # pi = c(p00, p01, p10, p11)
              prior= c(1-input$ETprior1-input$ETprior2+input$ETprior3,
                       input$ETprior2-input$ETprior3,
                       input$ETprior1-input$ETprior3,
                       input$ETprior3)*input$ETpriorSS

            }


            # oc.mat.efftox.temp= grid_search_cpp_minN(interm.eff=interm.eff,
            #                                     interm.tox=interm.tox,
            #                                     scenario=scen,
            #                                     r0 = r0, r1= r1,
            #                                     t0=t0, t1= t1,
            #                                     lambda_r=cut.seq,
            #                                     lambda_t= cut.seq,
            #                                     gamma_space=power.seq,
            #                                     typeI01=input$err_eff, typeI10=input$err_tox,
            #                                     typeI00=input$err_all,
            #                                     power_e = input$power_eff,
            #                                     power_t = input$power_tox,
            #                                     true_opt = input$truet1e,
            #                                     prior = scen[1,],
            #                                     cppfunction=GridSearchBiRcpp) # return a vector
            # oc.mat.efftox.temp= grid_search_cpp(interm.eff=interm.eff,
            #                                     interm.tox=interm.tox,
            #                                     scenario=scen,
            #                                     r0 = r0, r1= r1,
            #                                     t0=t0, t1= t1,
            #                                     lambda_r=cut.seq,
            #                                     lambda_t= cut.seq,
            #                                     gamma_space=power.seq,
            #                                     typeI01=input$err_eff, typeI10=input$err_tox,
            #                                     typeI00=input$err_all,
            #                                     true_opt = input$truet1e,
            #                                     prior = scen[1,],
            #                                     cppfunction=GridSearchBiRcpp) # return a vector

            if(input$e3n == (input$e1n*input$e2n)& input$e3a == (input$e1a*input$e2a) ){

              oc.mat.efftox =oc.mat.efftox.temp
              print("Independent")

            }else{
              print("Not independent")
              # optim.temp=oc.mat.efftox.temp[1:3]
              # lambda_E_space = seq(optim.temp[1]-0.2, min(1,optim.temp[1]+0.05 ),by=0.01)
              # lambda_T_space = seq(optim.temp[2]-0.2, min(1,optim.temp[2]+0.05 ),by=0.01)
              # gamma_space = seq(optim.temp[3]-0.2, min(1,optim.temp[3]+0.1 ),by=0.02)
              lambda_E_space = cut.seq
              lambda_E2_space = cut.seq
              lambda_T_space = cut.seq
              lambda_T2_space = cut.seq
              gamma_space = power.seq

              scen = hypotheses_corr(r0,r1,t0,t1,PA_ET=input$e3a, PN_ET = input$e3n)

              if(input$priorspec){
                prior = scen[1,]
              }else{
                # pi = c(p00, p01, p10, p11)
                prior= c(1-input$ETprior1-input$ETprior2+input$ETprior3,
                         input$ETprior2-input$ETprior3,
                         input$ETprior1-input$ETprior3,
                         input$ETprior3)*input$ETpriorSS

              }
              search_space = bond.unique(interm.eff=interm.eff, interm.tox=interm.tox, scenario=scen,
                                         lambda_r = lambda_E_space,
                                         lambda_t = lambda_T_space,
                                         gamma_space= gamma_space,
                                         r0=r0,t0=t0,
                                         prior = scen[1,])
              # search_space = bond.unique_2lambda(interm.eff=interm.eff, interm.tox=interm.tox, scenario=scen,
              #                            lambda_r = lambda_E_space,
              #                            lambda_t = lambda_T_space,
              #                            lambda_r2 = lambda_E_space,
              #                            lambda_t2 = lambda_T_space,
              #                            gamma_space= gamma_space,
              #                            r0=r0,t0=t0,
              #                            prior = scen[1,])

              oc.mat.efftox=optimization_corr_minN(interm.eff=interm.eff,
                                              interm.tox=interm.tox,
                                              search.temp = search_space,
                                              scenario=scen,
                                              r0=r0, t0=t0,
                                              typeI01=input$err_eff,
                                              typeI10=input$err_tox,
                                              typeI00=input$err_all,
                                              power11 = 0.8,
                                              lpow=0.3,
                                              ncore=10)

            }

            if(is.null(oc.mat.efftox) | length(oc.mat.efftox)==0) {
              globalPara$error="The power is too small (<5%) under the current setting. Please increase the type I error or the sample size."
              globalPara$calculated=FALSE
              globalPara$errFlag=TRUE
              globalPara$powertext=""
              globalPara$boundary=c()
              globalPara$prior=c()
            } else {
              optim=oc.mat.efftox
              #print(oc.mat.efftox[1:10,])
              # oc.mat.efftoxS = oc.mat.efftox[oc.mat.efftox[,3]<=input$err1,]
              # oc.mat.efftoxL = oc.mat.efftox[oc.mat.efftox[,3]>input$err1,]
              #
              # #print(oc.mat.efftoxS)
              # #print(oc.mat.efftoxL)
              #
              # if(is.null(oc.mat.efftoxS) | length(oc.mat.efftoxS)==0) {
              #   globalPara$error="The power is too small (<5%) under the current setting. Please increase the type I error or the sample size."
              #   globalPara$calculated=FALSE
              #   globalPara$errFlag=TRUE
              #   globalPara$powertext=""
              #   globalPara$boundary=c()
              #   globalPara$prior=c()
              # }
              #
              # else{
              #
              #   if(is.null(oc.mat.efftoxL) | length(oc.mat.efftoxL)==0){
              #     if(length(oc.mat.efftoxS)==4) optim <- oc.mat.efftoxS
              #     else optim=oc.mat.efftoxS[oc.mat.efftoxS[,4]==max(oc.mat.efftoxS[,4]),]
              #
              #   }else{
              #     if(length(oc.mat.efftoxS)==4) optimS <- oc.mat.efftoxS
              #     else optimS <- oc.mat.efftoxS[oc.mat.efftoxS[,4]==max(oc.mat.efftoxS[,4]),]
              #     if(length(oc.mat.efftoxL)==4) optimL <- oc.mat.efftoxL
              #     else optimL <- oc.mat.efftoxL[oc.mat.efftoxL[,3]==min(oc.mat.efftoxL[,3]),]
              #
              #     if(is.matrix(optimS)) optimS=optimS[1,]
              #     if(is.matrix(optimL)) optimL=optimL[1,]
              #
              #     #print(optimS)
              #     #print(optimL)
              #
              #     if(optimS[4]>=optimL[4]) {optim=optimS}
              #     else if(abs(optimS[3]-input$err1)>=abs(optimL[3]-input$err1) & (!input$truet1e)) {optim=optimL}
              #     else {optim=optimS}
              #   }
              #
              #   if(is.matrix(optim)) optim=optim[1,]


              globalPara$lambda_E <- optim$parameter[1]
              globalPara$lambda_T <- optim$parameter[2]
              globalPara$gamma <- optim$parameter[3]
              # globalPara$t1e <- optim$prob[6]
              # globalPara$power= optim[6]
              globalPara$a01= optim$prob["a01"]
              globalPara$a10= optim$prob["a10"]
              globalPara$a00= optim$prob["a00"]

              globalPara$power = optim$prob["power"]

              cat("Prior is " ,prior)
              boundary2=get_boundarycpp(interm.eff=interm.eff,
                                        interm.tox=interm.tox,
                                        prior=prior,
                                        lambda_e = globalPara$lambda_E, lambda_t = globalPara$lambda_T,
                                        gamma = globalPara$gamma, r0=r0, t0=t0)


              globalPara$boundary.eff = boundary2[[1]]
              globalPara$boundary.tox = boundary2[[2]]

              bound_temp=get_boundary_table(interm.eff=interm.eff,
                                            interm.tox=interm.tox,
                                            prior=prior,
                                            lambda_e = globalPara$lambda_E, lambda_t =globalPara$lambda_T,
                                            gamma = globalPara$gamma, r0=r0, t0=t0)

              # globalPara$power = BOP2_TE.OC(interm.eff=nobs.seq,
              #                               interm.tox=nobs.seqTX,
              #                               boundary.eff=globalPara$boundary.eff,
              #                               boundary.tox=globalPara$boundary.tox,
              #                               scenario=scen[4,])[1]
              # globalPara$boundary =cbind(unique(sort(c(interm.eff, interm.tox))),as.character(bound_temp[[1]]),
              #                            as.character(bound_temp[[2]]))

              globalPara$samplesize = optim$SampleSize

              globalPara$boundary =bound_temp
              # dimnames(globalPara$boundary) <- list(NULL,c("# patients treated","Stop if # response <="," OR # toxicity >="))
              globalPara$calculated=TRUE

              # temp2 = Calc$Getoc_Rcpp_Tox_Eff(input$seed, input$numOfSimForTiralSetting, matrix(c(1,1,0,0,0,1,0,1),
              #                                                                                   nrow=2,byrow=TRUE), nobs.seq,nobs.seqTX, optim[1], optim[2], prior, p.a, c(input$e1n,1-input$e2n), maxresp,resp_list_r);
              # globalPara$power=temp2[[3]]
              globalPara$powertext <- paste("The power of this trial is: ",round(globalPara$power,4),sep="")
              isSBcalculated$ET <- TRUE
              globalPara$prior=prior

            }
          }
          else{globalPara$error="Error: The value of Prob(Eff & Tox) in Prior Specification is invalid (either too large or too small).
                         Given the values of Prob(Eff) and Prob(Tox), Prob(Eff & Tox) must satisfy the constraint Prob(Eff) + Prob(Tox) -1 <= Prob(Eff & Tox) <= min{Prob(Eff), Prob(Tox)}.
                         The left-hand side of the constraint is to ensure that Prob(no Eff & no Tox)>=0."
          globalPara$calculated=FALSE
          globalPara$errFlag=TRUE
          globalPara$powertext=""
          globalPara$prior=c()
          globalPara$boundary=c()
          }
        }
        else {globalPara$error="Error: The value in Prior Specification is missing or invalid (either too large or too small)."
        globalPara$calculated=FALSE
        globalPara$errFlag=TRUE
        globalPara$powertext=""
        globalPara$prior=c()
        globalPara$boundary=c()}
      }
      else {globalPara$error="Error: The value of Pr(Eff & Tox) is invalid (either too large or too small).
                   Given the values of Pr(Eff) and Pr(Tox), Pr(Eff & Tox) must satisfy the constraint Pr(Eff) + Pr(Tox) -1 < Pr(Eff & Tox) < min{Pr(Eff), Pr(Tox)}.
                   The left-hand side of the constraint is to ensure that Pr(no Eff & no Tox)>0."
      globalPara$calculated=FALSE
      globalPara$errFlag=TRUE
      globalPara$powertext=""
      globalPara$prior=c()
      globalPara$boundary=c()}
    } else {globalPara$error="Error: The value of Pr(Eff)/Pr(Tox) in alternative hypothesis must be greater/smaller than the null."
    globalPara$calculated=FALSE
    globalPara$errFlag=TRUE
    globalPara$powertext=""
    globalPara$prior=c()
    globalPara$boundary=c()}
  }
  else {globalPara$error="Error: The value of null/alternative hypothesis is missing, too small (negative) or greater than 1"
  globalPara$calculated=FALSE
  globalPara$errFlag=TRUE
  globalPara$powertext=""
  globalPara$prior=c()
  globalPara$boundary=c()}

  return(globalPara)

}
# end function DoStoppingBoundaries
