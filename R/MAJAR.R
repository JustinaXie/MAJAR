require(MASS)
require(dplyr)
require(Rcpp)
require(RcppArmadillo)
require(parallel)
require(TMB)
require(optimParallel)
# Rcpp::sourceCpp("src/file.cpp")

Get_FDR<-function(ppr){
  P<-length(ppr)
  lfdr<-data.frame(index=1:P,Z=1-ppr)
  lfdr_order<-lfdr[order(lfdr$Z,decreasing = FALSE),]
  lfdr_order$FDR<-cumsum(lfdr_order$Z)/(1:P)
  tmp<-lfdr_order[order(lfdr_order$index,decreasing = FALSE),]
  tmp_FDR<-tmp$FDR
  return(tmp_FDR)
}
Inverse.Matrix<-function(mat){
  tryCatch({
    a = mat[1,1]; b = mat[1,2]; c = mat[2,1]; d = mat[2,2]
    s <- 1/(a*d-b*c)
    re <- s*matrix(c(d,-b,-c,a), nrow = 2, ncol = 2, byrow = TRUE)
    return(re)
  },error=function(x){
    warning("General Inverse Used")
    ginv(mat)
  })
}


MAJAR<-function(betajk_G, sjk2_G,
                betajk_GT, sjk2_GT,
                p=0.1,#Rj~Bernoulli(pi)
                lambda=0.1,#Ojk~Bernoulli(lambda)
                tau2=0.01,
                alpha=5,
                single_study=FALSE,
                rho=0.1,
                rel.eps=1e-8,
                param.eps = 1e-4,
                minIter = 100,
                maxIter=1e4,
                parcores=1,
                tmp=10^6,
                transform=TRUE,
                fast_implement=TRUE,
                verbose=FALSE){
  
  start <- Sys.time()
  
  betajk_G<-as.matrix(betajk_G)
  betajk_GT<-as.matrix(betajk_GT)
  sjk2_G<-as.matrix(sjk2_G)
  sjk2_GT<-as.matrix(sjk2_GT)
  
  MM<-nrow(betajk_G)#number of SNPs
  KK<-ncol(betajk_G)#number of study
  
  #Initialize mu_G and mu_GT, this is correct
  mu_G<-sapply(1:MM, function(j){
    sum(betajk_G[j,]/sjk2_G[j,])/(sum(1/sjk2_G[j,]) + 1/tau2)
  })
  
  mu_GT<-sapply(1:MM, function(j){
    sum(betajk_GT[j,]/sjk2_GT[j,])/(sum(1/sjk2_GT[j,]) + 1/tau2)
  })
  
  #Some known values
  os22_G<-sapply(1:MM, function(j){
    sum(1/sjk2_G[j,])
  })
  
  os22_GT<-sapply(1:MM, function(j){
    sum(1/sjk2_GT[j,])
  })
  
  b2s2_G<-sapply(1:MM, function(j){
    sum(betajk_G[j,]^2/sjk2_G[j,])
  })
  
  b2s2_GT<-sapply(1:MM, function(j){
    sum(betajk_GT[j,]^2/sjk2_GT[j,])
  })
  
  bs2_G<-sapply(1:MM, function(j){
    sum(betajk_G[j,]/sjk2_G[j,])
  })
  bs2_GT<-sapply(1:MM, function(j){
    sum(betajk_GT[j,]/sjk2_GT[j,])
  })
  
  ll<-c()
  lambda_old<-p_old<-tau2_old<-alpha_old<-rho_old<-1
  
  #intialize Vj_inv
  Vj11_inv<-os22_G+1/(tau2*(1-rho^2))
  Vj21_inv<--rho/(tau2*(1-rho^2))
  Vj22_inv<-os22_GT+1/(tau2*(1-rho^2))
  
  for (k in 1:maxIter) {
    if(verbose){
      print(paste("round", k, "start!"))
    }
    
    #E-step
    #R_j=1, need to consider the prior according to mamba R package
    
    #Rcpp: this is correct
    ll_p1<-unlist(mclapply(1:MM, function(j){
      ll_R1_j(bs2_j_G=bs2_G[j], bs2_j_GT=bs2_GT[j],
              os22_j_G=os22_G[j], os22_j_GT=os22_GT[j],
              b2s2_j_G=b2s2_G[j],  b2s2_j_GT=b2s2_GT[j],
              lambda,alpha,rho,tau2,
              sjk2_j_G=sjk2_G[j,],sjk2_j_GT=sjk2_GT[j,])
    }, mc.cores = parcores))-KK*log(2*pi)
    
    #R_j=0, O_jk=1, this is correct
    #Rcpp: this is correct
    ll_p2<-unlist(mclapply(1:MM, function(j){
      llbR0_j(betajk_j_G = betajk_G[j,],
              betajk_j_GT = betajk_GT[j,],
              sjk2_j_G = sjk2_G[j,],
              sjk2_j_GT = sjk2_GT[j,],
              alpha = alpha,
              lambda=lambda)
    }, mc.cores = parcores))
    
    
    #update Rj
    # b_gamma<-sapply(1:MM,function(j) max(log(p)+ll_p1[j],log(1-p)+ll_p2[j]))
    gamma_j<-exp(log(p)+ll_p1)/(exp(log(p)+ll_p1)+exp(log(1-p)+ll_p2))
    
    anyNA(gamma_j)
    
    #update Ojk|Rj=0
    #Rcpp version: this is correct
    delta_jk<-matrix(unlist(mclapply(1:MM,function(j){
      deltis(betajk_j_G=betajk_G[j,], betajk_j_GT=betajk_GT[j,],
             sjk2_j_G=sjk2_G[j,], sjk2_j_GT=sjk2_GT[j,],
             alpha=alpha, lambda=lambda)
    },mc.cores = parcores)),nrow = MM,byrow = TRUE)
    if(verbose){
      print("E-step completes")
    }
    
    #M-step
    p_old<-p
    lambda_old<-lambda
    alpha_old<-alpha
    tau2_old<-tau2
    rho_old<-rho
    
    #update pi
    p<-sum(gamma_j)/(MM+MM/tmp)
    
    #update lambda, this is correct
    lambda<-sum((1-gamma_j)*rowSums(delta_jk))/(sum(1-gamma_j)*KK)
    
    #update alpha, this is correct
    alpha<- sum((1-gamma_j)*rowSums(delta_jk*(betajk_G^2/(2*sjk2_G)+betajk_GT^2/(2*sjk2_GT))))/
      sum((1-gamma_j)*rowSums(delta_jk))
    if(single_study){
      alpha<-max(1,alpha)
    }
    if(fast_implement){
      #update rho and tau
      #C++ functions written by TMB
      if(!transform){
        nllk <- MakeADFun(data=list(bs2_G=bs2_G,bs2_GT=bs2_GT,os22_G=os22_G,os22_GT=os22_GT,gamma_j=gamma_j),
                          parameters=list(tau2=tau2,rho=rho),DLL="tau2_rho_f_v1",silent=TRUE)
        fit <- nlminb(start = nllk$par, objective=nllk$fn,gradient =nllk$gr,hessian = nllk$he,
                      lower =c(1e-64,-0.999), upper =c(1,0.999),
                      control=list(trace=0))
        
        tau2<-fit$par[1]
        # tau2<-max(tau2, 1e-17)
        rho<-fit$par[2]
        # rho<-ifelse(rho>0,min(rho,0.999),max(-0.999,rho))
      }else{
        nllk <- MakeADFun(data=list(bs2_G=bs2_G,bs2_GT=bs2_GT,os22_G=os22_G,os22_GT=os22_GT,gamma_j=gamma_j),
                          parameters=list(logTau2=log(tau2),atanhrho=atanh(rho)),DLL="tau2_rho_f_v2",silent=TRUE)
        fit <- nlminb(start = nllk$par, objective=nllk$fn,gradient =nllk$gr,hessian = nllk$he,
                      lower =c(-Inf,-Inf), upper =c(Inf,Inf),
                      control=list(trace=0))
        tau2<-unname(exp(fit$par[1]))
        #numerical issues
        tau2<-max(tau2, 1e-17)
        rho<-unname(tanh(fit$par[2]))
        #numerical issues
        rho<-ifelse(rho>0,min(rho,0.999),max(-0.999,rho))
      }
    }
    if(!fast_implement){
      # update rho and tau
      # rho_est,tau2_est
      ll_tau_rho<-function(params){
        rho_est<-params[1]
        tau2_est<-params[2]
        Vj11_inv<-os22_G+1/(tau2_est*(1-rho_est^2))
        Vj21_inv<--rho_est/(tau2_est*(1-rho_est^2))
        Vj22_inv<-os22_GT+1/(tau2_est*(1-rho_est^2))
        ll_j<-unlist(lapply(1:MM, function(j){
          V_inv<-matrix(c(Vj11_inv[j],Vj21_inv,Vj21_inv,Vj22_inv[j]),2,2)
          # m<-Inverse.Matrix(V_inv)%*%c(sum(betajk_G[j,]/sjk2_G[j,]),sum(betajk_GT[j,]/sjk2_GT[j,]))
          b_j<-c(bs2_G[j],bs2_GT[j])
          re<--2*log(tau2_est)-log(1-rho_est^2)-log(Vj11_inv[j]*Vj22_inv[j]-Vj21_inv^2)+t(b_j)%*%Inverse.Matrix(V_inv)%*%b_j
          re
        }))
        return(sum(gamma_j*ll_j))
      }
      if(verbose){
        print("M-step: tau2 and rho estimation start")
      }
      cl<-makeCluster(parcores)
      if(verbose){
        print(paste0(cl," of cores assigned!"))
      }
      setDefaultCluster(cl=cl)
      clusterExport(cl,c('Inverse.Matrix'))
      
      params_est<-optimParallel(par=c(rho,tau2), fn=ll_tau_rho,
                                lower = c(-0.999, 1e-16), upper = c(0.999, 1),
                                control = list(fnscale = -1))$par
      stopCluster(cl)
      
      rho<-params_est[1]
      tau2<-params_est[2]
      if(verbose){
        print("M-step completes")
      }
      
      #numerical issues
      rho<-ifelse(rho>0,min(rho,0.999),max(-0.999,rho))
      tau2<-max(tau2, 1e-17)
      
    }
    
    if(is.na(tau2)||is.na(rho)) {
      print("Numerical issues happened, the algorithm stops")
      break
    }
    if(verbose){
      print("M-step completes")
    }
    
    #update mu
    Vj11_inv<-os22_G+1/(tau2*(1-rho^2))
    Vj21_inv<--rho/(tau2*(1-rho^2))
    Vj22_inv<-os22_GT+1/(tau2*(1-rho^2))
    m<-matrix(unlist(lapply(1:MM, function(j){
      V_inv <- matrix(c(Vj11_inv[j],Vj21_inv,Vj21_inv,Vj22_inv[j]),2,2)
      tmp<-Inverse.Matrix(V_inv)%*%c(sum(betajk_G[j,]/sjk2_G[j,]),sum(betajk_GT[j,]/sjk2_GT[j,]))
      as.vector(tmp)
    })),nrow = MM,ncol=2,byrow = TRUE)
    
    mu_G<-m[,1]
    mu_GT<-m[,2]
    
    #Rcpp: this is correct
    llbR1<-unlist(mclapply(1:MM, function(j){
      ll_R1_j(bs2_j_G=bs2_G[j], bs2_j_GT=bs2_GT[j],
              os22_j_G=os22_G[j], os22_j_GT=os22_GT[j],
              b2s2_j_G=b2s2_G[j],  b2s2_j_GT=b2s2_GT[j],
              lambda,alpha,rho,tau2,
              sjk2_j_G=sjk2_G[j,],sjk2_j_GT=sjk2_GT[j,])
    }, mc.cores = parcores))-KK*log(2*pi)
    
    #R_j=0, O_jk=1, this is correct
    #Rcpp: this is correct
    llbR0<-unlist(mclapply(1:MM, function(j){
      llbR0_j(betajk_j_G = betajk_G[j,],
              betajk_j_GT = betajk_GT[j,],
              sjk2_j_G = sjk2_G[j,],
              sjk2_j_GT = sjk2_GT[j,],
              alpha = alpha,
              lambda=lambda)
    }, mc.cores = parcores))
    
    
    ll[k] <- sum(unlist(mclapply(1:MM, function(j){
      logsumexp(c(log(p) + llbR1[j], log(1-p) + llbR0[j]))
    }, mc.cores = parcores)))
    
    h2 <- tau2*(MM*p)
    #Debug
    if(verbose){
      print(paste0("tau2 estimate is ",tau2))
      print(paste0("rho estimate is ",rho))
      print(paste0("pi estimate is ",p))
      print(paste0("lambda estimate is ",lambda))
      print(paste0("alpha estimate is ",alpha))
    }
    
    #convergence criteria
    if(k>1){
      ll_rel <- (ll[k]-ll[k-1])/ll[k]
      ll_diff <- ll[k]-ll[k-1]
      if(verbose){
        print(paste0("relative increase of ll is ",ll_rel))
        print(paste0("increase of ll is ",ll_diff))
      }
      
      if((ll_rel<rel.eps && ll_diff<0.01)||
         # (abs(p-p_old)/(0.1+p)<param.eps)||
         (ll_diff < 0.001 && k > minIter) ||
         (ll_rel < rel.eps && k > minIter)){
        break
      }
    }
  }
  end <- Sys.time()
  
  Time <- difftime(end,start,units = "mins")
  
  FDR<-Get_FDR(gamma_j)
  
  return(list(ppr=gamma_j,
              outlier.mat=delta_jk,
              p=p,
              lambda=lambda,
              alpha=alpha,
              tau2=tau2,
              h2=h2,
              rho=rho,
              mu_G=mu_G,
              mu_GT=mu_GT,
              post.means_G=mu_G*gamma_j,
              post.means_GT=mu_GT*gamma_j,
              FDR = FDR,
              ll = ll[k],
              Time = Time))
}
