require(data.table)
require(parallel)
require(MASS)
require(dplyr)
require(mvtnorm)

generate_data_DGP<-function(M=NULL,K=NULL,data = c("random"),
                            p=NULL,lambda=NULL,h2=NULL,rho=NULL,alpha=NULL,
                            model=NULL, I2 = NULL){

  tau2 <- h2/(M*p)
  Rj<-rbinom(M, size=1, prob=p)
  muj<-mvrnorm(M,mu=rep(0,2),Sigma =matrix(c(tau2,rho*tau2,rho*tau2,tau2),2,2))

  if(data=="random"){
    data(res2ref_cpd, package="mamba",envir=environment())
    resref<-sqrt(res2ref_cpd$res2)
    set.seed(01)
    res.sd1<-lapply(1:M, function(j){sample(resref, K, replace = TRUE)})
    set.seed(02)
    res.sd2<-lapply(1:M, function(j){sample(resref, K, replace = TRUE)})

    sjk2_G<-lapply(1:M, function(j){res.sd1[[j]]^2 / rep(M,K)})
    sjk2_GT<-lapply(1:M, function(j){res.sd2[[j]]^2 / rep(M,K)})

    sjk2_G<-matrix(unlist(sjk2_G), ncol=K, byrow=TRUE)
    sjk2_GT<-matrix(unlist(sjk2_GT), ncol=K, byrow=TRUE)
  }


  Ojk<-lapply(1:M, function(j){rbinom(K, size=1, prob=lambda)})

  if(model == "MAJAR"){

    betajk_G<-lapply(1:M, function(j){
      Rj[j] * rnorm(K, muj[j,1], sd=sqrt(sjk2_G[j,])) +
        (1-Rj[j]) * ((1-Ojk[[j]])*rnorm(K, 0, sqrt(sjk2_G[j,])) +
                       (Ojk[[j]])*rnorm(K, 0, sqrt(alpha*sjk2_G[j,])))})

    betajk_GT<-lapply(1:M, function(j){
      Rj[j] * rnorm(K, muj[j,2], sd=sqrt(sjk2_GT[j,])) +
        (1-Rj[j]) * ((1-Ojk[[j]])*rnorm(K, 0, sqrt(sjk2_GT[j,])) +
                       (Ojk[[j]])*rnorm(K, 0, sqrt(alpha*sjk2_GT[j,])))})

    Ojk<-matrix(unlist(Ojk), ncol=K, byrow=TRUE)

    betajk_G<-matrix(unlist(betajk_G), ncol=K, byrow=TRUE)
    betajk_GT<-matrix(unlist(betajk_GT), ncol=K, byrow=TRUE)

    colnames(betajk_G) <- colnames(betajk_GT) <- paste0("Study",1:K)
    rownames(betajk_G) <- rownames(betajk_GT) <- paste0("SNP",1:M)
    betajk_G = round(betajk_G, 5); betajk_GT = round(betajk_GT, 5)
  }

  ## FE
  if(model == "FE"){
    muj_ind<-muj
    muj_ind[,1]<-Rj*muj[,1]
    muj_ind[,2]<-Rj*muj[,2]

    betajk_G<-matrix(unlist(lapply(1:M, function(j){
      rnorm(K, muj_ind[j,1], sd=sqrt(sjk2_G[j,]))
    })),ncol = K,byrow = TRUE)

    betajk_GT<-matrix(unlist(lapply(1:M, function(j){
      rnorm(K, muj_ind[j,2], sd=sqrt(sjk2_GT[j,]))
    })),ncol = K,byrow = TRUE)

    colnames(betajk_G) <- colnames(betajk_GT) <- paste0("Study",1:K)
    rownames(betajk_G) <- rownames(betajk_GT) <- paste0("SNP",1:M)
    betajk_G = round(betajk_G, 5); betajk_GT = round(betajk_GT, 5)
  }

  ## RE
  if(model == "RE"){
    muj_ind<-muj
    muj_ind[,1]<-Rj*muj[,1]
    muj_ind[,2]<-Rj*muj[,2]

    x2_g <- sapply(1:M, function(j){
      ((K-1)*sum(1/sjk2_G[j,]))/(sum(1/sjk2_G[j,])^2 - sum(1/sjk2_G[j,]^2))
    })
    x2_gt <- sapply(1:M, function(j){
      ((K-1)*sum(1/sjk2_GT[j,]))/(sum(1/sjk2_GT[j,])^2 - sum(1/sjk2_GT[j,]^2))
    })
    #Initialize value for I2
    # I2 = 0.3
    w2_g <- sapply(1:M, function(j){
      max(0,I2/(1-I2) * x2_g[j])
    })
    w2_gt <- sapply(1:M, function(j){
      max(0,I2/(1-I2) * x2_gt[j])
    })

    eta_g<-lapply(1:M, function(j){
      rnorm(K, mean = muj_ind[j,1], sd = sqrt(w2_g[j]))
    })

    eta_gt<-lapply(1:M, function(j){
      rnorm(K, mean = muj_ind[j,2], sd = sqrt(w2_gt[j]))
    })

    betajk_G<-matrix(unlist(lapply(1:M, function(j){
      rnorm(K, eta_g[[j]], sd=sqrt(sjk2_G[j,]))
    })),ncol = K,byrow = TRUE)

    betajk_GT<-matrix(unlist(lapply(1:M, function(j){
      rnorm(K, eta_gt[[j]], sd=sqrt(sjk2_GT[j,]))
    })),ncol = K,byrow = TRUE)

    colnames(betajk_G) <- colnames(betajk_GT) <- paste0("Study",1:K)
    rownames(betajk_G) <- rownames(betajk_GT) <- paste0("SNP",1:M)
    betajk_G = round(betajk_G, 5); betajk_GT = round(betajk_GT, 5)
  }

  ## RE2
  if(model == "RE2"){
    x2_g <- sapply(1:M, function(j){
      ((K-1)*sum(1/sjk2_G[j,]))/(sum(1/sjk2_G[j,])^2 - sum(1/sjk2_G[j,]^2))
    })
    x2_gt <- sapply(1:M, function(j){
      ((K-1)*sum(1/sjk2_GT[j,]))/(sum(1/sjk2_GT[j,])^2 - sum(1/sjk2_GT[j,]^2))
    })
    # I2 = 0.3
    w2_g <- sapply(1:M, function(j){
      max(0,I2/(1-I2) * x2_g[j])
    })
    w2_gt <- sapply(1:M, function(j){
      max(0,I2/(1-I2) * x2_gt[j])
    })

    eta_g<-lapply(1:M, function(j){
      Rj[j] * rnorm(K, muj[j,1], sd=sqrt(w2_g[j]))
    })
    eta_gt<-lapply(1:M, function(j){
      Rj[j] * rnorm(K, muj[j,2], sd=sqrt(w2_gt[j]))
    })

    betajk_G<-matrix(unlist(lapply(1:M, function(j){
      rnorm(K, eta_g[[j]], sd=sqrt(sjk2_G[j,]))
    })),ncol = K,byrow = TRUE)

    betajk_GT<-matrix(unlist(lapply(1:M, function(j){
      rnorm(K, eta_gt[[j]], sd=sqrt(sjk2_GT[j,]))
    })),ncol = K,byrow = TRUE)

    colnames(betajk_G) <- colnames(betajk_GT) <- paste0("Study",1:K)
    rownames(betajk_G) <- rownames(betajk_GT) <- paste0("SNP",1:M)
    betajk_G = round(betajk_G, 5); betajk_GT = round(betajk_GT, 5)
  }

  ## BE
  if(model=="BE"){
    Nej<-sapply(1:M, function(j){ Rj[j]*sample(1:K, 1) })
    bej<-lapply(1:M, function(j){
      be<-rep(0, K)
      if(Nej[j] > 0){
        be[sample(1:K, Nej[j])] <- 1
      }
      return(be)
    })

    betajk_G<-matrix(unlist(lapply(1:M, function(j){bej[[j]] * rnorm(K, muj[j,1], sd=sqrt(sjk2_G[j,])) +
        (1-bej[[j]]) *rnorm(K, 0, sqrt(sjk2_G[j,]))
    })),ncol=K,byrow = TRUE)

    betajk_GT<-matrix(unlist(lapply(1:M, function(j){bej[[j]] * rnorm(K, muj[j,2], sd=sqrt(sjk2_GT[j,])) +
        (1-bej[[j]]) *rnorm(K, 0, sqrt(sjk2_GT[j,]))
    })),ncol=K,byrow = TRUE)

    colnames(betajk_G) <- colnames(betajk_GT) <- paste0("Study",1:K)
    rownames(betajk_G) <- rownames(betajk_GT) <- paste0("SNP",1:M)
    betajk_G = round(betajk_G, 5); betajk_GT = round(betajk_GT, 5)

  }

  # P-value
  meta.z_G<-sapply(1:M,function(j) sum(betajk_G[[j]]/sjk2_G[j,]) / sqrt(sum(1/sjk2_G[j,])))
  meta.z_GT <- sapply(1:M,function(j) sum(betajk_GT[[j]]/sjk2_GT[j,]) / sqrt(sum(1/sjk2_GT[j,])))

  Pvalue_G <- matrix(unlist(lapply(1:M,function(j){
    2*(1-pnorm(abs(betajk_G[j,]/sqrt(sjk2_G[j,]))))
  })),ncol=K,byrow = TRUE)
  Pvalue_GT <- matrix(unlist(lapply(1:M,function(j){
    2*(1-pnorm(abs(betajk_GT[j,]/sqrt(sjk2_GT[j,]))))
  })),ncol=K,byrow = TRUE)

  return(list(betajk_G=betajk_G,
              betajk_GT=betajk_GT,
              meta.z_G=meta.z_G,
              meta.z_GT=meta.z_GT,
              Pvalue_G=Pvalue_G,
              Pvalue_GT=Pvalue_GT,
              Rj=Rj,
              Ojk=Ojk,
              muj_G=muj[,1],
              muj_GT=muj[,2],
              sjk2_G=sjk2_G,
              sjk2_GT=sjk2_GT))
}

Get_FDR<-function(ppr){
  P<-length(ppr)
  lfdr<-data.frame(index=1:P,Z=1-ppr)
  lfdr_order<-lfdr[order(lfdr$Z,decreasing = FALSE),]
  lfdr_order$FDR<-cumsum(lfdr_order$Z)/(1:P)
  tmp<-lfdr_order[order(lfdr_order$index,decreasing = FALSE),]
  tmp_FDR<-tmp$FDR
  return(tmp_FDR)
}
