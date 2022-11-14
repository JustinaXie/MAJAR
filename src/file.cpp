#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double logsumexp(NumericVector x) {
  return(log(sum(exp(x - max(x)))) + max(x));
}

// [[Rcpp::export]]
double llbR0_j(NumericVector betajk_j_G, NumericVector betajk_j_GT,
               NumericVector sjk2_j_G, NumericVector sjk2_j_GT,
               double lambda, double alpha) {
  int k_j=betajk_j_G.size();
  NumericVector cohortll0(k_j);
  for(int k=0; k < k_j; k++) {
    NumericVector llb0_inds(2L);
    llb0_inds[0]=log(1-lambda) + R::dnorm(betajk_j_G[k], 0.0, sqrt(sjk2_j_G[k]), true)+R::dnorm(betajk_j_GT[k], 0.0, sqrt(sjk2_j_GT[k]), true) ;
    llb0_inds[1]=log(lambda) + R::dnorm(betajk_j_G[k], 0.0, sqrt(alpha*sjk2_j_G[k]), true)+ R::dnorm(betajk_j_GT[k], 0.0, sqrt(alpha*sjk2_j_GT[k]), true);
    cohortll0[k] = logsumexp(llb0_inds) ;
  }
  return sum(cohortll0) ;
}

//[[Rcpp::export]]
NumericVector deltis(NumericVector betajk_j_G, NumericVector betajk_j_GT,
                     NumericVector sjk2_j_G,NumericVector sjk2_j_GT,
                     double lambda, double alpha) {
  int k_j=betajk_j_G.size();
  NumericVector cohortll0(k_j);
  for(int k=0; k<k_j; k++) {
    if(!Rcpp::NumericVector::is_na(betajk_j_G[k])){
      NumericVector llb0_inds(2L);
      llb0_inds[0]=log(1-lambda) + R::dnorm(betajk_j_G[k], 0.0, sqrt(sjk2_j_G[k]), true)+ R::dnorm(betajk_j_GT[k], 0.0, sqrt(sjk2_j_GT[k]), true) ;
      llb0_inds[1]=log(lambda) + R::dnorm(betajk_j_G[k], 0.0, sqrt(alpha*sjk2_j_G[k]), true)+ R::dnorm(betajk_j_GT[k], 0.0, sqrt(alpha*sjk2_j_GT[k]), true);
      cohortll0[k] = exp(llb0_inds[1]) / exp(logsumexp(llb0_inds));
      if(Rcpp::NumericVector::is_na(cohortll0[k])){
        double m=which_max(llb0_inds) ; 
        cohortll0[k]= m ;
      }
    }  else {
      cohortll0[k]=NA_REAL ;
    }
  }
  return cohortll0  ;
}

//[[Rcpp::export]]
double ll_R1_j(double bs2_j_G, double bs2_j_GT,double os22_j_G, double os22_j_GT,double b2s2_j_G, double b2s2_j_GT,
                double lambda, double alpha, double rho, double tau2,
                arma::vec sjk2_j_G,arma::vec sjk2_j_GT){
  int k_j=sjk2_j_G.size();
  arma::vec ll_j(1);
  arma::mat V(2,2);
  for(int k=0; k<k_j; k++) {
        // printf("The %dth study begins \n", k+1);
         arma::vec b_j = {bs2_j_G, bs2_j_GT};
        //Construct the Vj_inv matrix
        arma::mat V_inv = {{os22_j_G + 1/(tau2*(1-std::pow(rho,2))),- rho/(tau2*(1-std::pow(rho,2)))},
        {- rho/(tau2*(1-std::pow(rho,2))),os22_j_GT + 1/(tau2*(1-std::pow(rho,2)))}};
        // arma:: mat V = arma::inv(V_inv);
        try{
          V = arma::inv(V_inv); //cause an exception to throw
        }
        catch (const char* msg)
        {
          cerr << msg << endl;
          V = arma::pinv(V_inv);
        };
        ll_j = -0.5*sum(log(sjk2_j_G))-0.5*sum(log(sjk2_j_GT))-
              log(tau2)-0.5*log(1-std::pow(rho,2))-0.5*(b2s2_j_G+b2s2_j_GT)-
              0.5*log(arma::det(V_inv))+0.5*arma::trans(b_j)*V*b_j;
      }
  double ll_j_re = as_scalar(ll_j);
  return(ll_j_re) ;
}
// ll_p1<-sapply(1:MM, function(j){
//   b_j<-c(bs2_G[j,],bs2_GT[j,])
//   V_inv<-matrix(c(Vj11_inv[j],Vj21_inv,Vj21_inv,Vj22_inv[j]),2,2)
// # m<-Inverse.Matrix(V_inv)%*%c(sum(betajk_G[j,]/sjk2_G[j,]),sum(betajk_GT[j,]/sjk2_GT[j,]))
//   ll_re<- -0.5*sum(log(sjk2_G[j,]))-0.5*sum(log(sjk2_GT[j,]))-
//     log(tau2)-0.5*log(1-rho^2)-0.5*(b2s2_G[j]+b2s2_GT[j])-
//     0.5*log(det(V_inv))+0.5*t(b_j)%*%Inverse.Matrix(V_inv)%*%b_j
//   ll_re
// })-KK*log(2*pi)
//   anyNA(ll_p1)