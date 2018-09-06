#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


double cpp_posterior_density(const arma::vec& THETAi, const arma::vec& Yi, const arma::vec& MUi, const arma::mat& invS, const arma::vec& TAU, const arma::mat& LAMBDA, const int J, const int K) {

  int j;
  double lodd_j;
  double b;
  // PYi
  arma::vec LT = LAMBDA*THETAi;
  arma::vec PYi(J);
  for (j = 0; j < J; j++){
    lodd_j =  LT(j) - TAU(j);
    if(lodd_j<0){
      b = 0.0;
    } else {
      b = lodd_j;
    }
    PYi(j) = exp(lodd_j-b)/(exp(-b) + exp(lodd_j-b));
  }

  double ll = 0.0;
  for (j = 0; j < J; j++){
    if (Yi(j)==1L){
      ll += log(PYi(j));
    }
    if (Yi(j)==0L){
      ll += log(1.0-PYi(j));
    }
  }

  //Prior distriubtion
  arma::vec dMUi = THETAi-MUi;
  double twoprior = as_scalar(dMUi.t()*invS*dMUi);

  // Return result
  double dpost = -1.0*ll + 0.5*twoprior;
  return dpost;
}
