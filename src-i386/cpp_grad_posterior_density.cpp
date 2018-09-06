#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec cpp_grad_posterior_density(const arma::vec& THETAi, const arma::vec& Yi, const arma::vec& MUi, const arma::mat& invS, const arma::vec& TAU, const arma::mat& LAMBDA, const int J, const int K) {

  int j, k;

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

  arma::vec ddll(K);
  for (j = 0; j < J; j++){
    if (Yi(j)==1L || Yi(j)==0L){

      double YmP = (Yi(j)-PYi(j));

      for (k = 0; k<K; k++){
          ddll(k) += YmP*LAMBDA(j,k);
      }

    }
  }

  //Prior distriubtion
  arma::vec dMUi = THETAi-MUi;
  arma::vec ddprior = invS*dMUi;

  // Return result
  arma::vec ddpost = ddprior - ddll;
  return ddpost;
}
