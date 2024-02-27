// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double obj(double l2, double lr, arma::mat Q, double dn, double YTY,
           arma::rowvec YTZ, arma::mat ZTZ, arma::colvec ZTY) {
  arma::mat B = l2 * Q + lr * arma::eye(Q.n_rows, Q.n_cols);
  arma::mat BZTZ = B + ZTZ;
  arma::mat matprod = (YTZ * arma::inv(BZTZ) * ZTY);
  double mle = dn * log(abs(YTY - matprod(0,0))) + log(abs(arma::det(BZTZ)))
    - log(abs(arma::det(B)));
  return  mle;
}


