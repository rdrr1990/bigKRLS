#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigGaussKernel(const Mat<T>& X, Mat<T> out, const double sigma) {

int N = X.n_rows;

for(int i = 0; i < N; ++i){
  for(int j = i; j < N; ++j){
    double similarity = exp(-1 * sum(pow((X.row(i) - X.row(j)),2))/sigma);
    out(j,i) = similarity;
    out(i,j) = similarity;
  }
  // checking for user interrupt on the outer loop
  if(i % 501 == 0){
    Rcpp::checkUserInterrupt();
    Rprintf(".");
  }  
  }
}

// [[Rcpp::export]]
void BigGaussKernel(SEXP pA, SEXP pOut, const double sigma) {
  XPtr<BigMatrix> xpMat(pA);
  XPtr<BigMatrix> xpOut(pOut);

  xBigGaussKernel(
    arma::Mat<double>((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false),
    arma::Mat<double>((double *)xpOut->matrix(), xpOut->nrow(), xpOut->ncol(), false),
    sigma
  );
}
