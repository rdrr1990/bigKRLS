#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigGaussKernel(const Mat<T>& A, Mat<T> out, const double sigma) {

int W = A.n_rows;

for(int i = 0; i < W; i++){
for(int j = i; j < W; j++){
    double dist_val = exp(-1 * sum(pow((A.row(i) - A.row(j)),2))/sigma);
    out(j,i) = dist_val;
    out(i,j) = dist_val;

    }
  // checking for user interrupt on the outer loop
  if(i % 501 == 0){
    Rcpp::checkUserInterrupt();
    Rprintf("*");
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
