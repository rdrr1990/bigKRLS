#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigMultDiag(const Mat<T>& A, const arma::rowvec diag, Mat<T> out) {

  int W = A.n_cols;
  for(int i = 0; i < W; i++){
    out.col(i) = A.col(i)*diag[i];
    if(i % 501 == 0){
      Rcpp::checkUserInterrupt();
      Rprintf(".");
    }
  }
}

// [[Rcpp::export]]
void BigMultDiag(SEXP pA, const arma::rowvec diag, SEXP pOut) {

  XPtr<BigMatrix> xpA(pA);
  XPtr<BigMatrix> xpOut(pOut);

  xBigMultDiag(
    arma::Mat<double>((double *)xpA->matrix(), xpA->nrow(), xpA->ncol(), false),
    diag,
    arma::Mat<double>((double *)xpOut->matrix(), xpOut->nrow(), xpOut->ncol(), false)
  );
}
