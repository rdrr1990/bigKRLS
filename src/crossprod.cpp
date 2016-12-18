#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigCrossProd(const Mat<T>& A, const Mat<T>& B, Mat<T> out) {

  Mat<T> transposed = trans(A);
  out = transposed * B;

}

// [[Rcpp::export]]
void BigCrossProd(SEXP pA, SEXP pB, SEXP pOut) {

  XPtr<BigMatrix> xpA(pA);
  XPtr<BigMatrix> xpB(pB);
  XPtr<BigMatrix> xpOut(pOut);

  xBigCrossProd(
    arma::Mat<double>((double *)xpA->matrix(), xpA->nrow(), xpA->ncol(), false),
    arma::Mat<double>((double *)xpB->matrix(), xpB->nrow(), xpB->ncol(), false),
    arma::Mat<double>((double *)xpOut->matrix(), xpOut->nrow(), xpOut->ncol(), false)
  );

}
