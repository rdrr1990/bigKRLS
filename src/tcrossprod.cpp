#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigTCrossProd(const Mat<T>& A, const Mat<T>& B, Mat<T> out) {

  Mat<T> transposed = trans(B);
  out = A * transposed;

}

// [[Rcpp::export]]
void BigTCrossProd(SEXP pA, SEXP pB, SEXP pOut) {

  XPtr<BigMatrix> xpA(pA);
  XPtr<BigMatrix> xpB(pB);
  XPtr<BigMatrix> xpOut(pOut);

  xBigTCrossProd(
    arma::Mat<double>((double *)xpA->matrix(), xpA->nrow(), xpA->ncol(), false),
    arma::Mat<double>((double *)xpB->matrix(), xpB->nrow(), xpB->ncol(), false),
    arma::Mat<double>((double *)xpOut->matrix(), xpOut->nrow(), xpOut->ncol(), false)
  );
}
