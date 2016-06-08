#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigElementwise(const Mat<T>& inMatX, const Mat<T>& inMatY, Mat<T> outBigMat) {
  outBigMat = inMatX % inMatY;
}

// [[Rcpp::export]]
void BigElementwise(SEXP pInMatX, SEXP pInMatY, SEXP poutBigMat) {
  XPtr<BigMatrix> xpInMatX(pInMatX);
  XPtr<BigMatrix> xpInMatY(pInMatY);
  XPtr<BigMatrix> xpoutBigMat(poutBigMat);

  xBigElementwise(
    arma::Mat<double>((double *)xpInMatX->matrix(), xpInMatX->nrow(), xpInMatX->ncol(), false),
    arma::Mat<double>((double *)xpInMatY->matrix(), xpInMatY->nrow(), xpInMatY->ncol(), false),
    arma::Mat<double>((double *)xpoutBigMat->matrix(), xpoutBigMat->nrow(), xpoutBigMat->ncol(), false)
  );
}
