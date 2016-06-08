#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigEigen(const Mat<T>& inBigMat, Mat<T> valsBigMat, Mat<T> vecsBigMat) {
  vec valsBigMatVec;
  
  arma::eig_sym(valsBigMatVec, vecsBigMat, inBigMat);

  valsBigMat = valsBigMatVec;
}

// [[Rcpp::export]]
void BigEigen(SEXP pInBigMat, SEXP pValBigMat, SEXP pVecBigMat) {
  XPtr<BigMatrix> xpMat(pInBigMat);
  XPtr<BigMatrix> xpValMat(pValBigMat);
  XPtr<BigMatrix> xpVecMat(pVecBigMat);

  xBigEigen(
      arma::Mat<double>((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false),
      arma::Mat<double>((double *)xpValMat->matrix(), xpValMat->nrow(), xpValMat->ncol(), false),
      arma::Mat<double>((double *)xpVecMat->matrix(), xpVecMat->nrow(), xpVecMat->ncol(), false)
    );
}
