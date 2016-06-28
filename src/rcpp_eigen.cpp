#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigEigen(const Mat<T>& inBigMat, const Mat<T>& inEigtrunc, Mat<T> valsBigMat, Mat<T> vecsBigMat) {
  vec valsBigMatVec;
  int eigtrunc = as_scalar(inEigtrunc);
  
  // later, might implement partial eigen setup here
  arma::eig_sym(valsBigMatVec, vecsBigMat, inBigMat);
  //if(eigtrunc == inBigMat.n_rows){  
  //}
  //else{
  //  sp_mat sparse_data(inBigMat);
  //  arma::eigs_sym(valsBigMatVec, vecsBigMat, sparse_data, eigtrunc);
  //}
  
  valsBigMat = valsBigMatVec;
  valsBigMat = flipud(valsBigMat);
  vecsBigMat = fliplr(vecsBigMat);
}

// [[Rcpp::export]]
void BigEigen(SEXP pInBigMat, SEXP pInEigTrunc, SEXP pValBigMat, SEXP pVecBigMat) {
  XPtr<BigMatrix> xpMat(pInBigMat);
  XPtr<BigMatrix> xpInEigTrunc(pInEigTrunc);
  XPtr<BigMatrix> xpValMat(pValBigMat);
  XPtr<BigMatrix> xpVecMat(pVecBigMat);

  xBigEigen(
      arma::Mat<double>((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false),
      arma::Mat<double>((double *)xpInEigTrunc->matrix(), xpInEigTrunc->nrow(), xpInEigTrunc->ncol(), false),
      arma::Mat<double>((double *)xpValMat->matrix(), xpValMat->nrow(), xpValMat->ncol(), false),
      arma::Mat<double>((double *)xpVecMat->matrix(), xpVecMat->nrow(), xpVecMat->ncol(), false)
    );
}
