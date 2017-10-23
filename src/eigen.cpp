#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigEigen(const Mat<T>& A, const double eigtrunc, Mat<T> vals, Mat<T> vecs) {
  vec valsVec;
  
  // later, might implement partial eigen setup here
  arma::eig_sym(valsVec, vecs, A);
  //if(eigtrunc == A.n_rows){  
  //}
  //else{
  //  sp_mat sparse_data(A);
  //  arma::eigs_sym(valsVec, vecs, sparse_data, eigtrunc);
  //}
  
  vals = valsVec;
  vals = flipud(vals);
  vecs = fliplr(vecs);
}

// [[Rcpp::export]]
void BigEigen(SEXP pA, const double Eigtrunc, SEXP pValBigMat, SEXP pVecBigMat) {
  XPtr<SharedMemoryBigMatrix> xpMat(pA);
  XPtr<SharedMemoryBigMatrix> xpValMat(pValBigMat);
  XPtr<SharedMemoryBigMatrix> xpVecMat(pVecBigMat);

  xBigEigen(
      arma::Mat<double>((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false),
      Eigtrunc,
      arma::Mat<double>((double *)xpValMat->matrix(), xpValMat->nrow(), xpValMat->ncol(), false),
      arma::Mat<double>((double *)xpVecMat->matrix(), xpVecMat->nrow(), xpVecMat->ncol(), false)
    );
}
