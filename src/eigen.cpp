#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigEigen(const Mat<T>& A, const double Neig, Mat<T> vals, Mat<T> vecs) {
  
  vec temp_vals;
  
  if(Neig < A.n_rows){
    mat temp_vecs;
    sp_mat sparseA(A);
    arma::eigs_sym(temp_vals, temp_vecs, sparseA, Neig);
    vecs = temp_vecs; 
  }else{
    arma::eig_sym(temp_vals, vecs, A);
  }
  vals = temp_vals;
  
  vals = flipud(vals);
  vecs = fliplr(vecs);
}

// [[Rcpp::export]]
void BigEigen(SEXP pA, const double Neig, SEXP pValBigMat, SEXP pVecBigMat) {
  
  XPtr<SharedMemoryBigMatrix> xpMat(pA);   // A is typically the Kernel
  XPtr<SharedMemoryBigMatrix> xpValMat(pValBigMat);
  XPtr<SharedMemoryBigMatrix> xpVecMat(pVecBigMat);
  
  xBigEigen(
    arma::Mat<double>((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false),
    Neig,
    arma::Mat<double>((double *)xpValMat->matrix(), xpValMat->nrow(), xpValMat->ncol(), false),
    arma::Mat<double>((double *)xpVecMat->matrix(), xpVecMat->nrow(), xpVecMat->ncol(), false)
  );
}
