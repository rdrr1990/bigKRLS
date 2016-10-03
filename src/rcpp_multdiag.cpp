#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigMultDiag(const Mat<T>& toMult, const Mat<T>& vecMat, Mat<T> outMat) {

int W = toMult.n_cols;

for(int i = 0; i < W; i++){
  outMat.col(i) = toMult.col(i)*vecMat(i,0);
  }
}

// [[Rcpp::export]]
void BigMultDiag(SEXP pToMult, SEXP pVecMat, SEXP pOutMat) {

  XPtr<BigMatrix> xpToMult(pToMult);
  XPtr<BigMatrix> xpVecMat(pVecMat);
  XPtr<BigMatrix> xpOutMat(pOutMat);

  xBigMultDiag(
    arma::Mat<double>((double *)xpToMult->matrix(), xpToMult->nrow(), xpToMult->ncol(), false),
    arma::Mat<double>((double *)xpVecMat->matrix(), xpVecMat->nrow(), xpVecMat->ncol(), false),
    arma::Mat<double>((double *)xpOutMat->matrix(), xpOutMat->nrow(), xpOutMat->ncol(), false)
  );
}
