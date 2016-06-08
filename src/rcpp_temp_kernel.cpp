#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigTempKernel(const Mat<T>& inBigMat, Mat<T> outBigMat, Mat<T> sigma) {

int W = inBigMat.n_rows/3;

for(int i = 0; i < W; i++){
  for(int j = 0; j < W*2; j++){
    double val = exp(-1 * sum(pow((inBigMat.row(W*2 + i) - inBigMat.row(j)),2))/sigma(0,0));
    outBigMat(j,i) = val;
    }
  }
}

// [[Rcpp::export]]
void BigTempKernel(SEXP pInBigMat, SEXP pOutBigMat, SEXP pSigma) {
  XPtr<BigMatrix> xpMat(pInBigMat);
  XPtr<BigMatrix> xpOutMat(pOutBigMat);
  XPtr<BigMatrix> xpSigma(pSigma);

  xBigTempKernel(
    arma::Mat<double>((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false),
    arma::Mat<double>((double *)xpOutMat->matrix(), xpOutMat->nrow(), xpOutMat->ncol(), false),
    arma::Mat<double>((double *)xpSigma->matrix(), xpSigma->nrow(), xpSigma->ncol(), false)
  );
}
