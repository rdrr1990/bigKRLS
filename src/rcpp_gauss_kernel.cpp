#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigGaussKernel(const Mat<T>& inBigMat, Mat<T> outBigMat, Mat<T> sigma) {

int W = inBigMat.n_rows;

for(int i = 0; i < W; i++){
for(int j = i; j < W; j++){
    double dist_val = exp(-1 * sum(pow((inBigMat.row(i) - inBigMat.row(j)),2))/sigma(0,0));
    outBigMat(j,i) = dist_val;
    outBigMat(i,j) = dist_val;
    }
  }
}

// [[Rcpp::export]]
void BigGaussKernel(SEXP pInBigMat, SEXP pOutBigMat, SEXP pSigma) {
  XPtr<BigMatrix> xpMat(pInBigMat);
  XPtr<BigMatrix> xpOutMat(pOutBigMat);
  XPtr<BigMatrix> xpSigma(pSigma);

  xBigGaussKernel(
    arma::Mat<double>((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false),
    arma::Mat<double>((double *)xpOutMat->matrix(), xpOutMat->nrow(), xpOutMat->ncol(), false),
    arma::Mat<double>((double *)xpSigma->matrix(), xpSigma->nrow(), xpSigma->ncol(), false)
  );
}
