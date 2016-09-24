#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigTempKernel(const Mat<T>& inBigMatNew, const Mat<T>& inBigMatOld, Mat<T> outBigMat, Mat<T> sigma) {

int U = inBigMatNew.n_rows;
int V = inBigMatOld.n_rows;

for(int i = 0; i < U; i++){
  for(int j = 0; j < V; j++){
    double val = exp(-1 * sum(pow((inBigMatNew.row(i) - inBigMatOld.row(j)),2))/sigma(0,0));
    outBigMat(i,j) = val;
    }
  }
}

// [[Rcpp::export]]
void BigTempKernel(SEXP pInBigMatNew, SEXP pInBigMatOld, SEXP pOutBigMat, SEXP pSigma) {
  XPtr<BigMatrix> xpInMatNew(pInBigMatNew);
  XPtr<BigMatrix> xpInMatOld(pInBigMatOld);
  XPtr<BigMatrix> xpOutMat(pOutBigMat);
  XPtr<BigMatrix> xpSigma(pSigma);

  xBigTempKernel(
    arma::Mat<double>((double *)xpInMatNew->matrix(), xpInMatNew->nrow(), xpInMatNew->ncol(), false),
    arma::Mat<double>((double *)xpInMatOld->matrix(), xpInMatOld->nrow(), xpInMatOld->ncol(), false),
    arma::Mat<double>((double *)xpOutMat->matrix(), xpOutMat->nrow(), xpOutMat->ncol(), false),
    arma::Mat<double>((double *)xpSigma->matrix(), xpSigma->nrow(), xpSigma->ncol(), false)
  );
}
