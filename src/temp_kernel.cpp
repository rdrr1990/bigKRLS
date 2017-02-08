#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigTempKernel(const Mat<T>& A, const Mat<T>& B, Mat<T> out, const double sigma) {

int U = A.n_rows;
int V = B.n_rows;

for(int i = 0; i < U; i++){
  for(int j = 0; j < V; j++){
    double val = exp(-1 * sum(pow((A.row(i) - B.row(j)),2))/sigma);
    out(i,j) = val;
    }
  // checking for user interrupt on the outer loop
  if(i % 501 == 0){
    Rcpp::checkUserInterrupt();
    Rprintf("*");
  }
  }
}

// [[Rcpp::export]]
void BigTempKernel(SEXP pA, SEXP pB, SEXP pOut, const double sigma) {
  XPtr<BigMatrix> xpA(pA);
  XPtr<BigMatrix> xpB(pB);
  XPtr<BigMatrix> xpOut(pOut);
  
  xBigTempKernel(
    arma::Mat<double>((double *)xpA->matrix(), xpA->nrow(), xpA->ncol(), false),
    arma::Mat<double>((double *)xpB->matrix(), xpB->nrow(), xpB->ncol(), false),
    arma::Mat<double>((double *)xpOut->matrix(), xpOut->nrow(), xpOut->ncol(), false),
    sigma
  );
}
