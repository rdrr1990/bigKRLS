#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigCrossProd(const Mat<T>& A, const Mat<T>& B, Mat<T> out) {
  out = trans(A) * B;
}

// [[Rcpp::export]]
void BigCrossProd(SEXP pA, SEXP pB, SEXP pOut) {

  XPtr<BigMatrix> xpA(pA);
  XPtr<BigMatrix> xpB(pB);
  XPtr<BigMatrix> xpOut(pOut);

  xBigCrossProd(
    arma::Mat<double>((double *)xpA->matrix(), xpA->nrow(), xpA->ncol(), false),
    arma::Mat<double>((double *)xpB->matrix(), xpB->nrow(), xpB->ncol(), false),
    arma::Mat<double>((double *)xpOut->matrix(), xpOut->nrow(), xpOut->ncol(), false)
  );
}

template <typename T>
void xBigXtX(const Mat<T>& A, Mat<T> out) {
  out = trans(A) * A;
}

// [[Rcpp::export]]
void BigXtX(SEXP pA, SEXP pOut) {
  
  XPtr<BigMatrix> xpA(pA);
  XPtr<BigMatrix> xpOut(pOut);
  
  xBigXtX(
    arma::Mat<double>((double *)xpA->matrix(), xpA->nrow(), xpA->ncol(), false),
    arma::Mat<double>((double *)xpOut->matrix(), xpOut->nrow(), xpOut->ncol(), false)
  );
  
}


template <typename T>
void xBigTCrossProd(const Mat<T>& A, const Mat<T>& B, Mat<T> out) {
  out = A * trans(B);
}

// [[Rcpp::export]]
void BigTCrossProd(SEXP pA, SEXP pB, SEXP pOut) {

  XPtr<SharedMemoryBigMatrix> xpA(pA);
  XPtr<SharedMemoryBigMatrix> xpB(pB);
  XPtr<SharedMemoryBigMatrix> xpOut(pOut);

  xBigTCrossProd(
    arma::Mat<double>((double *)xpA->matrix(), xpA->nrow(), xpA->ncol(), false),
    arma::Mat<double>((double *)xpB->matrix(), xpB->nrow(), xpB->ncol(), false),
    arma::Mat<double>((double *)xpOut->matrix(), xpOut->nrow(), xpOut->ncol(), false)
  );
}

template <typename T>
void xBigXXt(const Mat<T>& A, Mat<T> out) {
  out = A * trans(A);
}

// [[Rcpp::export]]
void BigXXt(SEXP pA, SEXP pOut) {
  
  XPtr<SharedMemoryBigMatrix> xpA(pA);
  XPtr<SharedMemoryBigMatrix> xpOut(pOut);
  
  xBigXXt(
    arma::Mat<double>((double *)xpA->matrix(), xpA->nrow(), xpA->ncol(), false),
    arma::Mat<double>((double *)xpOut->matrix(), xpOut->nrow(), xpOut->ncol(), false)
  );
}



