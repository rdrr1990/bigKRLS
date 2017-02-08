#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
List xBigSolveForc(Mat<T> Eigenvectors, const colvec Eigenvalues, 
                   const colvec y, const double lambda){
  
  double Le = 0;
  
  int N = Eigenvectors.n_rows;
  
  colvec Ginv_diag(N); Ginv_diag.zeros();
  colvec coeffs(N); coeffs.zeros();
  
  Eigenvectors = trans(Eigenvectors);
  
  for(int i = 0; i < N; i++){
    colvec g(i+1);
    
    mat temp_eigen(Eigenvectors.memptr(), N, i+1, false);
    g = (Eigenvectors.col(i).t()/(Eigenvalues + lambda)) * temp_eigen;
    
    Ginv_diag[i] = g[i];
    coeffs(span(0,i-1)) += g * y[i];
    coeffs[i] += sum(g * y(span(0,i)));
    
    // checking for user interrupt
    if(i % 501 == 0){
      Rcpp::checkUserInterrupt();
      Rprintf("*");
    }
  }
  
  Eigenvectors = trans(Eigenvectors);
  
  for(int i = 0; i < N; i++){
    Le += pow((coeffs[i]/Ginv_diag[i]), 2);
  }
  
  List out(2);
  out[0] = Le;
  out[1] = coeffs;
  
  return out;
}

// [[Rcpp::export]]
List BigSolveForc(SEXP pEigenvectors, const arma::colvec Eigenvalues, 
                      const arma::colvec y, const double lambda) {
  
  XPtr<BigMatrix> xpEigenvectors(pEigenvectors);
  
  List out = xBigSolveForc(arma::Mat<double>((double *)xpEigenvectors->matrix(), 
                                             xpEigenvectors->nrow(), 
                                             xpEigenvectors->ncol(), false),
                                             Eigenvalues, y, lambda
  );
  return out;
}
