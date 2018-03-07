#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
List xBigSolveForc(Mat<T> Eigenvectors, const colvec Eigenvalues, 
                   const colvec y, const double lambda){ //, const int lastkeeper){
  
  
  // G is never inverted; G^{-1} is never actually constructed.
  // Instead relevant quantities from symmetric G^{-1} = f(Q)Q' are computed
  // where Q is the eigenvectors and f(Q) elementwise divides each vector by
  // Eigenvalues + lambda a la bMultDiag() 
  
  double Le = 0; // Leave one out error loss
  
  int N = Eigenvectors.n_rows;
  int K = Eigenvectors.n_cols; // accounts for possibility of eigentruncation
                               // as handled by bEigen() and eigen.cpp

  colvec Ginv_diag(N); Ginv_diag.zeros();
  colvec coeffs(N); coeffs.zeros();
  
  // memptr() requires access to columns (cannot access rows)
  // but trans is costless, in place
  Eigenvectors = trans(Eigenvectors);
  
  for(int i = 0; i < N; ++i){
    
    colvec ginv(i+1); 
    // vector containing part that would be in lower triangle of G^{-1}

    mat temp_eigen(Eigenvectors.memptr(), K, i+1, false);
    ginv = (Eigenvectors.col(i).t()/(Eigenvalues + lambda)) * temp_eigen;
      
    Ginv_diag[i] = ginv[i];
    coeffs(span(0, i-1)) += ginv * y[i];
    coeffs[i] += sum(ginv * y(span(0, i)));
      
    // checking for user interrupt
    if(i % 501 == 0){        
      Rcpp::checkUserInterrupt();
      Rprintf(".");
    }
  }  
  Eigenvectors = trans(Eigenvectors);

  for(int i = 0; i < N; ++i){
    Le += pow((coeffs[i]/Ginv_diag[i]), 2);
  }
  
  List out(2);
  out[0] = Le;
  out[1] = coeffs;
  
  return out;
}

// [[Rcpp::export]]
List BigSolveForc(SEXP pEigenvectors, const arma::colvec Eigenvalues, 
                  const arma::colvec y, const double lambda){
  
  XPtr<SharedMemoryBigMatrix> xpEigenvectors(pEigenvectors);
  
  List out = xBigSolveForc(arma::Mat<double>((double *)xpEigenvectors->matrix(), 
                                             xpEigenvectors->nrow(), 
                                             xpEigenvectors->ncol(), false),
                                             Eigenvalues, y, lambda);
  return out;
}
