#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigDerivMat(const Mat<T>& X, const Mat<T>& K, const Mat<T> VCovMatC, 
                  Mat<T> Derivatives, Mat<T> VarAvgDerivatives, 
                  const arma::colvec Xsd,  const arma::colvec coeffs, 
                  const double sigma) {
  
  int N = X.n_rows;
  int D = X.n_cols;
  
  vec unique_vals;
  int n_unique;
  
  for(int j = 0; j < D; j++){
    
    Rcpp::checkUserInterrupt();
    
    unique_vals = unique(X.col(j));
    n_unique = unique_vals.n_elem;
    // binary case
    if(n_unique == 2){
      
      // some initial setup and variable initialization for derivatives and varavgderiv
      double z0 = min(X.col(j));
      double z1 = max(X.col(j));
      double phi = -pow((z1-z0),2)/sigma;
        
      colvec KT_rowsums(N);
      colvec KC_rowsums(N);
      
      arma::Mat<double> adj_T(N,N);
      arma::Mat<double> adj_C(N,N);
      
      // BIG LOOP #1
      // this does 3 things:
      // - creates the derivative matrix
      // - creates adj_t and adj_c for the varavgderiv matrix
      // - creates rowsums for the treatment and control variance matrices
      for(int i = 0; i<N; i++){
        
        //constants for the derivative matrix
        int c1 = X.at(i,j) == z0;
        
        vec both_max = arma::conv_to<arma::vec>::from(X.col(j) + X.at(i,j) == 2*z1);
        vec both_min = arma::conv_to<arma::vec>::from(X.col(j) + X.at(i,j) == 2*z0);
        vec first_greater = arma::conv_to<arma::vec>::from(X.at(i,j) > X.col(j));
        vec second_greater = arma::conv_to<arma::vec>::from(X.at(i,j) < X.col(j));
        
        vec adj_T_local = 0*both_max + 1*both_min - 1*first_greater;
        vec adj_C_local = 1*both_max + 0*both_min - 1*second_greater;
        
        adj_T.row(i) = adj_T_local + first_greater - second_greater;
        adj_C.row(i) = adj_C_local - first_greater + second_greater;
        
        KT_rowsums[i] = as_scalar(exp(adj_T_local*phi) * K.col(i));
        KC_rowsums[i] = as_scalar(exp(adj_C_local*phi) * K.col(i));

        rowvec c2 = exp((-2*(both_max + both_min) + 1) * pow((z1 - z0), 2)/sigma);
        
        arma::mat val;
        rowvec kernel_vec = K.col(i);
        val = (Xsd[j] * pow(-1, c1) * (1 - c2) % kernel_vec) * coeffs;
        
        Derivatives.at(i,j) = as_scalar(val);
        
        // checking for user interrupt after each thousand observations
        if(i % 500 == 0){
          Rcpp::checkUserInterrupt();
          Rprintf(".");
        }
      }
      
      // BIG LOOP #2
      // this finishes the varavgderivative calculation      
      double vcv_sum = sum(sum((exp(adj_T*phi) % K * VCovMatC.t()), 0) % KT_rowsums +
                           sum((exp(adj_C*phi) % K * VCovMatC.t()), 0) % KC_rowsums -
                           2*sum((exp(adj_T*phi) % K * VCovMatC.t()), 0) % KC_rowsums);
      
      double vcv_hat =  2 * pow(Xsd[j], 2) * vcv_sum/pow(N, 2);
      VarAvgDerivatives[j] = vcv_hat;
    }
    
    // continuous case - much easier!
    else{
      arma::mat differences(N,N);
      arma::mat L(N,N);
      
      for(int i = 0; i< N; i++){
        differences.col(i) = X.col(j) - X.at(i,j);
        if(i % 500 == 0){
          Rcpp::checkUserInterrupt();
          Rprintf(".");
        }
      }
      
      L = differences % K;
      Derivatives.col(j) = (-2/sigma) * (L * coeffs);
      
      VarAvgDerivatives[j] = (1/pow(N, 2)) * pow((-2/sigma), 2) * sum(sum(L.t() * VCovMatC * L));
    }
  
  // checking for user interrupt after each variable
  Rcpp::checkUserInterrupt();
  }
}

// [[Rcpp::export]]
void BigDerivMat(SEXP pX, SEXP pK, SEXP pVCovMatC, SEXP pDerivatives,
                 SEXP pVarAvgDerivatives, const arma::colvec Xsd, 
                 const arma::colvec coeffs, const double sigma) {
                   
  XPtr<SharedMemoryBigMatrix> xpX(pX);
  XPtr<SharedMemoryBigMatrix> xpK(pK);
  XPtr<SharedMemoryBigMatrix> xpVCovMatC(pVCovMatC);
  XPtr<SharedMemoryBigMatrix> xpDerivatives(pDerivatives);
  XPtr<SharedMemoryBigMatrix> xpVarAvgDerivatives(pVarAvgDerivatives);

  xBigDerivMat(
    arma::Mat<double>((double *)xpX->matrix(), xpX->nrow(), xpX->ncol(), false),
    arma::Mat<double>((double *)xpK->matrix(), xpK->nrow(), xpK->ncol(), false),
    arma::Mat<double>((double *)xpVCovMatC->matrix(), xpVCovMatC->nrow(), xpVCovMatC->ncol(), false),
    arma::Mat<double>((double *)xpDerivatives->matrix(), xpDerivatives->nrow(), xpDerivatives->ncol(), false),
    arma::Mat<double>((double *)xpVarAvgDerivatives->matrix(), xpVarAvgDerivatives->nrow(), xpVarAvgDerivatives->ncol(), false),
    Xsd, coeffs, sigma
  );
}
