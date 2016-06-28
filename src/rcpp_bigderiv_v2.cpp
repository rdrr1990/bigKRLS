#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]

template <typename T>
void xBigDerivMat(const Mat<T>& inBigX, const Mat<T>& inSigma,
                  const Mat<T>& inKernel,const Mat<T>& inCoeffs,
                  const Mat<T> inVCovMatC, const Mat<T>& inBigXSD,
                  Mat<T> outDerivatives, Mat<T> outVarAvgDerivatives) {
  
  int N = inBigX.n_rows;
  int D = inBigX.n_cols;
  
  double sigma = as_scalar(inSigma);
  
  vec unique_vals;
  int n_unique;
  
  for(int j = 0; j < D; j++){
    unique_vals = unique(inBigX.col(j));
    n_unique = unique_vals.n_elem;
    
    // binary case
    if(n_unique == 2){
      
      // some initial setup and variable initialization for derivatives and varavgderiv
      double z0 = min(inBigX.col(j));
      double z1 = max(inBigX.col(j));
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
        int c1 = inBigX.at(i,j) == z0;
        
        vec both_max = arma::conv_to<arma::vec>::from(inBigX.col(j) + inBigX.at(i,j) == 2*z1);
        vec both_min = arma::conv_to<arma::vec>::from(inBigX.col(j) + inBigX.at(i,j) == 2*z0);
        vec first_greater = arma::conv_to<arma::vec>::from(inBigX.at(i,j) > inBigX.col(j));
        vec second_greater = arma::conv_to<arma::vec>::from(inBigX.at(i,j) < inBigX.col(j));
        
        vec adj_T_local = 0*both_max + 1*both_min - 1*first_greater;
        vec adj_C_local = 1*both_max + 0*both_min - 1*second_greater;
        
        adj_T.row(i) = adj_T_local + first_greater - second_greater;
        adj_C.row(i) = adj_C_local - first_greater + second_greater;
        
        KT_rowsums[i] = as_scalar(exp(adj_T_local*phi) * inKernel.col(i));
        KC_rowsums[i] = as_scalar(exp(adj_C_local*phi) * inKernel.col(i));

        rowvec c2 = exp((-2*(both_max + both_min) + 1) * pow((z1 - z0), 2)/sigma);
        
        arma::mat val;
        rowvec kernel_vec = inKernel.col(i);
        val = (inBigXSD[0,j] * pow(-1, c1) * (1 - c2) % kernel_vec) * inCoeffs;
        
        outDerivatives.at(i,j) = as_scalar(val);
      }
      

      
      // BIG LOOP #2
      // this finishes the varavgderivative calculation      
      double vcv_sum = sum(sum((exp(adj_T*phi) % inKernel * inVCovMatC.t()), 0) % KT_rowsums +
                           sum((exp(adj_C*phi) % inKernel * inVCovMatC.t()), 0) % KC_rowsums -
                           2*sum((exp(adj_T*phi) % inKernel * inVCovMatC.t()), 0) % KC_rowsums);
      
      double vcv_hat =  2 * pow(inBigXSD[j], 2) * vcv_sum/pow(N, 2);
      outVarAvgDerivatives[j] = vcv_hat;
      
    }
    // continuous case - much easier!
    else{
      arma::mat differences(N,N);
      arma::mat L(N,N);
      
      for(int i = 0; i< N; i++){
        differences.col(i) = inBigX.col(j) - inBigX.at(i,j);
      }
      
      L = differences % inKernel;
      outDerivatives.col(j) = (-2/sigma) * (L * inCoeffs);
      
      outVarAvgDerivatives[0,j] = (1/pow(N, 2)) * pow((-2/sigma), 2) * sum(sum(L.t() * inVCovMatC * L));
    }
  }
}

// [[Rcpp::export]]
void BigDerivMat(SEXP pinBigX, SEXP pInSigma, SEXP pInKernel,
                 SEXP pInCoeffs, SEXP pInVCovMatC, SEXP pInBigXSD,
                 SEXP pOutDerivatives, SEXP pOutVarAvgDerivatives) {
                   
  XPtr<BigMatrix> xpinBigX(pinBigX);
  XPtr<BigMatrix> xpInSigma(pInSigma);
  XPtr<BigMatrix> xpInKernel(pInKernel);
  XPtr<BigMatrix> xpInCoeffs(pInCoeffs);
  XPtr<BigMatrix> xpInVCovMatC(pInVCovMatC);
  XPtr<BigMatrix> xpInBigXSD(pInBigXSD);
  XPtr<BigMatrix> xpOutDerivatives(pOutDerivatives);
  XPtr<BigMatrix> xpOutVarAvgDerivatives(pOutVarAvgDerivatives);

  xBigDerivMat(
    arma::Mat<double>((double *)xpinBigX->matrix(), xpinBigX->nrow(), xpinBigX->ncol(), false),
    arma::Mat<double>((double *)xpInSigma->matrix(), xpInSigma->nrow(), xpInSigma->ncol(), false),
    arma::Mat<double>((double *)xpInKernel->matrix(), xpInKernel->nrow(), xpInKernel->ncol(), false),
    arma::Mat<double>((double *)xpInCoeffs->matrix(), xpInCoeffs->nrow(), xpInCoeffs->ncol(), false),
    arma::Mat<double>((double *)xpInVCovMatC->matrix(), xpInVCovMatC->nrow(), xpInVCovMatC->ncol(), false),
    arma::Mat<double>((double *)xpInBigXSD->matrix(), xpInBigXSD->nrow(), xpInBigXSD->ncol(), false),
    arma::Mat<double>((double *)xpOutDerivatives->matrix(), xpOutDerivatives->nrow(), xpOutDerivatives->ncol(), false),
    arma::Mat<double>((double *)xpOutVarAvgDerivatives->matrix(), xpOutVarAvgDerivatives->nrow(), 
                      xpOutVarAvgDerivatives->ncol(), false)
  );
}
