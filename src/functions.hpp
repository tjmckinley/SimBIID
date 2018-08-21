#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;

// set up code to pass custom cpp function
typedef SEXP (*funcPtr)(NumericVector pars, double tstart, double tstop, int tol,
           IntegerVector state, int counts);

// wrapper function to run custom Rcpp function           
template <typename T> 
IntegerVector core_processing(T func, NumericVector pars, double tstart, 
            double tstop, int tol, IntegerVector state, int counts) {
    return as<IntegerVector>(func(pars, tstart, tstop, tol, state, counts));
}

// function to compute Cholesky decomposition and to deal with
// cases where it fails
arma::mat cholArma(arma::mat sigma, double *scale);

// function to draw from a multivariate normal
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat cholSigma);

// function to calculate variance-covariance matrix from posterior samples
void calcPost(int i, int npars, 
                  arma::vec *tempmn, arma::mat *meanmat, 
                  arma::mat *meanmat1, NumericMatrix posterior, 
                  arma::mat *propcov);
                  
// function to update variance-covariance matrix based on current samples
void adaptUpdate(int i, int npars, 
                  arma::vec *tempmn, arma::mat *meanmat, 
                  arma::mat *meanmat1, NumericVector posterior, 
                  arma::mat *propcov);

// alive particle filter
double AlivePartFilter (int N, arma::vec pars, NumericMatrix dataset,
                        int tol, IntegerMatrix state, IntegerMatrix stateNew,
                        int *cumnt, int nmultskip, SEXP func_);

// a Metropolis-Hastings PMCMC algorithm for fitting time series models

// [[Rcpp::export]]
List PMCMC_cpp (NumericMatrix dataset, NumericMatrix priors, CharacterVector parnames, 
    NumericVector iniPars, NumericMatrix propVar_R,
    int niter, int npart, double scale, int tol, int nprintsum, int nmultskip, 
    int nupdate, int fixpars, int adapt, IntegerVector iniStates, SEXP func_);

#endif // __FUNCTIONS__
