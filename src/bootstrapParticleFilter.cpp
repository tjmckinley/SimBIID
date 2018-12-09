#include "functions.hpp"

// bootstrap particle filter
double bootstrapPartFilter (int N, arma::vec pars, IntegerMatrix state, IntegerMatrix stateNew, 
                            NumericVector weights, NumericVector weightsNew, NumericMatrix dataset, SEXP func_)
{
    // N is number of particles
    // pars is vector of parameters
    // state and stateNew are integer vectors
    // weights* are vector of particle weights
    // dataset is matrix of time and then time series counts
    // func_ is simulation function
    
    // initialise variables
    int j, k, t, r;
    double totweight = 0.0, u = 0.0;
    
    // initialise log-likelihood
    double LL = 0.0;
    
    // setup output
    NumericVector out(state.ncol() + 1);
    
    // extract function pointer
    funcPtr func = *XPtr<funcPtr>(func_);
    
    // setup data
    IntegerVector counts(dataset.ncol() - 1);
    
    // setup cumulative weights
    for(k = 0; k < weights.size(); k++){
        weights[k] = (double) k;
    }
    totweight = (double) weights.size();
    
    // loop over time series
    for(t = 0; t < (dataset.nrow() - 1); t++){
        
        // set data
        for(k = 0; k < counts.size(); k++) {
            counts[k] = (int) dataset(t + 1, k + 1);
        }
        
        // loop over particles
        for(k = 0; k < N; k++) {
            // simulate forward
            if(t == 0) {
                r = k;
            } else {
                // resample a particle from the previous time step
                r = 0;
                u = R::runif(0.0, totweight);
                while(u > weights[r]){
                    r++;
                }
            }
            out = core_processing<funcPtr>(func, as<NumericVector>(wrap(pars)), 
                    dataset(t, 0), dataset(t + 1, 0), state(r, _), counts);
            for(j = 0; j < state.ncol(); j++){
                stateNew(k, j) = (int) out[j + 1];
            }
            for(j = 0; j < out.size(); j++){
                Rprintf("out(%d) = %f ", j, out[j]);
            }
            Rprintf("\n");
            Rprintf("count(0) = %d\n", counts[0]);
            // set new weight
            weightsNew[k] = out[0];
        }
        
        // update states and weights
        state = stateNew;
        weights = weightsNew;
        totweight = sum(weights);
        
        // update log-likelihood
        LL += log(totweight) - log(N);
    }
    // return log-likelihood
    return(LL);
}
