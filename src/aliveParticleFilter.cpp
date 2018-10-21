#include "functions.hpp"

// alive particle filter
double AlivePartFilter (int N, arma::vec pars, IntegerMatrix state, IntegerMatrix stateNew,
                        IntegerVector tols, NumericMatrix dataset, IntegerVector whichind,
                        int *cumnt, int nmultskip, SEXP func_)
{
    // N is number of particles
    // pars is vector of parameters
    // state and stateNew are integer vectors
    // tols defines how close points have to match
    // dataset is vector of counts to match to
    // whichind is vector of indexes matching states to dataset
    // cumnt is cumulative number of simulations (to set upper bound for skipping)
    // nmultskip skips simulations
    // func_ is simulation function
    
    // initialise variables
    int j, k, t, r, matched;
    
    // initialise log-likelihood
    double logp = dataset.size() * log(N);
    
    // setup counter
    (*cumnt) = 0;
    int nts = 0;
    
    // setup output
    double LL = 0.0;
    IntegerVector out(state.ncol() + 1);
    
    // extract function pointer
    funcPtr func = *XPtr<funcPtr>(func_);
    
    // setup data
    IntegerVector counts(dataset.ncol() - 1);
    
    // loop over time series
    for(t = 0; t < (dataset.nrow() - 1); t++){
        // set counter
        nts = 0;
        
        // set data
        for(k = 0; k < counts.size(); k++) {
            counts[k] = dataset(t + 1, k + 1);
        }
        
        // loop over particles
        k = 0;
        while(k < (N + 1)){
            matched = 0;
            while(matched == 0){
                //check whether to skip out or not
                if(nmultskip > 0) {
                    if((*cumnt) > nmultskip) {
                        LL = NA_REAL;
                        return(LL);
                    }
                }
                // simulate forward
                if(t == 0) {
                    r = k;
                } else {
                    // resample a particle from the previous time step
                    r = (int) floor(R::runif(0.0, 1.0) * N);
                }
                out = core_processing<funcPtr>(func, as<NumericVector>(wrap(pars)), 
                        dataset(t, 0), dataset(t + 1, 0), state(r, _), tols,
                        counts, whichind);
                matched = out[0];
                stateNew(k, _) = out[Range(1, state.ncol())];
                // update counter
                nts++;
                (*cumnt)++;
            }
            k++;
        }
        // update states
        state = stateNew;
        
        // update log-likelihood
        logp -= log(nts - 1.0);
    }
    // return log-likelihood
    LL = logp;
    return(LL);
}
