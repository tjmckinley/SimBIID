#include "functions.hpp"

// a Metropolis-Hastings PMCMC algorithm for fitting time series models
List PMCMC_cpp (NumericMatrix dataset, NumericMatrix priors, CharacterVector parnames, 
    NumericVector iniPars, NumericMatrix propVar_R,
    int niter, int npart, double scale, IntegerVector tols, IntegerVector whichind, 
    int nprintsum, int nmultskip, 
    int nupdate, int fixpars, int adapt, IntegerVector iniState, SEXP func_)
{
    // 'dataset_R' is vector of removal times
    // 'iniPars' is a vector of initial values for the unknown parameters
    // 'priors' is an (npars x 2) matrix containing the shape and scale
    //    hyperparameters for the prior distributions
    // 'parnames' is vector of parameter names
    // 'propVar_R' is initial covariance matrix (on log or logistic scale) for parameters
    // 'niter' is the number of iterations over which to run the chain
    // 'npart' is number of particles to use for particle filter
    // 'scale' is mixing proportion for adaptive MCMC
    // 'tols' defines how close points have to match
    // 'whichind' matches states to dataset
    // 'nprintsum' determines how often summaries are printed to the screen
    // 'nmultskip' determines when to skip out of simulations
    // 'nupdate' determines when to start adaptive proposal
    // 'fixpars' is indicator determining whether to fix the parameters
    //      (useful for optimising variance)
    // 'adapt' determines whether to use adaptive proposal
    // 'iniState' is a vector of initial states
    // 'func_' is pointer to simulation function
    
    // initialise variables
    int i, j, k, l;
    double u;
    int cumnt, cumntProp;
    
    // calculate number of parameters
    int npars = iniPars.size();
    
    // print arguments to the screen
    Rprintf("Number of iterations: %d\n", niter);
    Rprintf("Number of particles: %d\n", npart);
    if(fixpars == 0) {
        Rprintf("Mixing proportion for proposal: %.2f\n", scale);
        if(adapt == 1) {
            Rprintf("Start adaptive proposal at iteration: %d\n", nupdate);
        }
    }
    Rprintf("Initial upper bound on simulations before rejecting: %d x %d = %d\n", nmultskip, npart, nmultskip * npart);
    // rescale since skipping in two parts
    // (integer division OK)
    int nmultskipini = nmultskip;
    nmultskip *= npart;
    
    Rprintf("\nNumber of parameters: %d\n\n", npars);
    Rprintf("Priors:\n");
    for(i = 0; i < priors.nrow(); i++) {
        if(priors(i, 0) == 1) {
            Rcpp::Rcout << parnames(i) << " ~ U(";
        }
        else {
            if(priors(i, 0) == 2) {
            Rcpp::Rcout << parnames(i) << " ~ N(";
            } else {
            Rcpp::Rcout << parnames(i) << " ~ G(";
            }
        }
        Rprintf("%.2f, %.2f)\n", priors(i, 1), priors(i, 2));
    }
    Rprintf("\n");
    
    Rprintf("Fixed tolerances for ABC:\n");
    for(i = 0; i < tols.size(); i++) {
        Rprintf("state[%d] tol = %d\n", whichind[i], tols[i]);
    }
    Rprintf("\n");
    int nclass = iniState.length();
    Rprintf("Number of classes: %d\n\n", nclass);
    Rprintf("Initial states:\n");
    for(i = 0; i < nclass; i++) {
        Rprintf("state[%d] = %d\n", i, iniState(i));
    }
    Rprintf("\n");
    
    // set up output matrix of length 'niter' to record chain
    // (append extra column for unnormalised posterior and
    // another two for no. of simulations)
    j = (fixpars == 0 ? (npars + 3):1);
    NumericMatrix output(niter, j);
    
    // set up temporary vectors for particle filtering
    IntegerMatrix state(npart + 1, nclass);
    IntegerMatrix stateNew(npart + 1, nclass);
    
    // monitor skiprate
    int nskip = 0, cumskip = 0;
    double skiprate = 0.0;
    
    // set variables for calculating log-likelihoods
    double LL, accCurr, accProp, acc;
    
    // initialise chain and set up vector to hold proposals
    arma::vec pars(npars), parsProp(npars);
    for(i = 0; i < npars; i++) {
        if(NumericVector::is_na(iniPars[i])){
            if(priors(i, 0) == 1) {
                pars(i) = R::runif(priors(i, 1), priors(i, 2));
            } else {
                if(priors(i, 0) == 2) {
                    pars(i) = R::rnorm(priors(i, 1), priors(i, 2));
                } else {
                    pars(i) = R::rgamma(priors(i, 1), priors(i, 2));
                }
            }   
        } else {
            pars(i) = iniPars[i];
        }
        if(priors(i, 0) == 1) {
            if(pars(i) < priors(i, 1) || pars(i) > priors(i, 2)) stop("Some initial values are not bounded correctly");
        } else {
            if(priors(i, 0) == 3) {
                if(pars(i) < 0.0) stop("Some initial values are not bounded correctly");
            }
        }
    }
    parsProp.zeros();
    
    // runs code multiple times for fixed parameters in
    // order to estimate variance of log-importance weight
    if(fixpars == 1) {
        for(k = 0; k < niter; k++) {
            if((k + 1) % 100 == 0) {
                Rprintf("k = %d\n", k + 1);
            }
            // check for user interruption
            R_CheckUserInterrupt();
            
            // set initial states
            for(j = 0; j < (npart + 1); j++) {
                for(l = 0; l < nclass; l++) {
                    state(j, l) = iniState(l);
                }
            }
            LL = AlivePartFilter(npart, pars, state, stateNew, tols, dataset, whichind, 
                                 &cumnt, nmultskip, func_);
            output(k, 0) = LL;
        }
        List outlist (2);
        outlist[0] = output;
        outlist[1] = pars;
        return(outlist);
    }
    
    // set up adaptive proposal distribution
    double cholScale = 0.0;
    arma::mat propVar(npars, npars);
    arma::mat propVarIni(npars, npars);
    arma::mat propVarChol(npars, npars);
    arma::mat propVarIniChol(npars, npars);
    propVar.zeros(); propVarIni.zeros();
    for(j = 0; j < npars; j++) {
        propVarIni(j, j) = pow(0.1, 2.0) / ((double) npars);
        for(k = 0; k < npars; k++) {
            propVar(j, k) = propVar_R(j, k);
        }
    }
    // calculate Cholesky decomposition for MVN sampling
    double adaptscale = pow(2.562, 2.0) / ((double) npars);
    propVarIniChol = cholArma(propVarIni, &cholScale);
    propVarChol = cholArma(propVar * adaptscale, &cholScale);
    arma::vec tempmn(npars);
    arma::mat meanmat(npars, npars);
    arma::mat meanmat1(npars, npars);
    
    // simulate initial conditions if required
    Rprintf("Initialising system...\n");
    int kmax = 1000;
    k = 0;
    while(k < kmax) {
        if( (k + 1) % 10 == 0) {
            Rprintf("k = %d\n", k + 1);
        }
        // check for user interruption
        R_CheckUserInterrupt();
        
        // run particle filter
        cumnt = 0;
        // set initial states
        for(j = 0; j < (npart + 1); j++) {
            for(l = 0; l < nclass; l++) {
                state(j, l) = iniState(l);
            }
        }
        LL = AlivePartFilter(npart, pars, state, stateNew, tols, dataset, whichind, 
                             &cumnt, nmultskip, func_);
        if(R_finite(LL) == 0) {
            // if initial values are provided, then reject
            if(all(!is_na(iniPars))) {
                Rprintf("Initial values do not produce valid estimate of the likelihood.\n");
                k = kmax + 2;
            } else {
                for(i = 0; i < npars; i++) {
                    if(NumericVector::is_na(iniPars[i])){
                        if(priors(i, 0) == 1) {
                            pars(i) = R::runif(priors(i, 1), priors(i, 2));
                        } else {
                            if(priors(i, 0) == 2) {
                                pars(i) = R::rnorm(priors(i, 1), priors(i, 2));
                            } else {
                                pars(i) = R::rgamma(priors(i, 1), priors(i, 2));
                            }
                        }   
                    } else {
                        pars(i) = iniPars[i];
                    }
                    if(priors(i, 0) == 1) {
                        if(pars(i) < priors(i, 1) || pars(i) > priors(i, 2)) stop("Some initial values are not bounded correctly");
                    } else {
                        if(priors(i, 0) == 3) {
                            if(pars(i) < 0.0) stop("Some initial values are not bounded correctly");
                        }
                    }
                }
            }
            //increment counter
            k++;
        } else {
            k = kmax + 1;
        }
    }
    if(R_finite(LL) == 0) {
        List outlist (1);
        outlist[0] = NA_REAL;
        if(k == (kmax + 1)) {
            Rprintf("\nCannot sample valid estimate of the likelihood.\n");
        }
        return(outlist);
    }
    Rprintf("Initialisation complete!\n\n");
    
    Rprintf("Initial parameter values:\n\n");
    for(i = 0; i < npars; i++) {
        Rcpp::Rcout << parnames(i) << " = " << pars(i);
        Rprintf("\n");
    }
    Rprintf("\n");
    
    // calculate log-likelihood – log-prior
    accCurr = LL;
    for(j = 0; j < npars; j++) {
        if(priors(j, 0) == 1) {
            accCurr += R::dunif(pars(j), priors(j, 1), priors(j, 2), 1);
        } else {
            if(priors(j, 0) == 2) {
                accCurr += R::dnorm(pars(j), priors(j, 1), priors(j, 2), 1);
            } else {
                accCurr += R::dgamma(pars(j), priors(j, 1), priors(j, 2), 1);
            }
        }  
    }
    
    //initialise timer
    Timer timer;
    int timer_cnt = 0;
    double prev_time = 0.0;
    
    // set up acceptance counters
    int nacc = 0, cumacc = 0;
    double accrate = 0.0;
    
    // start runs
    Rprintf("Starting runs...\n");
    for(i = 0; i < niter; i++) {
        //check for user interruptions
        if (i % 10 == 0) {
            R_CheckUserInterrupt();
        }
        
        // propose new value
        if(R::runif(0.0, 1.0) < scale) {
            parsProp = arma::conv_to<arma::vec>::from(mvrnormArma(1, pars, propVarIniChol));
        } else {
            parsProp = arma::conv_to<arma::vec>::from(mvrnormArma(1, pars, propVarChol));
        }
        // check validity
        k = 0;
        for(j = 0; j < npars; j++) {
            if(priors(j, 0) == 1) {
                k += (parsProp(j) <= priors(j, 1) || parsProp(j) >= priors(j, 2) ? 1:0);
            } else {
                if(priors(j, 0) == 3) {
                    k += (parsProp(j) <= 0.0 ? 1:0);
                }
            }  
        }
        // if valid then proceed with simulation
        if(k == 0) {
            // set prior information
            accProp = 0.0;
            for(j = 0; j < npars; j++) { 
                if(priors(j, 0) == 1) {
                    accProp += R::dunif(parsProp(j), priors(j, 1), priors(j, 2), 1);
                } else {
                    if(priors(j, 0) == 2) {
                        accProp += R::dnorm(parsProp(j), priors(j, 1), priors(j, 2), 1);
                    } else {
                        accProp += R::dgamma(parsProp(j), priors(j, 1), priors(j, 2), 1);
                    }
                }  
            }
            // set initial states
            for(j = 0; j < (npart + 1); j++) {
                for(l = 0; l < nclass; l++) {
                    state(j, l) = iniState(l);
                }
            }
            // run particle filter
            LL = AlivePartFilter(npart, parsProp, state, stateNew, tols, dataset, whichind, 
                                 &cumntProp, nmultskip, func_);
            
            //accept-reject proposal
            if(R_finite(LL) != 0) {
                accProp += LL;
                u = R::runif(0.0, 1.0);
                acc = accProp - accCurr;
                if(log(u) < acc) {
                    pars = parsProp;
                    cumnt = cumntProp;
                    accCurr = accProp;
                    nacc++;
                    cumacc++;
                } 
            } else {
                nskip++;
                cumskip++;
            }
        }
        
        // record number of simulations required in this iteration
        output(i, npars + 1) = cumnt;
        output(i, npars + 2) = cumntProp;
        
        // save current value of chain into output vector
        for(j = 0; j < npars; j++) {
            output(i, j) = pars(j);
        }
        output(i, npars) = accCurr;
        
        // print some output to the screen for book-keeping
        if ((i + 1) % nprintsum == 0) {	
            //calculate block run time
            timer.step("");
            NumericVector res(timer);
            
            //update skiprate
            skiprate = ((double) nskip) / ((double) nprintsum);
            nskip = 0;
            
            // update acceptance rate
            accrate = ((double) nacc) / ((double) nprintsum);
            nacc = 0;
            
            Rprintf("i = %d acc = %.2f skip = %.2f nsimskip = %d time = %.2f secs \n", i + 1, accrate, skiprate, nmultskip, (res[timer_cnt] / 1e9) - prev_time);
            
            //reset timer and acceptance rate counter
            prev_time = res[timer_cnt] / 1e9;
            timer_cnt++;
        }
        
        // calculations for adaptive proposal
        if((i + 1) >= nupdate && adapt == 1) {
            if((i + 1) == nupdate) {
                // check for acceptances
                accrate = ((double) cumacc) / ((double) nupdate);
                if(accrate > 0.0) {
                	calcPost(i, npars, &tempmn, &meanmat, &meanmat1, output, &propVar);
                } else {
                    List outlist (1);
                    outlist[0] = NA_REAL;
                    Rprintf("No initial acceptances.\n");
                    return(outlist);
                }
            } else {
            	adaptUpdate(i, npars, &tempmn, &meanmat, &meanmat1, output(i, _), &propVar);
            }
            propVarChol = cholArma(propVar * adaptscale, &cholScale);
        }
    }
    timer.step("");
    NumericVector res(timer);
    Rprintf("Final time = %.2f secs Final skip = %.2f\n", res[timer_cnt] / 1e9, cumskip / ((double) niter));
    
    // return MCMC chain
    List outlist (7);
    outlist[0] = output;
    outlist[1] = cumskip / ((double) niter);
    outlist[2] = ((double) cumacc) / ((double) niter);
    outlist[3] = nmultskipini;
    outlist[4] = npart;
    outlist[5] = res[timer_cnt] / 1e9;
    outlist[6] = propVar;
    return(outlist);
}

