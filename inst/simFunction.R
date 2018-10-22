Rcpp_ptr <- RcppXPtrUtils::cppXPtr('SEXP simFunction(NumericVector gdata, double tstart, double tstop, IntegerVector u, IntegerVector tols, IntegerVector counts, IntegerVector whichind) { 

    // initialise variables
    double tstar = 0.0, u_tmp = 0.0, totrate = 0.0;
    int i, j;

    RATELINES0
    
    // initialise time and rates of the system
    double t = tstart;
    IntegerVector uNew = u;
    RATELINES1 

    // set up output vector
    IntegerVector out(u.size() + 1);
    
    // sample next event time
    if(totrate > 0) {
        tstar = t + R::rexp(1.0 / totrate);
        while(tstar < tstop){
            // sample event type
            u_tmp = R::runif(0.0, totrate);
            RATELINES2 
            
            // update time
            t = tstar;
            
            // sample next event time
            if(totrate > 0) {
                tstar = t + R::rexp(1.0 / totrate);
            } else {
                tstar = tstop;
            }
            RATELINES3
        }
    }
    // check whether simulation matches data
    out[0] = 1;
    for(j = 0; j < counts.size(); j++) {
        if(fabs(uNew[whichind[j]] - counts[j]) <= tols[j]) {
            out[0] *= 1;
        } else {
            out[0] *= 0;
        }
    }
    out[Range(1, u.size())] = uNew;
    return out;
}')
