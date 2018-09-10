Rcpp_ptr <- RcppXPtrUtils::cppXPtr('SEXP simFunction(NumericVector gdata, double tstart, double tstop, int tol, IntegerVector u, int counts) {  
    // initialise variables
    double tstar = 0.0, u_tmp = 0.0, totrate = 0.0;
    int i, j;
    
    // initialise time and upRates of the system
    double t = tstart;
    IntegerVector uNew = u;
    
    totrate = sum(rates); 
    
    // sample next event time
    if(totrate > 0) {
        tstar = t + R::rexp(1.0 / totrate);
        while(tstar < tstop){
            // sample event type
            u_tmp = R::runif(0.0, totrate);
            
            totrate = sum(rates); 
            
            // update time
            t = tstar;
            
            // sample next event time
            if(totrate > 0) {
                tstar = t + R::rexp(1.0 / totrate);
            } else {
                tstar = tstop;
            }
        }
    }
    IntegerVector out(u.size() + 1);
    // check whether simulation matches data
    if(fabs(uNew[2] - counts) <= tol) {
        out[0] = 1;
    } else {
        out[0] = 0;
    }
    out[Range(1, u.size())] = uNew;
    return out;
}')
