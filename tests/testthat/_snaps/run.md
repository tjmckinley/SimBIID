# mparse/run works

    Rcpp_object <- Rcpp::cppFunction('NumericVector simFunction(NumericVector pars, double tstart, double tstop, IntegerVector u) { 
    
        // initialise variables
        double tstar = 0.0, u_tmp = 0.0, totrate = 0.0, cumrate = 0.0;
        int i = 0;
        
        // initialise time and rates of the system
        double t = tstart;
        NumericVector rates(2);
    
        // update rates
        rates[0] = pars[0]*u[0]*u[1];
        rates[1] = pars[1]*u[1];
        totrate = sum(rates);
        
        // check states
        for(i = 0; i < u.size(); i++){
            if(u[i] < 0) {
                stop("Some states less than zero");
            }
        }
        
        // check rates
        for(i = 0; i < rates.size(); i++){
            if(R_FINITE(rates[i])) {
                if(rates[i] < 0.0) {
                    stop("Some rates less than zero or non-finite");
                }
            } else {
                stop("Some rates are non-finite");
            }
        }
    
        // set up output vector
        NumericVector out(u.size() + 2);
        
        // sample next event time
        if(totrate > 0.0) {
            tstar = t + R::rexp(1.0 / totrate);
            while(tstar < tstop){
                // update event type
                u_tmp = R::runif(0.0, totrate);
                cumrate = rates[0];
                if(u_tmp < cumrate) {
                    u[0]--;
                    u[1]++;
                } else {
                    u[1]--;
                    u[2]++;
                }
            
                // update rates
                rates[0] = pars[0]*u[0]*u[1];
                rates[1] = pars[1]*u[1];
                totrate = sum(rates);
                
                // check states
                for(i = 0; i < u.size(); i++){
                    if(u[i] < 0) {
                        stop("Some states less than zero");
                    }
                }
                
                // check rates
                for(i = 0; i < rates.size(); i++){
                    if(R_FINITE(rates[i])) {
                        if(rates[i] < 0.0) {
                            stop("Some rates less than zero or non-finite");
                        }
                    } else {
                        stop("Some rates are non-finite");
                    }
                }
                
                // update time
                t = tstar;
                
                // sample next event time
                if(totrate > 0.0) {
                    tstar = t + R::rexp(1.0 / totrate);
                } else {
                    tstar = tstop;
                }
            }
        }
        
        // record final event time
        if(totrate > 0.0) {
            t = tstop;
        }
        
        // set output
        out[0] = (totrate == 0.0 ? 1:0);
        out[1] = t;
        out[Range(2, u.size() + 1)] = as<NumericVector>(u);
        
        // return output
        return out;
    }')

---

    'SimBIID_runs' object with n = 1 replicates.
    
    Output at final time point:
    # A tibble: 1 x 5
      completed     t     S     I     R
          <dbl> <dbl> <dbl> <dbl> <dbl>
    1         0    20   117     1     2

---

    Rcpp_object <- Rcpp::cppFunction('NumericVector simFunction(NumericVector pars, double tstart, double tstop, IntegerVector u) { 
    
        // initialise variables
        double tstar = 0.0, u_tmp = 0.0, totrate = 0.0, cumrate = 0.0;
        int i = 0;
        
        // initialise time and rates of the system
        double t = tstart;
        NumericVector rates(2);
    
        // update rates
        rates[0] = pars[0]*u[0]*u[1];
        rates[1] = pars[1]*u[1];
        totrate = sum(rates);
        
        // check states
        for(i = 0; i < u.size(); i++){
            if(u[i] < 0) {
                stop("Some states less than zero");
            }
        }
        
        // check rates
        for(i = 0; i < rates.size(); i++){
            if(R_FINITE(rates[i])) {
                if(rates[i] < 0.0) {
                    stop("Some rates less than zero or non-finite");
                }
            } else {
                stop("Some rates are non-finite");
            }
        }
    
        // set up output vector
        NumericVector out(u.size() + 2);
        
        // sample next event time
        if(totrate > 0.0) {
            tstar = t + R::rexp(1.0 / totrate);
            while(tstar < tstop){
                // update event type
                u_tmp = R::runif(0.0, totrate);
                cumrate = rates[0];
                if(u_tmp < cumrate) {
                    u[0]--;
                    u[1]++;
                    u[4]++;
                } else {
                    u[1]--;
                    u[2]++;
                    u[5]++;
                }
            
                // update rates
                rates[0] = pars[0]*u[0]*u[1];
                rates[1] = pars[1]*u[1];
                totrate = sum(rates);
                
                // check states
                for(i = 0; i < u.size(); i++){
                    if(u[i] < 0) {
                        stop("Some states less than zero");
                    }
                }
                
                // check rates
                for(i = 0; i < rates.size(); i++){
                    if(R_FINITE(rates[i])) {
                        if(rates[i] < 0.0) {
                            stop("Some rates less than zero or non-finite");
                        }
                    } else {
                        stop("Some rates are non-finite");
                    }
                }
                
                // update time
                t = tstar;
                
                // sample next event time
                if(totrate > 0.0) {
                    tstar = t + R::rexp(1.0 / totrate);
                } else {
                    tstar = tstop;
                }
            }
        }
        
        // record final event time
        if(totrate > 0.0) {
            t = tstop;
        }
        
        // set output
        out[0] = (totrate == 0.0 ? 1:0);
        out[1] = t;
        out[Range(2, u.size() + 1)] = as<NumericVector>(u);
        
        // return output
        return out;
    }')

---

    'SimBIID_runs' object with n = 1 replicates.
    
    Output at final time point:
    # A tibble: 1 x 8
      completed     t     S     I     R S_inc I_inc R_inc
          <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    1         0    20   117     1     2     0     3     2

---

    Rcpp_object <- Rcpp::cppFunction('List simFunction(NumericVector pars, double tstart, double tstop, IntegerVector u, NumericVector tspan) { 
    
        // initialise variables
        double tstar = 0.0, u_tmp = 0.0, totrate = 0.0, cumrate = 0.0;
        int i = 0;
        
        // initialise time and rates of the system
        double t = tstart;
        NumericVector rates(2);
    
        // update rates
        rates[0] = pars[0]*u[0]*u[1];
        rates[1] = pars[1]*u[1];
        totrate = sum(rates);
        
        // check states
        for(i = 0; i < u.size(); i++){
            if(u[i] < 0) {
                stop("Some states less than zero");
            }
        }
        
        // check rates
        for(i = 0; i < rates.size(); i++){
            if(R_FINITE(rates[i])) {
                if(rates[i] < 0.0) {
                    stop("Some rates less than zero or non-finite");
                }
            } else {
                stop("Some rates are non-finite");
            }
        }
    
        // set up output vector
        NumericMatrix out(tspan.size() + 1, u.size() + 2);
    
        // check tspan
        for(i = 1; i < tspan.size(); i++){
            if(tspan[i] <= tspan[i - 1]) {
                stop("tspan not ordered");
            }
            if(tstart >= tspan[0]){
                stop("tstart >= tspan[0]");
            }
            if(tstop < tspan[tspan.size() - 1]){
                stop("tstop < tspan[n]");
            }
        }
    
        // update tspan
        int k = 0;
        while(tspan[k] < t && k < tspan.size()) {
            out(k, 0) = NA_REAL;
            out(k, 1) = tspan[k];
            for(i = 0; i < u.size(); i++) {
                out(k, i + 2) = (double) u[i];
            }
            k++;
        }
        
        // sample next event time
        if(totrate > 0.0) {
            tstar = t + R::rexp(1.0 / totrate);
            while(tstar < tstop){
    
                // update tspan
                while(tspan[k] < tstar && k < tspan.size()) {
                    out(k, 0) = NA_REAL;
                    out(k, 1) = tspan[k];
                    for(i = 0; i < u.size(); i++) {
                        out(k, i + 2) = (double) u[i];
                    }
                    k++;
                }
                // update event type
                u_tmp = R::runif(0.0, totrate);
                cumrate = rates[0];
                if(u_tmp < cumrate) {
                    u[0]--;
                    u[1]++;
                } else {
                    u[1]--;
                    u[2]++;
                }
            
                // update rates
                rates[0] = pars[0]*u[0]*u[1];
                rates[1] = pars[1]*u[1];
                totrate = sum(rates);
                
                // check states
                for(i = 0; i < u.size(); i++){
                    if(u[i] < 0) {
                        stop("Some states less than zero");
                    }
                }
                
                // check rates
                for(i = 0; i < rates.size(); i++){
                    if(R_FINITE(rates[i])) {
                        if(rates[i] < 0.0) {
                            stop("Some rates less than zero or non-finite");
                        }
                    } else {
                        stop("Some rates are non-finite");
                    }
                }
                
                // update time
                t = tstar;
                
                // sample next event time
                if(totrate > 0.0) {
                    tstar = t + R::rexp(1.0 / totrate);
                } else {
                    tstar = tstop;
                }
            }
        }
        
        // record final event time
        if(totrate > 0.0) {
            t = tstop;
        }
        
        // set output
        while(k < tspan.size()) {
            out(k, 0) = NA_REAL;
            out(k, 1) = tspan[k];
            for(i = 0; i < u.size(); i++) {
                out(k, i + 2) = (double) u[i];
            }
            k++;
        }
        out(tspan.size(), 0) = (totrate == 0.0 ? 1:0);
        out(tspan.size(), 1) = t;
        for(i = 0; i < u.size(); i++) {
            out(tspan.size(), i + 2) = (double) u[i];
        }
        
        // return output
        List out1(2);
        out1[0] = out(tspan.size(), _);
        out1[1] = out(Range(0, tspan.size() - 1), Range(1, out.ncol() - 1));
        return out1;
    }')

---

    'SimBIID_runs' object with n = 1 replicates.
    
    Output at final time point:
    # A tibble: 1 x 5
      completed     t     S     I     R
          <dbl> <dbl> <dbl> <dbl> <dbl>
    1         0    20   117     1     2
    
    Time-series counts:
    # A tibble: 20 x 4
           t     S     I     R
       <dbl> <dbl> <dbl> <dbl>
     1     1   119     1     0
     2     2   118     2     0
     3     3   118     2     0
     4     4   118     2     0
     5     5   118     2     0
     6     6   118     2     0
     7     7   117     2     1
     8     8   117     1     2
     9     9   117     1     2
    10    10   117     1     2
    11    11   117     1     2
    12    12   117     1     2
    13    13   117     1     2
    14    14   117     1     2
    15    15   117     1     2
    16    16   117     1     2
    17    17   117     1     2
    18    18   117     1     2
    19    19   117     1     2
    20    20   117     1     2

---

    Rcpp_object <- Rcpp::cppFunction('List simFunction(NumericVector pars, double tstart, double tstop, IntegerVector u, NumericVector tspan) { 
    
        // initialise variables
        double tstar = 0.0, u_tmp = 0.0, totrate = 0.0, cumrate = 0.0;
        int i = 0;
        
        // initialise time and rates of the system
        double t = tstart;
        NumericVector rates(2);
    
        // update rates
        rates[0] = pars[0]*u[0]*u[1];
        rates[1] = pars[1]*u[1];
        totrate = sum(rates);
        
        // check states
        for(i = 0; i < u.size(); i++){
            if(u[i] < 0) {
                stop("Some states less than zero");
            }
        }
        
        // check rates
        for(i = 0; i < rates.size(); i++){
            if(R_FINITE(rates[i])) {
                if(rates[i] < 0.0) {
                    stop("Some rates less than zero or non-finite");
                }
            } else {
                stop("Some rates are non-finite");
            }
        }
    
        // set up output vector
        NumericMatrix out(tspan.size() + 1, u.size() + 2);
    
        // check tspan
        for(i = 1; i < tspan.size(); i++){
            if(tspan[i] <= tspan[i - 1]) {
                stop("tspan not ordered");
            }
            if(tstart >= tspan[0]){
                stop("tstart >= tspan[0]");
            }
            if(tstop < tspan[tspan.size() - 1]){
                stop("tstop < tspan[n]");
            }
        }
    
        // update tspan
        int k = 0;
        while(tspan[k] < t && k < tspan.size()) {
            out(k, 0) = NA_REAL;
            out(k, 1) = tspan[k];
            for(i = 0; i < u.size(); i++) {
                out(k, i + 2) = (double) u[i];
            }
            for(i = u.size() / 2; i < u.size(); i++) {
                u[i] = 0;
            }
            k++;
        }
        
        // sample next event time
        if(totrate > 0.0) {
            tstar = t + R::rexp(1.0 / totrate);
            while(tstar < tstop){
    
                // update tspan
                while(tspan[k] < tstar && k < tspan.size()) {
                    out(k, 0) = NA_REAL;
                    out(k, 1) = tspan[k];
                    for(i = 0; i < u.size(); i++) {
                        out(k, i + 2) = (double) u[i];
                    }
                    for(i = u.size() / 2; i < u.size(); i++) {
                        u[i] = 0;
                    }
                    k++;
                }
                // update event type
                u_tmp = R::runif(0.0, totrate);
                cumrate = rates[0];
                if(u_tmp < cumrate) {
                    u[0]--;
                    u[1]++;
                    u[4]++;
                } else {
                    u[1]--;
                    u[2]++;
                    u[5]++;
                }
            
                // update rates
                rates[0] = pars[0]*u[0]*u[1];
                rates[1] = pars[1]*u[1];
                totrate = sum(rates);
                
                // check states
                for(i = 0; i < u.size(); i++){
                    if(u[i] < 0) {
                        stop("Some states less than zero");
                    }
                }
                
                // check rates
                for(i = 0; i < rates.size(); i++){
                    if(R_FINITE(rates[i])) {
                        if(rates[i] < 0.0) {
                            stop("Some rates less than zero or non-finite");
                        }
                    } else {
                        stop("Some rates are non-finite");
                    }
                }
                
                // update time
                t = tstar;
                
                // sample next event time
                if(totrate > 0.0) {
                    tstar = t + R::rexp(1.0 / totrate);
                } else {
                    tstar = tstop;
                }
            }
        }
        
        // record final event time
        if(totrate > 0.0) {
            t = tstop;
        }
        
        // set output
        while(k < tspan.size()) {
            out(k, 0) = NA_REAL;
            out(k, 1) = tspan[k];
            for(i = 0; i < u.size(); i++) {
                out(k, i + 2) = (double) u[i];
            }
            for(i = u.size() / 2; i < u.size(); i++) {
                u[i] = 0;
            }
            k++;
        }
        out(tspan.size(), 0) = (totrate == 0.0 ? 1:0);
        out(tspan.size(), 1) = t;
        for(i = 0; i < u.size(); i++) {
            out(tspan.size(), i + 2) = (double) u[i];
        }
        
        // return output
        List out1(2);
        out1[0] = out(tspan.size(), _);
        out1[1] = out(Range(0, tspan.size() - 1), Range(1, out.ncol() - 1));
        return out1;
    }')

---

    'SimBIID_runs' object with n = 1 replicates.
    
    Output at final time point:
    # A tibble: 1 x 8
      completed     t     S     I     R S_inc I_inc R_inc
          <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    1         0    20   117     1     2     0     0     0
    
    Time-series counts:
    # A tibble: 20 x 7
           t     S     I     R S_inc I_inc R_inc
       <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
     1     1   119     1     0     0     1     0
     2     2   118     2     0     0     1     0
     3     3   118     2     0     0     0     0
     4     4   118     2     0     0     0     0
     5     5   118     2     0     0     0     0
     6     6   118     2     0     0     0     0
     7     7   117     2     1     0     1     1
     8     8   117     1     2     0     0     1
     9     9   117     1     2     0     0     0
    10    10   117     1     2     0     0     0
    11    11   117     1     2     0     0     0
    12    12   117     1     2     0     0     0
    13    13   117     1     2     0     0     0
    14    14   117     1     2     0     0     0
    15    15   117     1     2     0     0     0
    16    16   117     1     2     0     0     0
    17    17   117     1     2     0     0     0
    18    18   117     1     2     0     0     0
    19    19   117     1     2     0     0     0
    20    20   117     1     2     0     0     0

