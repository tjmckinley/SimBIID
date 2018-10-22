

    // initialise variables
    double tstar = 0.0, u_tmp = 0.0, totrate = 0.0, cumrate = 0.0;
    int i, j;
    
    // initialise time and rates of the system
    double t = tstart;
    RATELINES1 

    // set up output vector
    MATCHCRIT0
    
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
    MATCHCRIT1
    return out;
}')
