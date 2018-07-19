## function to simulate outputs
runRef <- function(i, priors, func) {
    
    ## set up vector
    pars <- rep(NA, nrow(priors))
    
    ## sample from prior
    for(j in 1:nrow(priors)) {
        pars[j] <- do.call(priors$dist[j], list(1, priors$p1[j], priors$p2[j]))
    }
    
    ## simulate from model
    out <- func(pars)
    
    ## return output
    list(pars = pars, out = out)
}      
