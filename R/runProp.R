## function to generate valid proposal
runProp <- function(i, t, priors, prevWeights, prevPars, propCov, tols, data, func) {
    accrate <- 0
    valid <- 0
    while(valid == 0) {
        if(t == 1) {
            pars <- rep(NA, nrow(priors))
            ## sample from prior
            for(j in 1:nrow(priors)) {
                pars[j] <- runif(1, priors[, 1], priors[, 2])
            }
        } else {
            ## sample from previous generation
            k <- sample(1:length(prevWeights), 1, prob = prevWeights)
            pars <- rmvnorm(1, 
                            mean = prevPars[k, ], 
                            sigma = propCov
            )
            pars <- as.vector(pars)
        }
        
        if(all(apply(cbind(pars, priors), 1, function(x) {
            x[1] > x[2] & x[1] < x[3]
        }))) {
            ## simulate from model
            out <- func(pars)
            
            ## check matching
            if(all(!is.na(out))) {
                if(all(abs(out - data) < tols)) {
                    valid <- 1
                }
            }
        }
        
        ## update counter
        accrate <- accrate + 1
    }
    
    if(t == 1) {
        weightsNew <- 1
    } else {
        ## calculate unnormalised weight
        weightsNew <- prod(apply(cbind(pars, priors), 1, function(x) {
            dunif(x[1], x[2], x[3])
        }))
        weightsNew <- weightsNew / sum(prevWeights * apply(prevPars, 1, function(x, pars, propCov) {
            dmvnorm(pars, mean = x, sigma = propCov)
        }, pars = pars, propCov = propCov))
    }
    list(pars = pars, out = out, weightsNew = weightsNew, accrate = accrate)
}      
