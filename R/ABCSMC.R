#' @title Runs ABC-SMC algorithm
#'
#' @description Runs ABC-SMC algorithm of Toni et al. (2009)
#'
#' @export
#'
#' @param npart     The number of particles (must be a positive integer).
#' @param tols 		A \code{matrix} of tolerances, with the number of rows defining
#'                  the number of generations of required, and the number of columns
#'                  defining the number of summary statistics / data points to match to.
#' @param priors    A \code{matrix} containing the lower and upper bounds of uniform
#'                  priors; with number of rows equal to the number of parameters, and
#'                  two columns containing the lower and upper bounds respectively.
#' @param func      A function taking a single argument \code{pars} that runs the simulator
#'                  and returns the simulated summary measures against which to compare.
#' @param data      A \code{vector} containing the observed summary statistics to match to.
#'
#' @return An \code{ABCSMC} object, essentially a \code{list} containing:
#' \itemize{
#' \item{\code{pars}:}{ a \code{list} of \code{matrix} objects containing the accepted
#'                      particles. Each element of the list corresponds to a generation 
#'                      of ABC-SMC, with each matrix being of dimension 
#'                      \code{npart} x \code{npars};}
#' \item{\code{output}:}{ a \code{list} of \code{matrix} objects containing the simulated
#'                      summary statistics. Each element of the list corresponds to a
#'                      generation of ABC-SMC, with each matrix being of dimension 
#'                      \code{npart} x \code{length(data)};}
#' \item{\code{weights}:}{ a \code{list} of \code{vector} objects containing the particle
#'                      weights. Each element of the list corresponds to a
#'                      generation of ABC-SMC, with each vector being of length
#'                      \code{npart};}
#' \item{\code{accrate}:}{ a \code{vector} of length \code{nrow(tols)} containing the
#'                      acceptance rates for each generation of ABC.}
#' }
#'

ABC_SMC <- function(npart, tols, priors, func, data) {

    ## set up output objects
    accrate <- rep(NA, nrow(tols))
    pars <- list(); out <- list(); weights <- list();
    
    ## run sequential algorithm
    for(t in 1:nrow(tols)) {
        ## print progress to the screen
        cat(paste0("Generation ", t, ":\n"))
        
        ## set up inputs and outputs
        pars[[t]] <- matrix(NA, npart, nrow(priors))
        out[[t]] <- matrix(NA, npart, length(data))
        accrate[t] <- 0
        weightsNew <- rep(NA, npart)
        
<<<<<<< HEAD
        ## simulate particles
        for(i in 1:npart) {
            valid <- 0
            while(valid == 0) {
                if(t == 1) {
                    ## sample from prior
                    for(j in 1:nrow(priors)) {
                        pars[[t]][i, j] <- runif(1, priors[, 1], priors[, 2])
                    }
                } else {
                    ## sample from previous generation
                    k <- sample(1:npart, 1, prob = weights[[t - 1]])
                    pars[[t]][i, ] <- rmvnorm(1, 
                        mean = pars[[t - 1]][k, ], 
                        sigma = propCov
                    )
                }
                
                if(all(apply(cbind(pars[[t]][i, ], priors), 1, function(x) {
                    x[1] > x[2] & x[1] < x[3]
                }))) {
                    ## simulate from model
                    out[[t]][i, ] <- func(pars[[t]][i, ])
                    
                    ## check matching
                    if(all(!is.na(out[[t]][i, ]))) {
                        if(all(abs(out[[t]][i, ] - data) < tols[t, ])) {
                            valid <- 1
                        }
                    }
                }
                
                ## update counter
                accrate[t] <- accrate[t] + 1
            }
            
            if(t == 1) {
                weightsNew[i] <- 1
            } else {
                ## calculate unnormalised weight
                weightsNew[i] <- prod(apply(cbind(pars[[t]][i, ], priors), 1, function(x) {
                    dunif(x[1], x[2], x[3])
                }))
                weightsNew[i] <- weightsNew[i] / sum(weights[[t - 1]] * apply(pars[[t - 1]], 1, function(x) {
                    dmvnorm(pars[[t]][i, ], mean = x, sigma = propCov)
                }))
            }
            if(i %% 10 == 0) {
                cat(paste0("i = ", i, ", accrate = ", signif(i / accrate[t], 2), "\n"))
            }
=======
        ## set timer
        ptm <- proc.time()
        
        ## run generation
        if(!require(parallel)) {
            temp <- lapply(1:npart, runProp,
                t = t, priors = priors, 
                prevWeights = tempWeights, prevPars = tempPars, 
                propCov = propCov, tols = tols[t, ], data = data, func = func)
        } else  {
            temp <- mclapply(1:npart, runProp,
                t = t, priors = priors, 
                prevWeights = tempWeights, prevPars = tempPars, 
                propCov = propCov, tols = tols[t, ], data = data, func = func, mc.cores = mc.cores)
>>>>>>> 83c7044... Amended print rounding and set timer
        }
        
        ## set accrate
        accrate[t] <- npart / accrate[t] 
        
        ## set weights
        weights[[t]] <- weightsNew / sum(weightsNew)
        
        ## set proposal covariance
        propCov <- cov(pars[[t]]) * 2
<<<<<<< HEAD
=======
        
        ## stop timer
        ptm1 <- proc.time() - ptm
        
        ## print progress to the screen
        cat(paste0("Generation ", t, ", accrate = ", signif(accrate[t], 2), 
                   ", time = ", signif(ptm1[3], 2), " secs\n"))
>>>>>>> 83c7044... Amended print rounding and set timer
    }
    
    ## output results
    output <- list(pars = pars, output = out, weights = weights, accrate = accrate)
    class(output) <- "ABCSMC"
    output
}

        
<<<<<<< HEAD
=======
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
>>>>>>> 83c7044... Amended print rounding and set timer
