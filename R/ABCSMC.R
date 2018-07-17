#' @title Runs ABC-SMC algorithm
#'
#' @description Runs ABC-SMC algorithm of Toni et al. (2009)
#'
#' @export
#'
#' @param x         An \code{ABCSMC} object.
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
#'                      acceptance rates for each generation of ABC;}
#' \item{\code{tols}:}{ a copy of the \code{tols} input;}
#' \item{\code{priors}:}{ a copy of the \code{priors} input;}
#' \item{\code{data}:}{ a copy of the \code{data} input;}
#' \item{\code{func}:}{ a copy of the \code{func} input.}
#' }
#' @rdname ABCSMC

ABCSMC <- function(x, ...) {
    UseMethod("ABCSMC")
}

#' @rdname ABCSMC
#' @export

ABCSMC.ABCSMC <- function(x, tols) {
    
    ## check inputs
    stopifnot(class(x) == "ABCSMC")
    
    ## extract tolerances and check against new tolerances
    stopifnot(ncol(tols) == length(x$data))
    if(sum(apply(rbind(x$tols[nrow(x$tols), ], tols), 2, function(x) {
        sum(x[-1] >= x[1])
    })) > 0) {
        stop("New tolerances not less than original tolerances")
    }
    
    ## run ABC-SMC
    temp <- ABCSMC.default(nrow(x$pars[[1]]), tols, x$priors, x$func, x$data, 
                           prevPars = x$pars[[length(x$pars)]], 
                           prevWeights = x$weights[[length(x$weights)]])
    
    ## combine with original runs
    x$pars <- c(x$pars, temp$pars)
    x$output <- c(x$output, temp$output)
    x$weights <- c(x$weights, temp$weights)
    x$accrate <- c(x$accrate, temp$accrate)
    x$tols <- rbind(x$tols, tols)
    
    ## return new object
    x
}

#' @rdname ABCSMC
#' @export

ABCSMC.default <- function(npart, tols, priors, func, data, ...) {
    
    ## check inputs
    checkInput(npart, "numeric", 1, int = T)
    checkInput(tols, c("numeric", "matrix"))
    checkInput(priors, c("numeric", "matrix"), ncol = 2)
    checkInput(func, "function", 1)
    checkInput(data, c("vector", "numeric"))
    stopifnot(length(data) == ncol(tols))
    fargs <- formals(func)
    stopifnot(length(fargs) == 1)
    stopifnot(names(fargs) == "pars")
    stopifnot(all(apply(tols, 2, function(x) {
        all(diff(x) < 0)
    })))
    stopifnot(npart > 1)
    stopifnot(all(tols > 0))
    stopifnot(all(apply(priors, 1, diff) > 0))
    
    ## set timer
    ptm_ini <- proc.time()

    ## set up output objects
    accrate <- rep(NA, nrow(tols))
    pars <- list(); out <- list(); weights <- list();
    
    ## extract ellipsis arguments
    args <- list(...)
    
    if(exists("args$prevPars")) {
        stopifnot(exists("args$prevWeights"))
        pars[[1]] <- args$prevPars
        weights[[1]] <- args$prevWeights
        tols <- rbind(rep(NA, ncol(tols)), tols)
        accrate <- rep(NA, nrow(tols))
        out[[1]] <- NA
        propCov <- cov(pars[[1]]) * 2
        init <- 2
    } else {
        init <- 1
    }
    
    ## run sequential algorithm
    for(t in init:nrow(tols)) {
        ## print progress to the screen
        cat(paste0("Generation ", t, ":\n"))
        
        ## set up inputs and outputs
        pars[[t]] <- matrix(NA, npart, nrow(priors))
        out[[t]] <- matrix(NA, npart, length(data))
        accrate[t] <- 0
        weightsNew <- rep(NA, npart)
        
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
        }
        
        ## set accrate
        accrate[t] <- npart / accrate[t] 
        
        ## set weights
        weights[[t]] <- weightsNew / sum(weightsNew)
        
        ## set proposal covariance
        propCov <- cov(pars[[t]]) * 2
        
        ## stop timer
        ptm1 <- proc.time() - ptm
        
        ## print progress to the screen
        cat(paste0("Generation ", t, ", accrate = ", signif(accrate[t], 2), 
                   ", time = ", signif(ptm1[3], 2), " secs\n"))
    }
    
    ## stop timer
    ptm1 <- proc.time() - ptm_ini
    cat(paste0("\nFinal run time = ", signif(ptm1[3], 2), " secs\n"))
    
    # remove extraneous components if extending runs
    if(init > 1) {
        pars <- pars[-1]
        out <- out[-1]
        tols <- tols[-1, ]
        weights <- weights[-1]
        accrate <- accrate[-1]
    }
    
    ## output results
    output <- list(pars = pars, output = out, weights = weights, accrate = accrate,
                   tols = tols, priors = priors, data = data, func = func)
    class(output) <- "ABCSMC"
    output
}

