#' @title Runs ABC-SMC algorithm
#'
#' @description Runs ABC-SMC algorithm of Toni et al. (2009)
#'
#' @export
#'
#' @param x         An \code{ABCSMC} object.
#' @param npart     The number of particles (must be a positive integer).
#' @param tols 		A \code{matrix} of tolerances, with the number of rows defining
#'                  the number of generations required, and the number of columns
#'                  defining the number of summary statistics / data points to match to.
#' @param priors    A \code{data.frame} containing columns \code{parnames}, \code{dist}, \code{p1} and 
#'                  \code{p2}, with number of rows equal to the number of parameters. The column
#'                  \code{parname} simply gives names to each parameter for plotting and summarising.
#'                  Each entry in the \code{dist} column must contain one of \code{c("unif", "norm", "gamma")}, 
#'                  and the corresponding \code{p1} and \code{p2} entries relate to the hyperparameters 
#'                  (lower and upper bounds in the uniform case; mean and standard deviation in the 
#'                  normal case; and shape and rate in the gamma case).
#'                  two columns containing the lower and upper bounds respectively.
#' @param func      A function taking a single argument \code{pars} that runs the simulator
#'                  and returns the simulated summary measures against which to compare. The output from
#'                  the function must be a vector with length equal to \code{nrow(data)} and with entries
#'                  in the same order as the rows of \code{data}.
#' @param data      A \code{data.frame} with two columns containing the observed summary statistics to match to. 
#'                  The first column must be called \code{outnames} and contain the output names, and the
#'                  second column must be called \code{values} and contain the observations.
#' @param parallel  A \code{logical} determining whether to use parallel processing or not.
#' @param mc.cores  Number of cores to use if using parallel processing.
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

ABCSMC.ABCSMC <- function(x, tols, parallel = F, mc.cores = NA) {
    
    ## check inputs
    stopifnot(class(x) == "ABCSMC")
    
    ## extract tolerances and check against new tolerances
    stopifnot(ncol(tols) == nrow(x$data))
    if(sum(apply(rbind(x$tols[nrow(x$tols), ], tols), 2, function(x) {
        sum(x[-1] >= x[1])
    })) > 0) {
        stop("New tolerances not less than original tolerances")
    }
    
    ## run ABC-SMC
    temp <- ABCSMC.default(nrow(x$pars[[1]]), tols, x$priors, x$func, x$data, 
                           parallel = parallel,
                           mc.cores = mc.cores, 
                           prevPars = x$pars[[length(x$pars)]], 
                           prevWeights = x$weights[[length(x$weights)]],
                           genstart = nrow(x$tols) + 1)
    
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

ABCSMC.default <- function(npart, tols, priors, func, data, parallel = F, mc.cores = NA, ...) {
    
    ## check inputs
    stopifnot(checkInput(parallel, c("vector", "logical"), 1))
    if(parallel) {
        if(!require(parallel)) {
            stop("Must have 'parallel' package installed to use parallelisation")
        }
        nc <- detectCores()
        nc <- ifelse(is.na(nc), 1, nc)
        if(!is.na(mc.cores[1])) {
            stopifnot(checkInput(mc.cores, "numeric", 1, int = T))
            mc.cores <- min(nc, mc.cores)
        } else {
            mc.cores <- nc
        }
        parallel <- (mc.cores > 1)
        cat(paste0("Number of cores: ", mc.cores, "\n"))
    }
    stopifnot(checkInput(npart, "numeric", 1, int = T))
    stopifnot(checkInput(tols, c("numeric", "matrix")))
    stopifnot(checkInput(priors, "data.frame", ncol = 4))
    stopifnot(checkInput(func, "function", 1))
    stopifnot(checkInput(data, "data.frame"))
    stopifnot(all(sort(match(colnames(data), c("outnames", "values"))) - 1:2 == 0))
    data <- select(data, outnames, values)
    stopifnot(checkInput(data$outnames, "character"))
    stopifnot(checkInput(data$values, "numeric"))
    stopifnot(nrow(data) == ncol(tols))
    fargs <- formals(func)
    stopifnot(length(fargs) == 1)
    stopifnot(names(fargs) == "pars")
    stopifnot(all(apply(tols, 2, function(x) {
        all(diff(x) < 0)
    })))
    stopifnot(npart > 1)
    stopifnot(all(tols > 0))
    
    ## check priors
    stopifnot(all(sort(match(colnames(priors), c("parnames", "dist", "p1", "p2"))) - 1:4 == 0))
    priors <- select(priors, parnames, dist, p1, p2)
    stopifnot(checkInput(priors$parnames, "character"))
    stopifnot(checkInput(priors$dist, "character"))
    stopifnot(checkInput(priors$p1, "numeric"))
    stopifnot(checkInput(priors$p2, "numeric"))
    stopifnot(all(priors$dist %in% c("unif", "norm", "gamma")))
    temp <- priors[priors$dist == "unif", , drop = F]
    if(nrow(temp) > 0) {
        ## check uniform bounds correct
        stopifnot(all(apply(temp[, 3:4, drop = F], 1, diff) > 0))
    }
    temp <- priors[priors$dist == "norm", , drop = F]
    if(nrow(temp) > 0) {
        ## check normal hyperparameters correct
        stopifnot(all(temp$p2 > 0))
    }
    temp <- priors[priors$dist == "gamma", , drop = F]
    if(nrow(temp) > 0) {
        ## check gamma bounds correct
        stopifnot(all(temp$p1 > 0))
        stopifnot(all(temp$p2 > 0))
    }
    orig_priors <- priors
    priors$ddist <- paste0("d", priors$dist)
    priors$dist <- paste0("r", priors$dist)    
    
    ## set timer
    ptm_ini <- proc.time()

    ## set up output objects
    accrate <- rep(NA, nrow(tols))
    pars <- list(); out <- list(); weights <- list();
    
    ## extract ellipsis arguments
    args <- list(...)
    
    if(exists("prevPars", where = args)) {
        stopifnot(exists("prevWeights", where = args))
        stopifnot(exists("genstart", where = args))
        ## set up dummy objects (removed at the end)
        pars[[1]] <- args$prevPars
        weights[[1]] <- args$prevWeights
        tols <- rbind(rep(NA, ncol(tols)), tols)
        out[[1]] <- NA
        propCov <- cov(pars[[1]]) * 2
        init <- 2
        genstart <- args$genstart - 2
    } else {
        init <- 1
        genstart <- 0
    }
    
    ## run sequential algorithm
    for(t in init:nrow(tols)) {  
        ## set up arguments
        if(t == 1) {
            tempWeights <- rep(1, npart)
            tempPars <- rep(1, npart)
        } else {
            tempWeights <- weights[[t - 1]]
            tempPars <- pars[[t - 1]]
        }
        
        ## set timer
        ptm <- proc.time()
        
        ## run generation
        if(!parallel) {
            temp <- lapply(1:npart, runProp,
                t = t, priors = priors, 
                prevWeights = tempWeights, prevPars = tempPars, 
                propCov = propCov, tols = tols[t, ], data = data$values, func = func)
        } else  {
            temp <- mclapply(1:npart, runProp,
                t = t, priors = priors, 
                prevWeights = tempWeights, prevPars = tempPars, 
                propCov = propCov, tols = tols[t, ], data = data$values, func = func, mc.cores = mc.cores)
        }
        
        ## extract relative components
        weights[[t]] <- map_dbl(temp, "weightsNew")
        weights[[t]] <- weights[[t]] / sum(weights[[t]])
        pars[[t]] <- map(temp, "pars")
        pars[[t]] <- do.call("rbind", pars[[t]])
        out[[t]] <- map(temp, "out")
        out[[t]] <- do.call("rbind", out[[t]])
        accrate[t] <- npart / sum(map_dbl(temp, "accrate"))
        
        ## set names
        colnames(pars[[t]]) <- priors$parnames
        colnames(out[[t]]) <- data$outnames
        
        ## set proposal covariance
        propCov <- cov(pars[[t]]) * 2
        
        ## stop timer
        ptm1 <- proc.time() - ptm
        
        ## print progress to the screen
        cat(paste0("Generation ", t + genstart, ", accrate = ", signif(accrate[t], 2), 
                   ", time = ", signif(ptm1[3], 2), " secs\n"))
    }
    
    ## stop timer
    ptm1 <- proc.time() - ptm_ini
    cat(paste0("\nFinal run time = ", signif(ptm1[3], 2), " secs\n"))
    
    ## remove extraneous components if extending runs
    if(init > 1) {
        pars <- pars[-1]
        out <- out[-1]
        tols <- tols[-1, ]
        weights <- weights[-1]
    }
    
    ## output results
    output <- list(pars = pars, output = out, weights = weights, accrate = accrate,
                   tols = tols, priors = orig_priors, data = data, func = func)
    class(output) <- "ABCSMC"
    output
}

