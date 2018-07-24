#' @title Produces ABC reference table
#'
#' @description Runs ABC reference table
#'
#' @export
#'
#' @param npart     The number of particles (must be a positive integer).
#' @param priors    A \code{data.frame} containing columns \code{parnames}, \code{dist}, \code{p1} and 
#'                  \code{p2}, with number of rows equal to the number of parameters. The column
#'                  \code{parname} simply gives names to each parameter for plotting and summarising.
#'                  Each entry in the \code{dist} column must contain one of \code{c("unif", "norm", "gamma")}, 
#'                  and the corresponding \code{p1} and \code{p2} entries relate to the hyperparameters 
#'                  (lower and upper bounds in the uniform case; mean and standard deviation in the 
#'                  normal case; and shape and rate in the gamma case).
#' @param func      Function that runs the simulator. The first argument must be \code{pars}. The function
#'                  must return a \code{vector} of simulated summary measures, or a missing value (\code{NA})
#'                  if there is an error. The output from the function must be a vector with length equal 
#'                  to \code{nrow(data)} and with entries in the same order as the rows of \code{data}.
#' @param data      A \code{data.frame} with a single row and columns containing the observed summary
#'                  statistics in order to extract output names.
#' @param parallel  A \code{logical} determining whether to use parallel processing or not.
#' @param mc.cores  Number of cores to use if using parallel processing.
#' @param ...       Extra arguments to be passed to \code{func}.
#'
#' @return An \code{data.frame} object with \code{npart} rows, where the first \code{p} columns correspond to 
#'         the proposed parameters, and the remaining columns correspond to the simulated outputs.

ABCRef <- function(npart, priors, func, data, parallel = F, mc.cores = NA, ...) {
    
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
    stopifnot(checkInput(priors, "data.frame", ncol = 4))
    stopifnot(checkInput(func, "function", 1))
    stopifnot(checkInput(data, "data.frame", nrow = 1))
    stopifnot(all(sapply(data, is.numeric)))
    fargs <- formals(func)
    stopifnot(length(fargs) >= 1)
    stopifnot(names(fargs)[1] == "pars")
    stopifnot(npart > 1)
    
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
    priors$ddist <- paste0("d", priors$dist)
    priors$dist <- paste0("r", priors$dist)   
    
    ## extract arguments for "func"
    fargs <- fargs[is.na(match(names(fargs), "pars"))]
    if(length(fargs) > 0) {
        args <- list(...)
        fargs <- match(names(fargs), names(args))
        stopifnot(all(!is.na(fargs)))
        fargs <- args[fargs]
    }
    
    ## set timer
    ptm_ini <- proc.time()
    ptm <- ptm_ini

    ## set up output objects
    pars <- list(); out <- list(); 
    
    ## run generation
    if(!parallel) {
        temp <- lapply(1:npart, runRef, priors = priors, func = func, func_args = fargs)
    } else  {
        temp <- mclapply(1:npart, runRef, priors = priors, func = func, func_args = fargs, 
                         mc.cores = mc.cores)
    }
    
    ## extract relative components
    out <- map(temp, "out")
    pars <- map(temp, "pars")
    pars <- pars[!map_lgl(out, function(x) is.na(x[1]))]
    out <- out[!map_lgl(out, function(x) is.na(x[1]))]
    
    while(length(out) != npart) {
        
        ## stop timer
        ptm1 <- proc.time() - ptm
        ptm <- ptm1
        
        ## print progress
        cat(paste0("Current pars = ", length(out), ", time = ", signif(ptm1[3], 2), " secs\n"))
        
        ## run generation
        if(!parallel) {
            temp <- lapply(1:npart, runRef, priors = priors, func = func)
        } else  {
            temp <- mclapply(1:npart, runRef, priors = priors, func = func, mc.cores = mc.cores)
        }
        ## extract relative components
        out <- c(out, map(temp, "out"))
        pars <- c(pars, map(temp, "pars"))
        pars <- pars[!map_lgl(out, function(x) is.na(x[1]))]
        out <- out[!map_lgl(out, function(x) is.na(x[1]))]
    }
    out <- do.call("rbind", out)
    pars <- do.call("rbind", pars)
    out <- as.data.frame(out)
    pars <- as.data.frame(pars)
    
    ## set names
    colnames(pars) <- priors$parnames
    stopifnot(ncol(out) == ncol(data))
    colnames(out) <- colnames(data)
    
    ## stop timer
    ptm1 <- proc.time() - ptm_ini
    
    ## print progress to the screen
    cat(paste0("Final run time = ", signif(ptm1[3], 2), " secs\n"))
    
    return(cbind(pars, out))
}

