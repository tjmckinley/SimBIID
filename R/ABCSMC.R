#' @title Runs ABC-SMC algorithm
#'
#' @description     Runs ABC-SMC algorithm of Toni et al. (2009) for fitting 
#'                  infectious disease models to time series count data 
#'
#' @export
#'
#' @param x         An \code{ABCSMC} object or the number of particles (must be a positive integer).
#' @param data      A \code{data.frame} with a single row and columns containing the observed summary statistics
#'                  to match to. Columns must match to `tols`.
#' @param priors    A \code{data.frame} containing columns \code{parnames}, \code{dist}, \code{p1} and 
#'                  \code{p2}, with number of rows equal to the number of parameters. The column
#'                  \code{parname} simply gives names to each parameter for plotting and summarising.
#'                  Each entry in the \code{dist} column must contain one of \code{c("unif", "norm", "gamma")}, 
#'                  and the corresponding \code{p1} and \code{p2} entries relate to the hyperparameters 
#'                  (lower and upper bounds in the uniform case; mean and standard deviation in the 
#'                  normal case; and shape and rate in the gamma case).
#' @param func      Function that runs the simulator and checks whether the simulation matches the data. 
#'                  The first four arguments must be \code{pars}, \code{data}, \code{tols} and 
#'                  \code{iniStates}. If the simulations do not match the data then the function must 
#'                  return an \code{NA}, else it must returns a \code{vector} of simulated summary measures. 
#'                  In this latter case the output from the function must be a vector with length equal to 
#'                  \code{ncol(data)} and with entries in the same order as the columns of \code{data}.
#' @param iniStates A numerical vector of initial states for the infectious disease model.
#' @param tols 		A \code{data.frame} of tolerances, with the number of rows defining
#'                  the number of generations required, and columns defining the summary statistics
#'                  to match to. The columns must match to those in `data`.
#' @param parallel  A \code{logical} determining whether to use parallel processing or not.
#' @param mc.cores  Number of cores to use if using parallel processing.
#' @param ...       Further arguments to pass to \code{func}. (Not used if extending runs.)
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
#'                      \code{npart} x \code{ncol(data)};}
#' \item{\code{weights}:}{ a \code{list} of \code{vector} objects containing the particle
#'                      weights. Each element of the list corresponds to a
#'                      generation of ABC-SMC, with each vector being of length
#'                      \code{npart};}
#' \item{\code{accrate}:}{ a \code{vector} of length \code{nrow(tols)} containing the
#'                      acceptance rates for each generation of ABC;}
#' \item{\code{tols}:}{ a copy of the \code{tols} input;}
#' \item{\code{priors}:}{ a copy of the \code{priors} input;}
#' \item{\code{data}:}{ a copy of the \code{data} input;}
#' \item{\code{func}:}{ a copy of the \code{func} input;}
#' \item{\code{iniStates}}{ a vector of initial states for the infectious disease model;}
#' \item{\code{addargs}:}{ a copy of the \code{...} inputs.}
#' }
#' @rdname ABCSMC

ABCSMC <- function(x, ...) {
    UseMethod("ABCSMC")
}

#' @rdname ABCSMC
#' @export

ABCSMC.ABCSMC <- function(x, tols, parallel = F, mc.cores = NA) {
    
    ## check inputs
    if(class(x) != "ABCSMC"){
        stop("'x' not ABCSMC object")
    }
    
    ## extract tolerances and check against new tolerances
    if(ncol(tols) != ncol(x$data)){
        stop("ncol(tols) != ncol(x$data)")
    }
    if(!all(match(colnames(tols), colnames(x$data)) - 1:ncol(x$data) == 0)){
        stop("colnames(tols) does not match colnames(data)")
    }
    if(sum(apply(rbind(x$tols[nrow(x$tols), ], tols), 2, function(x) {
        sum(x[-1] >= x[1])
    })) > 0) {
        stop("New tolerances not less than original tolerances")
    }
    
    ## collect arguments
    tempargs <- list(
        x = nrow(x$pars[[1]]), 
        tols = tols, 
        priors = x$priors, 
        func = x$func, 
        data = x$data, 
        iniStates = x$iniStates,
        parallel = parallel,
        mc.cores = mc.cores, 
        prevPars = x$pars[[length(x$pars)]], 
        prevWeights = x$weights[[length(x$weights)]],
        genstart = nrow(x$tols) + 1
    )
    tempargs <- c(tempargs, x$addargs)
    
    ## run ABC-SMC
    temp <- do.call("ABCSMC.default", tempargs)
    
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

ABCSMC.default <- function(x, data, priors, func, iniStates, tols, parallel = F, mc.cores = NA, ...) {
    
    
    ## check missing arguments
    if(missing(x) & missing(data)){
        stop("'x' and 'data' can't both be missing")
    }
    if(missing(priors)){
        stop("'priors' argument missing")
    }
    if(missing(func)){
        stop("'func' argument missing")
    }
    if(missing(iniStates)){
        stop("'iniStates' argument missing")
    }
    if(missing(tols)){
        stop("'tols' argument missing")
    }
    
    ## check inputs
    npart <- x
    checkInput(parallel, c("vector", "logical"), 1)
    if(parallel) {
        if(!require(parallel)) {
            stop("Must have 'parallel' package installed to use parallelisation")
        }
        nc <- detectCores()
        nc <- ifelse(is.na(nc), 1, nc)
        if(!is.na(mc.cores[1])) {
            checkInput(mc.cores, "numeric", 1, int = T)
            mc.cores <- min(nc, mc.cores)
        } else {
            mc.cores <- nc
        }
        parallel <- (mc.cores > 1)
        cat(paste0("Number of cores: ", mc.cores, "\n"))
    }
    checkInput(npart, "numeric", 1, int = T, gt = 1)
    checkInput(tols, "data.frame", gte = 0)
    if(!all(sapply(tols, is.numeric))) {
        stop("'tols' must be numeric")
    }
    checkInput(priors, "data.frame", ncol = 4)
    checkInput(func, "function", 1)
    checkInput(data, "data.frame", nrow = 1)
    if(!all(sapply(data, is.numeric))) {
        stop("'data' must be numeric")
    }
    if(ncol(data) != ncol(tols)){
        stop("Number of columns of 'data' and 'tols' must match")
    }
    if(!all(match(colnames(data), colnames(tols)) - 1:ncol(data) == 0)){
        stop("colnames(tols) !- colnames(data)")
    }
    fargs <- formals(func)
    if(length(fargs) < 3){
        stop("Number of arguments of 'func' must be at least 4")
    }
    if(!all(match(names(fargs)[1:3], c("pars", "data", "tols", "iniStates")) - 1:4 == 0)){
        stop("First four arguments of 'func' must be: 'pars', 'data', 'tols' and 'iniStates'")
    }
    if(!all(apply(tols, 2, function(x) {
        all(diff(x) < 0)
    }))){
        stop("'tols' cannot increase")
    }
    ## check iniStates
    checkInput(iniStates, "numeric", int = T, gte = 0)
    checkInput(sum(iniStates), "numeric", int = T, gt = 1)
    
    ## check priors
    if(!all(sort(match(colnames(priors), c("parnames", "dist", "p1", "p2"))) - 1:4 == 0)){
        stop("Column names of 'priors' must be: 'parnames', 'dist', 'p1' and 'p2'")
    }
    priors <- select(priors, parnames, dist, p1, p2)
    checkInput(priors$parnames, "character")
    checkInput(priors$dist, "character")
    checkInput(priors$p1, "numeric")
    checkInput(priors$p2, "numeric")
    if(!all(priors$dist %in% c("unif", "norm", "gamma"))){
        stop("'priors' must be of form: 'unif', 'norm' or 'gamma'")
    }
    temp <- priors[priors$dist == "unif", , drop = F]
    if(nrow(temp) > 0) {
        ## check uniform bounds correct
        if(!all(apply(temp[, 3:4, drop = F], 1, diff) > 0)){
            stop("Priors: uniform bounds in wrong order")
        }
    }
    temp <- priors[priors$dist == "norm", , drop = F]
    if(nrow(temp) > 0) {
        ## check normal hyperparameters correct
        if(!all(temp$p2 > 0)){
            stop("Priors: normal variances must be > 0")
        }
    }
    temp <- priors[priors$dist == "gamma", , drop = F]
    if(nrow(temp) > 0) {
        ## check gamma bounds correct
        if(!all(temp$p1 > 0) | !all(temp$p2 > 0)){
            stop("Priors: gamma hyperparameters must be > 0")
        }
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
    
    ## extract arguments for "func"
    fargs <- fargs[is.na(match(names(fargs), c("pars", "data", "tols", "iniStates")))]
    if(length(fargs) > 0) {
        fargs1 <- match(names(fargs), names(args))
        if(!all(!is.na(fargs1))){
            stop(paste0("Need to include: ", paste(fargs[is.na(fargs1)], collapse = ", "), 
                        " arguments in function call"))
        }
        fargs <- args[fargs1]
    }
    
    if(exists("prevPars", where = args)) {
        if(!exists("prevWeights", where = args)){
            stop("'prevweights' does not exist")
        }
        if(!exists("genstart", where = args)){
            stop("'genstart' doe not exist")
        }
        ## set up dummy objects (removed at the end)
        pars[[1]] <- args$prevPars
        weights[[1]] <- args$prevWeights
        tols <- rbind(rep(NA, ncol(tols)), tols)
        out[[1]] <- NA
        propCov <- cov(pars[[1]]) * 2
        init <- 2
        genstart <- args$genstart - 2
        ## remove additional arguments
        args <- args[-match("prevPars", names(args))]
        args <- args[-match("prevWeights", names(args))]
        args <- args[-match("genstart", names(args))]
    } else {
        init <- 1
        genstart <- 0
    }
    
    ## save tols
    orig_tols <- tols
    tols <- as.matrix(tols)
    
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
                propCov = propCov, tols = tols[t, ], data = as.numeric(data[1, ]), 
                iniStates = iniStates,
                func = func, func_args = fargs)
        } else  {
            temp <- mclapply(1:npart, runProp,
                t = t, priors = priors, 
                prevWeights = tempWeights, prevPars = tempPars, 
                propCov = propCov, tols = tols[t, ], data = as.numeric(data[1, ]), 
                iniStates = iniStates,
                func = func, func_args = fargs, 
                mc.cores = mc.cores)
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
        colnames(out[[t]]) <- colnames(data)
        
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
        orig_tols <- orig_tols[-1, , drop = F]
        weights <- weights[-1]
    }
    
    ## output results
    output <- list(pars = pars, output = out, weights = weights, accrate = accrate,
                   tols = orig_tols, priors = orig_priors, data = data,
                   func = func, iniStates = iniStates, addargs = args)
    class(output) <- "ABCSMC"
    output
}

