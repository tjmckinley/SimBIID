#' @title Runs particle MCMC algorithm
#'
#' @description Runs particle MCMC algorithm using alive particle filter for fitting 
#'              infectious disease models to time series count data
#'
#' @export
#'
#' @param x 		    A \code{PMCMC} object, or a \code{data.frame} containing time series count data, with the first column called
#'                      \code{time}, followed by columns of time-series counts.
#' @param priors        A \code{data.frame} containing columns \code{parnames}, \code{dist}, \code{p1} and 
#'                      \code{p2}, with number of rows equal to the number of parameters. The column
#'                      \code{parname} simply gives names to each parameter for plotting and summarising.
#'                      Each entry in the \code{dist} column must contain one of \code{c("unif", "norm", "gamma")}, 
#'                      and the corresponding \code{p1} and \code{p2} entries relate to the hyperparameters 
#'                      (lower and upper bounds in the uniform case; mean and standard deviation in the 
#'                      normal case; and shape and rate in the gamma case).
#' @param func          A \code{SimBIID_model} object or an \code{XPtr} to simulation function. If the latter, 
#'                      then this function must take the following arguments in order: 
#'                      \itemize{
#'                      \item{\code{NumericVector pars}:}{ a vector of parameters;}
#'                      \item{\code{double tstart}:}{ the start time;}
#'                      \item{\code{double tstop}:}{ the end time;}
#'                      \item{\code{IntegerVector u}:}{ a vector of states;}
#'                      \item{\code{IntegerVector tols}:}{ a vector of tolerances;}
#'                      \item{\code{IntegerVector counts}:}{ a vector of observed states;}
#'                      \item{\code{IntegerVector whichind}:}{ a vector the same length as \code{counts}, 
#'                      indicating which elements of \code{u} correspond to which elements of \code{counts}.
#'                      Must index from 0 NOT 1.}}
#' @param u             A \code{data.frame} of initial states, with a single row, and columns defining the 
#'                      compartments.
#' @param npart         An integer specifying the number of particles for the alive particle filter.
#' @param tols          Tolerances for matching data during ABC. Defaults to matching every count column
#'                      of \code{data} exactly.
#' @param whichind      Not needed if \code{func} is a \code{SimBIID_model} object (in this case it matches
#'                      column names of \code{data} to column names of \code{u}. Otherwise this must be a
#'                      vector relating which elements of \code{u} match to columns \code{2:ncol(data)}
#'                      of `data`. If left as \code{NULL} then defaults to elements \code{1:length(u)}
#'                      or returns an error. Must be same length as \code{tols} and must index from 1 (it's 
#'                      converted to C indexing internally).
#' @param iniPars       Vector of initial parameter values. If left unspecified, then these are 
#'                      sampled from the prior distributions.
#' @param fixpars       A logical determining whether to fix the input parameters (useful for 
#'                      determining the variance of the marginal likelihood estimates).
#' @param niter         An integer specifying the number of iterations to run the MCMC.
#' @param nprintsum     Prints summary of MCMC to screen every \code{nprintsum} iterations. 
#'                      Also defines how often adaptive scaling of proposal variances occur.
#' @param nmultskip     Upper bound for how many simulations to run before rejecting proposal.
#' @param adapt         Logical determining whether to use adaptive proposal or not.
#' @param propVar       A numeric (npars x npars) matrix with log (or logistic) covariances to use
#'                      as (initial) proposal matrix. If left unspecified then defaults to 
#'                      \code{diag(nrow(priors)) * (0.1 ^ 2) / nrow(priors)}.
#' @param adaptmixprop  Mixing proportion for adaptive proposal.
#' @param nupdate       Controls when to start adaptive update.
#'
#' @details             Function runs the alive particle filter for a given model. 
#'                      If running with \code{fixpars = T} then this runs \code{niter} simulations
#'                      using fixed parameter values. This can be used to optimise the number of 
#'                      particles after a training run.
#'
#' @return If the code throws an error, then it returns a missing value (\code{NA}). If 
#'         \code{fixpars = T} it returns a list of length 2 containing:
#' \itemize{
#'      \item{\code{output}:}{ a matrix with two columns. The first contains the simulated
#'          log-likelihood, and the second is a binary indicator relating to whether the
#'          simulation was 'skipped' or not (1 = skipped, 0 = not skipped);}
#'      \item{\code{pars}:}{ a vector of parameters used for the simulations.}
#' }
#' If \code{fixpars = F}, the routine returns a \code{PMCMC} object, essentially a 
#'          \code{list} containing:
#' \itemize{
#'  \item{\code{pars}:}{ an \code{mcmc} object containing posterior samples for the parameters;}
#'  \item{\code{tols}:}{ tolerances;}
#'  \item{\code{whichind}:}{ matching indicators;}
#'  \item{\code{u}:}{ a copy of the \code{u} input;}
#'  \item{\code{skiprate}:}{ the cumulative skip rate;}
#'  \item{\code{accrate}:}{ the cumulative acceptance rate;}
#'  \item{\code{nmultskip}:}{ the chosen value of \code{nmultskip};}
#'  \item{\code{npart}:}{ the chosen number of particles;}
#'  \item{\code{time}:}{ the time taken to run the routine (in seconds);}
#'  \item{\code{propVar}:}{ the proposal covariance for the parameter updates;}
#'  \item{\code{data}:}{ a copy of the \code{x} input;}
#'  \item{\code{priors}:}{ a copy of the \code{priors} input;}
#'  \item{\code{func}:}{ a copy of the \code{func} input.}
#' }
#' @rdname PMCMC
#' 

PMCMC <- function(x, ...) {
    UseMethod("PMCMC")
}

#' @rdname PMCMC
#' @export

PMCMC.PMCMC <- function(x, niter = 1000, nprintsum = 100, 
                        adapt = T, adaptmixprop = 0.05, 
                        nupdate = 100) {
    ## check object
    if(class(x) != "PMCMC") {
        stop("'x' not a PMCMC object")
    }
    if(length(x) <= 2){
        stop("'x' not a valid PMCMC object to continue running...")
    }
    
    ## collect arguments
    tempargs <- list(
        x = x$data, 
        priors = x$priors, 
        func = x$func, 
        u = x$u, 
        npart = x$npart, 
        tols = x$tols, 
        whichind = x$whichind, 
        iniPars = x$pars[nrow(x$pars), -c((ncol(x$pars) - 2):ncol(x$pars))], 
        fixpars = F, 
        niter = niter, 
        nprintsum = nprintsum, 
        nmultskip = x$nmultskip, 
        adapt = adapt, 
        propVar = x$propVar, 
        adaptmixprop = adaptmixprop, 
        nupdate = nupdate
    )
    
    ## run PMCMC
    temp <- do.call("PMCMC.default", tempargs)
    
    ## combine with original runs
    x$pars <- as.mcmc(rbind(as.matrix(x$pars), as.matrix(temp$pars)))
    
    ## return new object
    x
}
    
#' @rdname PMCMC
#' @export

PMCMC.default <- function(
    x, priors, func, u, npart = 100,
    tols = rep(0, ncol(x) - 1), whichind = NULL, 
    iniPars = NA, fixpars = F, 
    niter = 1000, nprintsum = 100, nmultskip = 1000, 
    adapt = T, propVar = NA, adaptmixprop = 0.05, nupdate = 100
) {
    
    ## check inputs are present
    if(missing(x)){
        stop("'x' argument missing")
    }
    if(missing(priors)){
        stop("'priors' argument missing")
    }
    if(missing(u)){
        stop("'u' argument missing")
    }
    if(missing(func)){
        stop("'func' argument missing")
    } 
    
    ## check data set
    data <- x
    checkInput(data, "data.frame")
    if(colnames(data)[1] != "time"){
        stop("First column of 'data' must be 'time'")
    }
    if(ncol(data) < 1) {
        stop("Must have at least one count column in 'data'")
    }
    checkInput(data$time, "numeric")
    for(j in 2:ncol(data)) {
        checkInput(data[, j, drop = T], "numeric", int = T)
    }   
    
    ## check priors 
    checkInput(priors, "data.frame", ncol = 4)
    if(!all(sort(match(colnames(priors), c("parnames", "dist", "p1", "p2"))) - 1:4 == 0)){
        stop("'priors' must have column names: 'parnames', 'dist', 'p1', 'p2'")
    }
    priors <- select(priors, parnames, dist, p1, p2)
    checkInput(priors$parnames, "character")
    checkInput(priors$dist, "character")
    checkInput(priors$p1, "numeric")
    checkInput(priors$p2, "numeric")
    checkInput(priors$dist, inSet = c("unif", "norm", "gamma"))
    temp <- priors[priors$dist == "unif", , drop = F]
    if(nrow(temp) > 0) {
        ## check uniform bounds correct
        if(!all(apply(temp[, 3:4, drop = F], 1, diff) > 0)) {
            stop("Priors: 'uniform' bounds not in correct order")
        }
    }
    temp <- priors[priors$dist == "norm", , drop = F]
    if(nrow(temp) > 0) {
        ## check normal hyperparameters correct
        if(!all(temp$p2 > 0)){
            stop("Priors: 'normal' variance not positive")
        }
    }
    temp <- priors[priors$dist == "gamma", , drop = F]
    if(nrow(temp) > 0) {
        ## check gamma bounds correct
        if(!all(temp$p1 > 0) | !all(temp$p2 > 0)){
            stop("Priors: 'gamma' hyperparameters not positive")
        }
    }
    orig_priors <- priors
    priors$parnames <- NULL
    priors$dist <- match(priors$dist, c("unif", "norm", "gamma"))
    priors <- as.matrix(priors)
    
    ## check function
    funcorig <- func
    if(class(func) != "XPtr" & class(func) != "SimBIID_model"){
        stop("'func' not a 'SimBIID_model' object or an 'XPtr' object")
    }
    if(class(func) == "XPtr"){
        checkXPtr(func, "SEXP", c("NumericVector", "double", "double", "IntegerVector",
           "IntegerVector", "IntegerVector", "IntegerVector"))
    } else {
        ## check model
        if(func$tspan){
            warning("'SimBIID_model' object will have 'tspan' set to F")
        }
        ## generate model
        func <- mparseRcpp(
            transitions = func$transitions,
            compartments = func$compartments,
            pars = func$pars,
            matchCrit = T,
            addVars = func$addVars,
            stopCrit = func$stopCrit,
            tspan = F,
            runFromR = F
        )
        ## compile model
        func <- compileRcpp(func)
        
        ## set whichind
        checkInput(colnames(data)[-1], inSet = colnames(u))
        whichind <- match(colnames(data)[-1], colnames(u))
    }
    
    ## check initial conditions
    if(!any(is.na(iniPars))) {
        checkInput(iniPars, c("numeric", "vector"), nrow(priors), naAllow = T)
    } else {
        iniPars <- rep(NA, nrow(priors))
    }
    ## check u
    checkInput(u, "data.frame", nrow = 1)
    for(j in ncol(u)){
        checkInput(u[, j], c("vector", "numeric"), int = T, gte = 0)
    }
    checkInput(sum(u[1, ]), "numeric", int = T, gt = 1)
    uorig <- u
    u <- unlist(u)
    
    ## check proposal variances
    if(is.na(propVar[1])) {
        propVar <- diag(nrow(priors))
        ## adjust for scaling parameter for initial iterations
        propVar <- propVar * ((0.1 ^ 2) / nrow(propVar))
        propVar <- propVar / ((2.562 ^ 2) / nrow(propVar))
    } else {
        checkInput(propVar, c("numeric", "matrix"), nrow = nrow(priors), ncol = nrow(priors))
    }
    
    ## check tolerance argument
    checkInput(tols, c("numeric", "vector"), ncol(data) - 1)
    if(!all(tols >= 0)) {
        stop("'tols' must be >= 0")
    }
    
    ## check whichind
    whichindorig <- whichind
    if(!is.null(whichind)) {
        checkInput(
            whichind, c("numeric", "vector"), ncol(data) - 1, 
            int = T, inSet = 1:length(u), uni = T
        )
        whichind <- whichind - 1
    } else {
        whichind <- 1:length(u) - 1
        if((ncol(data) - 1) == length(u)){
            stop("length of 'u' does not match number of columns of 'data' (-1)")
        }
    }
    
    ## check runtime arguments
    checkInput(fixpars, c("logical", "vector"), 1)
    if(fixpars & any(is.na(iniPars))) {
        stop("Must input initial parameters if fixing parameters.")
    }
    checkInput(niter, c("numeric", "vector"), 1, int = T, gt = 0)
    checkInput(npart, c("numeric", "vector"), 1, int = T, gt = 0)
    checkInput(nprintsum, c("numeric", "vector"), 1, int = T, gt = 0)
    checkInput(nmultskip, c("numeric", "vector"), 1, int = T, gt = 1)
    
    ## check adaptive update and proposal covariance matrices
    checkInput(adapt, c("logical", "vector"), 1)
    checkInput(adaptmixprop, c("numeric", "vector"), 1, gt = 0, lt = 1)
    checkInput(nupdate, c("numeric", "vector"), 1, int = T, gt = 0)
    
    ## run function
    output <- PMCMC_cpp(as.matrix(data), priors, orig_priors$parnames, iniPars, propVar, niter, npart, 
                    adaptmixprop, tols, whichind, nprintsum, nmultskip, nupdate, as.numeric(fixpars), 
                    as.numeric(adapt), u, func)
    
    ## check to see if code has stopped afer initialisation
    if(length(output[[1]]) == 1) {
        return(NA)
    }
    if(length(output) == 2) {
        names(output) <- c("output", "pars")
        return(output)
    }
    
    ## convert output into correct format
    colnames(output[[1]]) <- c(orig_priors$parnames, "logPost", "nsims", "nsimsProp")
    output[[1]] <- as.mcmc(output[[1]])
    
    ## finalise output and set names
    output <- c(output[1], list(tols), list(whichindorig), list(uorig), output[-1], list(data), list(orig_priors), list(funcorig))
    names(output) <- c("pars", "tols", "whichind", "u", "skiprate", "accrate", 
        "nmultskip", "npart", "time", "propVar", "data", "priors", "func")
        
    ## export class and object
    class(output) <- "PMCMC"
    output
}

