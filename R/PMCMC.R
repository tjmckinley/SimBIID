#' @title Runs particle MCMC algorithm
#'
#' @description Runs particle MCMC algorithm using alive particle filter for fitting 
#'              infectious disease models to time series count data
#'
#' @export
#'
#' @param dataset 		Data frame containing time series count data, of form (time, counts*).
#' @param priors        A \code{data.frame} containing columns \code{parnames}, \code{dist}, \code{p1} and 
    #'                  \code{p2}, with number of rows equal to the number of parameters. The column
    #'                  \code{parname} simply gives names to each parameter for plotting and summarising.
    #'                  Each entry in the \code{dist} column must contain one of \code{c("unif", "norm", "gamma")}, 
    #'                  and the corresponding \code{p1} and \code{p2} entries relate to the hyperparameters 
    #'                  (lower and upper bounds in the uniform case; mean and standard deviation in the 
    #'                  normal case; and shape and rate in the gamma case).
#' @param iniStates     A numerical vector of initial states for the infectious disease model.
#' @param func          SEXP pointer to simulation function.
#' @param iniPars       Vector of initial parameter values. If left unspecified, then these are 
#'                      sampled from the prior distributions.
#' @param tol           Tolerance for matching data during ABC.
#' @param fixpars       A logical determining whether to fix the input parameters (useful for 
#'                      determining the variance of the marginal likelihood estimates).
#' @param niter         An integer specifying the number of iterations to run the MCMC.
#' @param npart         An integer specifying the number of particles for the alive particle filter.
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
#'      \item{\code{output}}{ a matrix with two columns. The first contains the simulated
#'          log-likelihood, and the second is a binary indicator relating to whether the
#'          simulation was 'skipped' or not (1 = skipped, 0 = not skipped);}
#'      \item{\code{pars}}{ a vector of parameters used for the simulations.}
#' }
#' If \code{fixpars = F}, the routine returns a \code{PMCMC} object, essentially a 
#'          \code{list} containing:
#' \itemize{
#'  \item{\code{pars}}{an \code{mcmc} object containing posterior samples for the parameters;}
#'  \item{\code{tol}}{tolerance for the NTG birds;}
#'  \item{\code{iniStates}}{a vector of initial states for the infectious disease model;}
#'  \item{\code{skiprate}}{the cumulative skip rate;}
#'  \item{accrate}{the cumulative acceptance rate;}
#'  \item{nmultskip}{the chosen value of \code{nmultskip};}
#'  \item{npart}{the chosen number of particles;}
#'  \item{time}{the time taken to run the routine (in seconds);}
#'  \item{\code{propVar}}{the proposal covariance for the parameter updates;}
#'  \item{\code{dataset}}{data frame containing time series count data data, of form (group, count*);}
#'  \item{\code{priors}:}{ a copy of the \code{priors} input.}
#' }
#'

PMCMC <- function(dataset, priors, iniStates, func, iniPars = NA, 
    tol = 0, fixpars = F, 
    niter = 1000, npart = 100, nprintsum = 1000, nmultskip = 1000, 
    adapt = T, propVar = NA, adaptmixprop = 0.05, nupdate = 100) {
    
    ## check inputs are present
    stopifnot(!missing(dataset) & !missing(priors) & !missing(iniStates) & !missing(func))
    
    ## check data set
    stopifnot(checkInput(dataset, "data.frame"))
    stopifnot(colnames(dataset)[1] == "time")
    stopifnot(ncol(dataset) > 1)
    stopifnot(checkInput(dataset$time, "numeric"))
    for(j in 2:ncol(dataset)) {
        stopifnot(checkInput(dataset[, j, drop = T], "numeric", int = T))
    }   
    
    ## check priors 
    stopifnot(checkInput(priors, "data.frame", ncol = 4))
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
    priors$parnames <- NULL
    priors$dist <- match(priors$dist, c("unif", "norm", "gamma"))
    priors <- as.matrix(priors)
    
    ## check function
    stopifnot(class(func) == "XPtr")
    checkXPtr(func, "SEXP", c("NumericVector", "double", "double", "int",
           "IntegerVector", "int"))
    
    ## check initial conditions
    if(!any(is.na(iniPars))) {
        stopifnot(checkInput(iniPars, c("numeric", "vector"), nrow(priors), naAllow = T))
    } else {
        iniPars <- rep(NA, nrow(priors))
    }
    stopifnot(checkInput(iniStates, c("numeric", "vector"), int = T))
    
    ## check proposal variances
    if(is.na(propVar[1])) {
        propVar <- diag(nrow(priors))
        ## adjust for scaling parameter for initial iterations
        propVar <- propVar * ((0.1 ^ 2) / nrow(propVar))
        propVar <- propVar / ((2.562 ^ 2) / nrow(propVar))
    } else {
        stopifnot(checkInput(propVar, c("numeric", "matrix"), nrow = nrow(priors), ncol = nrow(priors)))
    }
    
    ## check tolerance argument
    stopifnot(checkInput(tol, c("numeric", "vector"), 1))
    stopifnot(tol >= 0)
    
    ## check runtime arguments
    stopifnot(checkInput(fixpars, c("logical", "vector"), 1))
    if(fixpars & any(is.na(iniPars))) {
        stop("Must input initial parameters if fixing parameters.")
    }
    stopifnot(checkInput(niter, c("numeric", "vector"), 1, int = T))
    stopifnot(checkInput(npart, c("numeric", "vector"), 1, int = T))
    stopifnot(checkInput(nprintsum, c("numeric", "vector"), 1, int = T))
    stopifnot(checkInput(nmultskip, c("numeric", "vector"), 1, int = T))
    stopifnot(niter > 0 & npart > 0 & nprintsum > 0 & nmultskip > 1)
    
    ## check adaptive update and proposal covariance matrices
    stopifnot(checkInput(adapt, c("logical", "vector"), 1))
    stopifnot(checkInput(adaptmixprop, c("numeric", "vector"), 1))
    stopifnot(checkInput(nupdate, c("numeric", "vector"), 1, int = T))
    stopifnot(adaptmixprop > 0 & adaptmixprop < 1 & nupdate > 0)
    
    ## run function
    output <- PMCMC_cpp(as.matrix(dataset), priors, orig_priors$parnames, iniPars, propVar, niter, npart, 
                    adaptmixprop, tol, nprintsum, nmultskip, nupdate, as.numeric(fixpars), 
                    as.numeric(adapt), iniStates, func)
    
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
    output <- c(output[1], tol = list(tol), iniStates = iniStates, output[-1], list(dataset), list(orig_priors))
    names(output) <- c("pars", "tol", "iniStates", "skiprate", "accrate", 
        "nmultskip", "npart", "time", "propVar", "dataset", "priors")
        
    ## export class and object
    class(output) <- "PMCMC"
    output
}

