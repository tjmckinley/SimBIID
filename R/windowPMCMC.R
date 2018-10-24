#' @title Time windows for \code{PMCMC} objects
#' 
#' @description \code{window} method for class \code{PMCMC}
#' 
#' @details Acts as a wrapper function for \code{\link[coda]{window.mcmc}} 
#'          from the \code{coda} package
#' 
#' @param x a \code{PMCMC} object, usually as a result of a call to
#'          \code{PMCMC}.
#' @param \dots arguments to pass to \code{\link{window.mcmc}}
#' @return a \code{PMCMC} object
#'
#' @export

window.PMCMC <- function(x, ...) {
    if(class(x) != "PMCMC"){
        stop("'x' is not a PMCMC object")
    }
    
    #extract 'mcmc' object
    y <- x$pars
    
    #extract subset
    y <- window(y, ...)
    
    #generate new 'bayesLog' x
    x$pars <- y
    x
}
