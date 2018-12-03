#' @title Summarises \code{PMCMC} objects
#'
#' @description Summary method for \code{PMCMC} objects.
#'
#' @export
#'
#' @param object    A \code{PMCMC} object.
#' @param transfunc Is a \code{function} object where the arguments to the function must
#'                  match all or a subset of the parameters in the model. This function needs 
#'                  to return a \code{data.frame} object with columns containing the transformed
#'                  parameters.
#'
#' @return          A \code{summary.mcmc} object.
#'

summary.PMCMC <- function(object, transfunc = NA) {
    
    ## check object is a PMCMC object
    if(class(object) != "PMCMC"){
        stop("'object' is not a PMCMC object")
    }
    
    pars <- as.matrix(object$pars)
    
    ## check transformations
    stopifnot(length(transfunc) == 1)
    if(is.function(transfunc)) {
    
        ## check function arguments
        fargs <- formals(transfunc)
        stopifnot(all(names(fargs) %in% colnames(pars)))
        
        ## perform transformations if required
        temppars <- pars[, match(names(fargs), colnames(pars))]
        temppars <- as.data.frame(temppars)
        temppars <- as.list(temppars)
        names(temppars) <- names(fargs)
        temp <- do.call("transfunc", temppars)
        stopifnot(checkInput(temp, "data.frame", nrow = nrow(pars)))
        temp <- as.matrix(temp)
        
        ## bind to current posterior samples
        pars <- cbind(pars, temp)
    }
    summary(as.mcmc(pars))
}
    
    
    
