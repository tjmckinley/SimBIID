#' @title Summarises \code{PMCMC} objects
#'
#' @description Summary method for \code{PMCMC} objects.
#'
#' @param object    A \code{PMCMC} object.
#' @param transfunc Is a \code{function} object where the arguments to the function must
#'                  match all or a subset of the parameters in the model. This function needs 
#'                  to return a \code{data.frame} object with columns containing the transformed
#'                  parameters.
#' @param ...       Not used here.
#'
#' @return          A \code{summary.mcmc} object.
#'
#' @rdname summaryPMCMC
#'
#' @export

print.PMCMC <- function(x, ...) {
    ## check object is a PMCMC object
    if(class(x) != "PMCMC"){
        stop("'x' is not a PMCMC object")
    }
    
    cat("An object of class: 'PMCMC'\n")
    cat(paste0("Chain consists of ", nrow(x$pars), " iterations with ", nrow(x$priors), " parameters.\n"))
    
    ## print prior information
    temp <- x$priors %>%
        mutate(p1 = as.character(signif(p1, 2))) %>%
        mutate(p2 = as.character(signif(p2, 2))) %>%
        mutate(temp = ifelse(dist == "unif", paste0("U(lower = ", p1, ", upper = ", p2, ")"), NA)) %>%
        mutate(temp = ifelse(dist == "gamma", paste0("G(shape = ", p1, ", rate = ", p2, ")"), temp)) %>%
        mutate(temp = ifelse(dist == "norm", paste0("N(mean = ", p1, ", sd = ", p2, ")"), temp)) %>%
        mutate(temp = paste0(parnames, " ~ ", temp)) %>%
        select(temp)
    colnames(temp) <- ""
    cat("\nPriors:\n")
    print(temp, row.names = F, col.names = F, quote = F)
}

#' @rdname summaryPMCMC
#' @export

summary.PMCMC <- function(object, transfunc = NA, ...) {
    
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
        checkInput(temp, "data.frame", nrow = nrow(pars))
        temp <- as.matrix(temp)
        
        ## bind to current posterior samples
        pars <- cbind(pars, temp)
    }
    summary(as.mcmc(pars))
}
    
    
    
