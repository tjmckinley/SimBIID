#' @title Summarises \code{ABCSMC} objects
#'
#' @description Summary method for \code{ABCSMC} objects.
#'
#' @export
#'
#' @param object    An \code{ABCSMC} object.
#' @param gen       The generation of ABC that you want to extract. If left missing then
#'                  defaults to final generation.
#' @param transfunc Is a \code{function} object where the arguments to the function must
#'                  match parameters in the model. This function needs to return a \code{data.frame}
#'                  object with a single column containing the transformed parameters.
#'
#' @return          A \code{data.frame} with weighted posterior means and variances.
#'

summary.ABCSMC <- function(object, gen = NA, transfunc = NA) {
    
    ## check x is an ABCSMC object
    stopifnot(class(object) == "ABCSMC")
    
    ## check gen is valid
    stopifnot(length(gen) == 1)
    gen <- ifelse(is.na(gen), length(object$pars), gen)
    checkInput(gen, "numeric", 1, int = T)
    stopifnot(gen > 0)
    stopifnot(gen <= length(object$pars))
    
    ## extract relevant parts of the object
    weights <- object$weights[[gen]]
    pars <- object$pars[[gen]]
    pars <- as.data.frame(pars)
    
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
        stopifnot(checkInput(temp, "data.frame", nrow = nrow(pars), ncol = 1))
        
        ## bind to current posterior samples
        pars <- cbind(pars, temp)
    }
    
    ## extract parameter names
    parnames <- colnames(pars)
    
    ## calculate weighted mean
    postmn <- apply(cbind(weights, pars), 1, function(x) x[-1] * x[1])
    if(is.null(nrow(postmn))) postmn <- matrix(postmn, nrow = 1)
    postmn <- apply(postmn, 1, sum)
    
    ## calculate weighted variance
    postvar <- apply(rbind(postmn, pars), 2, function(x) (x[-1] - x[1])^2)
    postvar <- apply(cbind(weights, postvar), 1, function(x) x[-1] * x[1])
    if(is.null(nrow(postvar))) postvar <- matrix(postvar, nrow = 1)
    postvar <- apply(postvar, 1, sum)
    
    ## return summary
    postsum <- data.frame(Mean = postmn, SD = sqrt(postvar))
    rownames(postsum) <- parnames
    postsum
}
    
    
    
