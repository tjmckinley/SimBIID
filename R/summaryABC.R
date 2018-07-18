#' @title Summarises \code{ABCSMC} objects
#'
#' @description Summary method for \code{ABCSMC} objects.
#'
#' @export
#'
#' @param object    An \code{ABCSMC} object.
#' @param gen       The generation of ABC that you want to extract. If left missing then
#'                  defaults to final generation.
#'
#' @return          A \code{data.frame} with weighted posterior means and variances.
#'

summary.ABCSMC <- function(object, gen = NA) {
    
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
    rownames(postsum) <- priors$parnames
    postsum
}
    
    
    
