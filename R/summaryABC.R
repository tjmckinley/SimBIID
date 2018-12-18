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
#'                  match all or a subset of the parameters in the model. This function needs 
#'                  to return a \code{data.frame} object with columns containing the transformed
#'                  parameters.
#'
#' @return          A \code{data.frame} with weighted posterior means and variances.
#'
#' @rdname summaryABCSMC

print.ABCSMC <- function(x, ...) {
    ## check object is a ABCSMC object
    if(class(x) != "ABCSMC"){
        stop("'x' is not a ABCSMC object")
    }
    
    cat("An object of class: 'ABCSMC'\n")
    cat(paste0("Consists of ", nrow(x$tols), " generations with ", nrow(x$priors), " parameters.\n"))
    
    ## print data information
    cat("\nData:\n\n")
    print(x$data, row.names = F)
    
    ## print tolerance information
    temp <- x$tols %>%
        as.data.frame() %>%
        mutate(Generation = 1:n()) %>%
        select(Generation, everything())
    cat("\n\nTolerances:\n\n")
    print(temp, row.names = F)
    
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
    cat("\n\nPriors:\n")
    print(temp, row.names = F, col.names = F, quote = F)
}

#' @rdname summaryABCSMC
#' @export

summary.ABCSMC <- function(object, gen = NA, transfunc = NA) {
    
    ## check x is an ABCSMC object
    if(class(object) != "ABCSMC"){
        stop("'object' not of type ABCSMC")
    }
    
    ## check gen is valid
    if(length(gen) != 1){
        stop("'gen' must be of length 1")
    }
    gen <- ifelse(is.na(gen), length(object$pars), gen)
    checkInput(gen, "numeric", 1, int = T, gt = 0, lte = length(object$pars))
    
    ## extract relevant parts of the object
    weights <- object$weights[[gen]]
    pars <- object$pars[[gen]]
    pars <- as.data.frame(pars)
    
    ## check transformations
    if(length(transfunc) != 1){
        stop("'transfunc' must be of length 1")
    }
    if(is.function(transfunc)) {
    
        ## check function arguments
        fargs <- formals(transfunc)
        checkInput(names(fargs), inSet = colnames(pars))
        
        ## perform transformations if required
        temppars <- pars[, match(names(fargs), colnames(pars))]
        temppars <- as.data.frame(temppars)
        temppars <- as.list(temppars)
        names(temppars) <- names(fargs)
        temp <- do.call("transfunc", temppars)
        checkInput(temp, "data.frame", nrow = nrow(pars))
        
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
    
    
    
