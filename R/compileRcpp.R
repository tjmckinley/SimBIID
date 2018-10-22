#' @title Compile parsed Rcpp model
#'
#' @description Compiles an object of class \code{parsedRcpp} into an
#'              \code{XPtr} object for use in Rcpp functions, or an
#'              object of class \code{function} for calling directly from R.
#'
#' @export
#'
#' @param parsed: An object of class \code{parsedRcpp}.
#'
#' @return An object of class \code{XPtr} that points to the compiled function, or
#'         an R \code{function} object for calling directly from R.

compileRcpp <- function(
    parsed = NULL
) {
    if(missing(parsed)) {
        stop("'parsed' object missing")
    }
    if(class(parsed) != "parsedRcpp") {
        stop("'parsed' object not of class 'parsedRcpp'")
    }
    ## write to temporary file
    filename <- tempfile()
    writeLines(parsed, filename)
    
    ## compile into external pointer
    source(filename)
    
    ## return pointer
    Rcpp_object
}

