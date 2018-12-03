#' @title Compile parsed Rcpp model
#'
#' @description Compiles an object of class \code{SimBIID_model} into an
#'              \code{XPtr} object for use in Rcpp functions, or an
#'              object of class \code{function} for calling directly from R.
#'
#' @export
#'
#' @param parsed: An object of class \code{SimBIID_model}.
#'
#' @return An object of class \code{XPtr} that points to the compiled function, or
#'         an R \code{function} object for calling directly from R.

compileRcpp <- function(
    parsed = NULL
) {
    if(missing(parsed)) {
        stop("'parsed' object missing")
    }
    if(class(parsed) != "SimBIID_model") {
        stop("'parsed' object not of class 'SimBIID_model'")
    }
    ## write to temporary file
    filename <- tempfile()
    writeLines(parsed$code, filename)
    
    ## compile into external pointer or function
    source(filename)
    
    ## return pointer or function
    Rcpp_object
}

