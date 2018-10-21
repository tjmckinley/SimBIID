#' @title Compile parsed Rcpp model
#'
#' @description Compiles an object of class \code{parsedRcpp} into an
#'              \code{XPtr} object for use in \code{ABCSMC} functions.
#'
#' @export
#'
#' @param parsed: An object of class \code{parsedRcpp}.
#'
#' @return An object of class \code{XPtr} that points to the compiled function.

mparse_Rcpp <- function(
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
    writeLines(Rcpp_code, filename)
    
    ## compile into external pointer
    source(filename)
    
    ## return pointer
    Rcpp_ptr
}

