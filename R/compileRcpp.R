#' @title Compiles \code{SimBIID_model} object
#'
#' @description Compiles an object of class \code{SimBIID_model} into an
#'              \code{XPtr} object for use in Rcpp functions, or an
#'              object of class \code{function} for calling directly from R.
#'
#' @export
#'
#' @param model: An object of class \code{SimBIID_model}.
#'
#' @return An object of class \code{XPtr} that points to the compiled function, or
#'         an R \code{function} object for calling directly from R.

compileRcpp <- function(
    model = NULL
) {
    if(missing(model)) {
        stop("'model' object missing")
    }
    if(class(model) != "SimBIID_model") {
        stop("'model' object not of class 'SimBIID_model'")
    }
    ## write to temporary file
    filename <- tempfile()
    writeLines(model$code, filename)
    
    ## compile into external pointer or function
    source(filename)
    
    ## return pointer or function
    Rcpp_object
}

