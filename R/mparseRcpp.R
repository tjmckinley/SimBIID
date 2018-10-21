#' @title Parse custom model using \code{SimInf} style markup
#'
#' @description Parse custom model using \code{SimInf} style markup and 
#' create code to compile to an \code{XPtr} object to use in \code{ABCSMC} functions.
#' Does not have full functionality of \code{mparse}. Currently only supports
#' simulations on a single node.
#'
#' @export
#'
#' @param transitions: character vector containing transitions on the form \code{"X ->
#'          ... -> Y"}. The left (right) side is the initial (final)
#'          state and the propensity is written in between the
#'          \code{->}-signs. The special symbol \code{@} is reserved for the empty
#'          set. For example, \code{transitions = c("S -> k1*S*I -> I", "I ->
#'          k2*I -> R")} expresses a SIR model.
#'
#' @param compartments: contains the names of the involved compartments, for
#'          example, \code{compartments = c("S", "I", "R")}.
#'
#' @param pars: named vector of parameters.
#'
#' @return An object of class \code{parsedRcpp} that contains code to compile
#'         into an \code{XPtr} object.

mparseRcpp <- function(
    transitions = NULL, 
    compartments = NULL,
    pars = NULL
) {
    ## Check transitions
    if (!is.atomic(transitions) || !is.character(transitions) || any(nchar(transitions) == 0)) {
        stop("'transitions' must be specified in a character vector.")
    }

    ## Check compartments
    if (!is.atomic(compartments) || !is.character(compartments) ||
        any(duplicated(compartments)) || any(nchar(compartments) == 0)) {
        stop("'compartments' must be specified in a character vector.")
    }

    ## Check pars
    pars_names <- NULL
    if (!is.null(pars)) {
        if (is.data.frame(pars)) {
            pars_names <- colnames(pars)
        } else {
            if (is.atomic(pars) && is.numeric(pars)) {
                pars_names <- names(pars)
            } else {
                stop("'pars' must either be a 'data.frame' or a 'numeric' vector.")
            }
        }

        if (is.null(pars_names) || any(duplicated(pars_names)) || any(nchar(pars_names) == 0)) {
            stop("'pars' must have non-duplicated parameter names.")
        }
    }

    if (any(duplicated(c(compartments, pars_names)))) {
        stop("'pars' and 'compartments' have names in common.")
    }

    ## Parse transitions
    transitions <- SimInf:::parse_transitions(
        transitions, compartments, NULL, pars_names, NULL
    )

    ## write Rcpp code to file
    Rcpp_code <- Rcpp_mparse(transitions)
    ## replace "gdata" with "pars"
    Rcpp_code <- gsub("gdata", "pars", Rcpp_code)
    class(Rcpp_code) <- "parsedRcpp"
    Rcpp_code
}

## print function for parsedRcpp object
#' @export
print.parsedRcpp <- function(x, ...) {
    writeLines(x)
}
