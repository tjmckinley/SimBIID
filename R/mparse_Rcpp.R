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
#' @param gdata: optional data that are common to all nodes in the model. Can
#'          be specified either as a named numeric vector or as as a
#'          one-row data.frame. The names are used to identify the
#'          parameters in the transitions. The global data vector is
#'          passed as an argument to the transition rate functions and
#'          the post time step function.
#'
#' @return An object of class \code{parsedRcpp} that contains code to compile
#'         into an \code{XPtr} object.

mparse_Rcpp <- function(
    transitions = NULL, 
    compartments = NULL,
    gdata = NULL
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

    ## Check gdata
    gdata_names <- NULL
    if (!is.null(gdata)) {
        if (is.data.frame(gdata)) {
            gdata_names <- colnames(gdata)
        } else {
            if (is.atomic(gdata) && is.numeric(gdata)) {
                gdata_names <- names(gdata)
            } else {
                stop("'gdata' must either be a 'data.frame' or a 'numeric' vector.")
            }
        }

        if (is.null(gdata_names) || any(duplicated(gdata_names)) || any(nchar(gdata_names) == 0)) {
            stop("'gdata' must have non-duplicated parameter names.")
        }
    }

    if (any(duplicated(c(compartments, gdata_names)))) {
        stop("'gdata' and 'compartments' have names in common.")
    }

    ## Parse transitions
    transitions <- SimInf:::parse_transitions(
        transitions, compartments, NULL, gdata_names, NULL
    )

    ## write Rcpp code to file
    Rcpp_code <- Rcpp_mparse(transitions)
    class(Rcpp_code) <- "parsedRcpp"
    Rcpp_code
}

