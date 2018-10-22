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
#' @param matchCrit: \code{logical} determining whether to implement match criteria or not.
#' 
#' @param addVars: a named vector where the names specify the additional variables and the
#'                 values specify the values of the variable. These can be used to specify variables 
#'                 that can be used for stopping criteria.
#' 
#' @param stopCrit: A \code{character} vector including additional stopping criteria for rejecting
#'                  simulations early. These will be inserted within \code{if(CRIT){out[0] = 0; return out;}} statements
#'                  within the underlying Rcpp code, which a return value of 0 corresponds to rejecting
#'                  the simulation. Variables in \code{CRIT} must match either those in \code{compartments}
#'                  and/or \code{addVars}.
#'                  
#' @param runFromR: \code{logical} determining whether code is to be compiled to run directly in R,
#'                  or whether to be compiled as an \code{XPtr} object for use in Rcpp.
#'
#' @return An object of class \code{parsedRcpp} that contains code to compile
#'         into an \code{XPtr} object.

mparseRcpp <- function(
    transitions = NULL, 
    compartments = NULL,
    pars = NULL,
    matchCrit = F,
    addVars = NULL,
    stopCrit = NULL,
    runFromR = T
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

    ## check element names
    if (any(duplicated(c(compartments, pars_names)))) {
        stop("'pars' and 'compartments' have names in common.")
    }
    
    ## check matchCrit
    stopifnot(checkInput(matchCrit, "logical", 1))
    if(matchCrit) {
        tn <- paste(rep(" ", 4), collapse = "")
        tn1 <- paste(rep(" ", 8), collapse = "")
        matchCrit <- c(
            "// check whether simulation matches data",
            "out[0] = 1;",
            "for(j = 0; j < counts.size(); j++) {",
            paste0(tn, "if(fabs(uNew[whichind[j]] - counts[j]) <= tols[j]) {"),
            paste0(tn1, "out[0] *= 1;"),
            paste0(tn, "} else {"),
            paste0(tn1, "out[0] *= 0;"),
            paste0(tn, "}"),
            "}",
            "out[Range(1, u.size())] = uNew;"
        )
        matchCrit <- paste0(tn, matchCrit)
    } else {
        matchCrit <- NULL
    }
    
    ## check addVars
    if(!is.null(addVars)) {
        stopifnot(checkInput(addVars, "character"))
        addVars <- paste(paste0("double ", addVars), collapse = ", ")
        addVars <- paste0(", ", addVars)
    }
    
    ## check stopCrit
    if(!is.null(stopCrit)) {
        stopifnot(checkInput(stopCrit, "character"))
        tn <- paste(rep(" ", 12), collapse = "")
        tn1 <- paste(rep(" ", 16), collapse = "")
        stopCrit <- lapply(stopCrit, function(x, tn, tn1) {
            x <- paste0(tn, "if(", x, "){")
            x <- c(x, paste0(tn1, "out[0] = 0;"))
            x <- c(x, paste0(tn1, "return out;"))
            x <- c(x, paste0(tn, "}"))
        }, tn = tn, tn1 = tn1)
        stopCrit <- do.call("c", stopCrit)
        stopCrit <- c(paste0(tn, "// early stopping criteria"), stopCrit)
    }
    
    ## check run from R
    stopifnot(checkInput(runFromR, "logical", 1))

    ## Parse transitions
    transitions <- SimInf:::parse_transitions(
        transitions, compartments, NULL, pars_names, NULL
    )

    ## write Rcpp code to file
    Rcpp_code <- Rcpp_mparse(transitions, matchCrit, addVars, stopCrit, runFromR)
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
