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
#' @param stopCrit: A \code{character} vector including additional stopping criteria for rejecting
#'                  simulations early. These will be inserted within \code{if(CRIT){out[0] = 0; return out;}} statements
#'                  within the underlying Rcpp code, which a return value of 0 corresponds to rejecting
#'                  the simulation. Variables in \code{CRIT} must match either those in \code{compartments}
#'                  and/or \code{addVars}.
#' 
#' @param addVars: a character vector where the names specify the additional variables to be added to the 
#'                 function call. These can be used to specify variables that can be used for 
#'                 e.g. additional stopping criteria.
#'                  
#' @param tspan: A \code{logical} determining whether to return time series counts or not.
#'                  
#' @param runFromR: \code{logical} determining whether code is to be compiled to run directly in R,
#'                  or whether to be compiled as an \code{XPtr} object for use in Rcpp.
#'
#' @return An object of class \code{parsedRcpp}, which is essentially a \code{list} 
#'         containing elements:
#'         \itemize{
#'             \item{code:}{ parsed code to compile;}
#'             \item{matchCrit:}{ copy of \code{matchCrit} argument;}
#'             \item{stopCrit:}{ copy of \code{stopCrit} argument;}
#'             \item{addVars:}{ copy of \code{addVars} argument;}
#'             \item{tspan:}{ copy of \code{tspan} argument;}
#'             \item{runFromR:}{ copy of \code{runFromR} argument.}
#'         }
#'         This can be compiled into an \code{XPtr} or \code{function} object
#'         using \code{compileRcpp()}.

mparseRcpp <- function(
    transitions = NULL, 
    compartments = NULL,
    pars = NULL,
    matchCrit = F,
    addVars = NULL,
    stopCrit = NULL,
    tspan = F,
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
    checkInput(matchCrit, "logical", 1)
    if(matchCrit) {
        tn <- paste(rep(" ", 4), collapse = "")
        tn1 <- paste(rep(" ", 8), collapse = "")
        matchCrit <- c(
            "// check whether simulation matches data",
            "out[0] = 1;",
            "for(j = 0; j < counts.size(); j++) {",
            paste0(tn, "if(fabs(u[whichind[j]] - counts[j]) <= tols[j]) {"),
            paste0(tn1, "out[0] *= 1;"),
            paste0(tn, "} else {"),
            paste0(tn1, "out[0] *= 0;"),
            paste0(tn, "}"),
            "}",
            "out[Range(1, u.size())] = u;"
        )
        matchCrit <- paste0(tn, matchCrit)
    } else {
        matchCrit <- NULL
    }
    
    ## check addVars
    if(!is.null(addVars)) {
        checkInput(addVars, "character")
        addVars <- paste(paste0("double ", addVars), collapse = ", ")
        addVars <- paste0(", ", addVars)
    }
    
    ## check stopCrit
    if(!is.null(stopCrit)) {
        checkInput(stopCrit, "character")
        tn <- paste(rep(" ", 12), collapse = "")
        tn1 <- paste(rep(" ", 16), collapse = "")
        stopCrit <- lapply(stopCrit, function(x, tn, tn1) {
            x <- paste0(tn, "if(", x, "){")
            x <- c(x, paste0(tn1, "out[0] = 0;"))
            if(is.null(matchCrit)) {
                x <- c(x, paste0(tn1, "out[1] = t;"))
                x <- c(x, paste0(tn1, "out[Range(2, u.size() + 1)] = u;"))
            } else {
                x <- c(x, paste0(tn1, "out[Range(1, u.size())] = u;"))
            }
            x <- c(x, paste0(tn1, "return out;"))
            x <- c(x, paste0(tn, "}"))
        }, tn = tn, tn1 = tn1)
        stopCrit <- do.call("c", stopCrit)
        stopCrit <- c(paste0(tn, "// early stopping criteria"), stopCrit)
    }
    
    ## check tspan
    checkInput(tspan, "logical", 1)
    if(tspan) {
        if(!is.null(matchCrit)) {
            stop("'tspan' and 'matchCrit' can't be specified together")
        }
    } else {
        tspan <- NULL
    }
    
    ## check run from R
    checkInput(runFromR, "logical", 1)

    ## Parse transitions
    transitions <- SimInf:::parse_transitions(
        transitions, compartments, NULL, pars_names, NULL
    )

    ## write Rcpp code to file
    Rcpp_code <- Rcpp_mparse(transitions, matchCrit, addVars, stopCrit, tspan, runFromR)
    ## replace "gdata" with "pars"
    Rcpp_code <- gsub("gdata", "pars", Rcpp_code)
    
    ## set up output list
    output <- list(
        code = Rcpp_code,
        matchCrit = matchCrit,
        stopCrit = stopCrit,
        addVars = addVars,
        tspan = tspan,
        runFromR = runFromR
    )
    class(output) <- "parsedRcpp"
    output
}

## print function for parsedRcpp object
#' @export
print.parsedRcpp <- function(x, ...) {
    writeLines(x$code)
}
