#' @title Parse custom model using \code{SimInf} style markup
#'
#' @description Parse custom model using \code{SimInf} style markup.
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
#' @param pars: a \code{character} vector containing the names of the parameters.
#' 
#' @param obsProcess: \code{data.frame} determining the observation process. Columns must be in the order:
#'                    \code{dataNames}, \code{dist}, \code{p1}, \code{p2}. \code{dataNames} is a \code{character}
#'                    denoting the relevant compartment name to place the observation process onto; \code{dist} 
#'                    is a \code{character} specifying the distribution of the observation process (must be one of 
#'                    \code{"unif"}, \code{"pois"} or \code{"binom"} at the current time); \code{p1} is the first parameter 
#'                    (the lower bound in the case of \code{"unif"}, the rate in the case of \code{"pois"}, or the \code{size} 
#'                    in the case of \code{"binom"}); and finally \code{p2} is the second parameter (the upper bound in the
#'                    case of \code{"unif"}, \code{NA} in the case of \code{"pois"}, and \code{prob} in the case of 
#'                    \code{"binom"}).
#' 
#' @param stopCrit: A \code{character} vector including additional stopping criteria for rejecting
#'                  simulations early. These will be inserted within \code{if(CRIT){out[0] = 0; return out;}} statements
#'                  within the underlying Rcpp code, which a return value of 0 corresponds to rejecting
#'                  the simulation. Variables in \code{CRIT} must match either those in \code{compartments}
#'                  and/or \code{addVars}.
#' 
#' @param addVars: a \code{character} vector where the names specify the additional variables to be added to the 
#'                 function call. These can be used to specify variables that can be used for 
#'                 e.g. additional stopping criteria.
#'                  
#' @param tspan: A \code{logical} determining whether to return time series counts or not.
#' 
#' @param afterTstar: A \code{character} containing code to insert after each new event time is
#'                    generated. 
#'                  
#' @param runFromR: \code{logical} determining whether code is to be compiled to run directly in R,
#'                  or whether to be compiled as an \code{XPtr} object for use in Rcpp.
#'
#' @return An object of class \code{SimBIID_model}, which is essentially a \code{list} 
#'         containing elements:
#'         \itemize{
#'             \item{code:}{ parsed code to compile;}
#'             \item{transitions:}{ copy of \code{transitions} argument;}
#'             \item{compartments:}{ copy of \code{compartments} argument;}
#'             \item{pars:}{ copy of \code{pars} argument;}
#'             \item{obsProcess:}{ copy of \code{obsProcess} argument;}
#'             \item{stopCrit:}{ copy of \code{stopCrit} argument;}
#'             \item{addVars:}{ copy of \code{addVars} argument;}
#'             \item{tspan:}{ copy of \code{tspan} argument;}
#'             \item{afterTstar:}{ copy of \code{afterTstar} argument;}
#'             \item{runFromR:}{ copy of \code{runFromR} argument.}
#'         }
#'         This can be compiled into an \code{XPtr} or \code{function} object
#'         using \code{compileRcpp()}.

mparseRcpp <- function(
    transitions = NULL, 
    compartments = NULL,
    pars = NULL,
    obsProcess = NULL,
    addVars = NULL,
    stopCrit = NULL,
    tspan = F,
    afterTstar = NULL,
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
        if (!is.atomic(pars) | !is.character(pars)) {
            stop("'pars' must be a 'character' vector of parameter names.")
        }
        if (any(nchar(pars) == 0)) {
            stop("'pars' must have non-empty parameter names.")
        }
        if (any(duplicated(pars))) {
            stop("'pars' must have non-duplicated parameter names.")
        }
    } else {
        stop("Must input 'pars' as character vector of parameter names.")
    }

    ## check element names
    if (any(duplicated(c(compartments, pars)))) {
        stop("'pars' and 'compartments' have names in common.")
    }
    
    ## check obsProcess
    obsProcess_orig <- obsProcess
    if(!is.null(obsProcess[1])){
        checkInput(obsProcess, "data.frame", ncol = 4, naAllow = T)
        checkInput(colnames(obsProcess), inSet = c("dataNames", "dist", "p1", "p2"))
        # checkInput(obsProcess$dataNames, "character", inSet = compartments)
        checkInput(obsProcess$dist, "character", inSet = c("unif", "pois", "binom"))
        checkInput(obsProcess$p1, "character")
        if(!all(is.na(obsProcess$p2))){
            checkInput(obsProcess$p2, "character", naAllow = T)
        }
        ## set up compiled column
        obsProcess$compiled <- NA
        
        for(i in 1:nrow(obsProcess)){
            ## check for missing inputs
            if(obsProcess$dist[i] == "unif" | obsProcess$dist[i] == "binom" ){
                if(any(is.na(obsProcess[i, 3:4]))){
                    stop(paste0("'obsProcess' parameters can't be missing for '", obsProcess$dist[i], "'"))
                }
            } else {
                if(is.na(obsProcess$p1[i])){
                    stop(paste0("'obsProcess' p1 can't be missing for '", obsProcess$dist[i], "'"))
                }
                if(!is.na(obsProcess$p2[i])){
                    stop(paste0("'obsProcess' p2 must be missing for '", obsProcess$dist[i], "'"))
                }
            }
            ## parse character and map to compartments and parameters
            temp <- SimInf:::parse_transitions(
                paste0(compartments[1], " -> ", obsProcess$p1[i], " -> ", compartments[1]), 
                compartments, NULL, pars, NULL)
            obsProcess$p1[i] <- temp[[1]]$propensity
            
            temp <- SimInf:::parse_transitions(
                paste0(compartments[1], " -> ", obsProcess$p2[i], " -> ", compartments[1]),
                compartments, NULL, pars, NULL)
            obsProcess$p2[i] <- temp[[1]]$propensity
            
            if(obsProcess$dist[i] == "unif" | obsProcess$dist[i] == "binom" ){
                obsProcess$compiled[i] <- paste0("out[0] += R::d", obsProcess$dist[i], 
                    "(counts[", i - 1, "], ", obsProcess$p1[i], 
                    ", ", obsProcess$p2[i], ", 1);")
            } else {
                obsProcess$compiled[i] <- paste0("out[0] += R::d", obsProcess$dist[i], 
                     "(counts[", i - 1, "], ", obsProcess$p1[i], 
                     ", 1);")
            }
        }
        ## merge observation processes
        compObsProcess <- c("out[0] = 0.0;", obsProcess$compiled, 
            "out[Range(1, u.size())] = as<NumericVector>(u);")
        compObsProcess <- paste("    ", compObsProcess)
    } else {
        compObsProcess <- NULL
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
            if(!tspan) {
                x <- c(x, paste0(tn1, "out[0] = 0;"))
                if(is.null(obsProcess[1])) {
                    x <- c(x, paste0(tn1, "out[1] = t;"))
                    x <- c(x, paste0(tn1, "out[Range(2, u.size() + 1)] = u;"))
                } else {
                    x <- c(x, paste0(tn1, "out[Range(1, u.size())] = u;"))
                }
            } else {
                x <- c(x, paste0(tn1, "out(tspan.size(), 0) = 0;"))
                # if(is.null(obsProcess)) {
                    x <- c(x, paste0(tn1, "out(tspan.size(), 1) = t;"))
                    x <- c(x, paste0(tn1, "for(int m = 0; m < u.size(); m++){"))
                    x <- c(x, paste0(tn1, "    out(tspan.size(), m + 2) = u[m];"))
                    x <- c(x, paste0(tn1, "}"))
                # } else {
                #     x <- c(x, paste0(tn1, "out[tspan.size(), Range(1, u.size())] = u;"))
                # }
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
        if(!is.null(obsProcess[1])) {
            stop("'tspan' and 'obsProcess' can't be specified together")
        }
    }
    
    ## check afterTstar
    if(!is.null(afterTstar)) {
        checkInput(afterTstar, "character", 1)
        warning("No consistency check on 'afterTstar': you might want to check parsed code before compiling.")
    }
    
    ## check run from R
    checkInput(runFromR, "logical", 1)

    ## Parse transitions
    transitions1 <- SimInf:::parse_transitions(
        transitions, compartments, NULL, pars, NULL
    )

    ## write Rcpp code to file
    Rcpp_code <- Rcpp_mparse(transitions1, compObsProcess, addVars, stopCrit, tspan, afterTstar, runFromR)
    
    ## replace "gdata" with "pars"
    Rcpp_code <- gsub("gdata", "pars", Rcpp_code)
    
    ## set up output list
    output <- list(
        code = Rcpp_code,
        transitions = transitions,
        compartments = compartments,
        pars = pars,
        obsProcess = obsProcess_orig,
        stopCrit = stopCrit,
        addVars = addVars,
        tspan = tspan,
        afterTstar = afterTstar,
        runFromR = runFromR
    )
    class(output) <- "SimBIID_model"
    output
}

## print function for SimBIID_model object
#' @export
print.SimBIID_model <- function(x, ...) {
    writeLines(x$code)
}
