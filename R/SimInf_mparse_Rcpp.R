## function to parse SimInf model and turn into Rcpp code
## (code based on that of SimInf_mparse in SimInf package)

SimInf_mparse_Rcpp <- function(
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
    ldata_names <- NULL
    transitions <- SimInf:::parse_transitions(
        transitions, compartments, ldata_names, gdata_names
    )

    ## write Rcpp code to file
    Rcpp_code <- Rcpp_mparse(transitions)
    filename <- tempfile()
    writeLines(Rcpp_code, filename)
    source(filename)
    Rcpp_ptr
}

