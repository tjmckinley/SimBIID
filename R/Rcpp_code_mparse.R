## Cpp code: Generate Rcpp code for mparse
## (based on idea in SimInf_mparse in SimInf package)
Rcpp_code_mparse <- function(transitions) {
    Rcpp_code <- readLines(system.file("", "simFunction.R", package = "ABCSMC"))
    
    ## number of transitions
    nrates <- length(transitions)
    lines <- paste0("    NumericVector rates(", nrates, ");")
    lines <- c(lines, paste0("    NumericVector cumrates(", nrates, ");"), "    ")
    currline <- 9
    Rcpp_code <- c(Rcpp_code[1:currline], lines, Rcpp_code[(currline + 1):length(Rcpp_code)])
    currline <- currline + length(lines)
    
    ## update rates
    upRates <- "    // update rates"
    for(i in 1:length(transitions)) {
        temp <- transitions[[i]]$propensity
        temp <- paste0("    rates[", i - 1, "] = ", temp, ";")
        upRates <- c(upRates, temp)
    }
    upRates <- c(upRates, "    cumrates[0] = rates[0];")
    upRates <- c(upRates, paste0("    for(j = 1; j < ", nrates, "; j++) {"))
    upRates <- c(upRates, "        cumrates[j] = cumrates[j - 1] + rates[j];")
    upRates <- c(upRates, "    }")
    Rcpp_code <- c(Rcpp_code[1:currline], upRates, Rcpp_code[(currline + 1):length(Rcpp_code)])
    currline <- currline + length(upRates)
    
    ## update states
    upStates <- character()
    tempnSpace <- 12
    tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
    for(i in 1:length(transitions)) {
        if(i < length(transitions)) {
            upStates <- c(upStates, paste0(tempSpace, "if(u_tmp < cumrates[", i - 1, "]) {"))
            tempnSpace <- tempnSpace + 4
            tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
        }
        ## update negative states
        temp <- which(transitions[[i]]$S < 0)
        stopifnot(length(temp) <= 1)
        if(length(temp) > 0) {
            upStates <- c(upStates, paste0(tempSpace, "u[", temp - 1, "]--;"))
        }
        ## update positive states
        temp <- which(transitions[[i]]$S > 0)
        stopifnot(length(temp) <= 1)
        if(length(temp) > 0) {
            upStates <- c(upStates, paste0(tempSpace, "u[", temp - 1, "]++;"))
        }
        tempnSpace <- tempnSpace - 4
        tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
        if(i < length(transitions)) {
            upStates <- c(upStates, paste0(tempSpace, "} else {"))
            tempnSpace <- tempnSpace + 4
            tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
        }
    }
    for(i in 1:(length(transitions) - 1)) {
        upStates <- c(upStates, paste0(tempSpace, "}"))
        tempnSpace <- tempnSpace - 4
        tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
    }
    
    ## update rates
    upRates <- sapply(as.list(upRates), function(x) {
        paste0("        ", x)
    })
    currline <- 18 + currline - 9
    Rcpp_code <- c(Rcpp_code[1:currline], upStates, upRates, Rcpp_code[(currline + 1):length(Rcpp_code)])
    Rcpp_code
}

