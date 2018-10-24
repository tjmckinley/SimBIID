## Cpp code: Generate Rcpp code for mparse
## (based on idea in SimInf_mparse in SimInf package)
Rcpp_mparse <- function(transitions, matchCrit, addVars, stopCrit, tspan, runFromR) {
    
    ## read source code
    Rcpp_code <- readLines(system.file("", "simFunction.R", package = "ABCSMC"))
    
    ## set compilation type
    if(runFromR) {
        compType <- "Rcpp_object <- Rcpp::cppFunction"
        ## set output type
        if(is.null(matchCrit)) {
            if(is.null(tspan)){
                compType <- paste0(compType, "('NumericVector ")
            } else {
                compType <- paste0(compType, "('NumericMatrix ")
            }
        } else {
            compType <- paste0(compType, "('IntegerVector ")
        }
    } else {
        compType <- "Rcpp_object <- RcppXPtrUtils::cppXPtr('SEXP "
    }
    
    ## set matching critera
    if(!is.null(matchCrit)) {
        Rcpp_code[1] <- paste0(compType, "simFunction(NumericVector gdata, double tstart, double tstop, IntegerVector u, IntegerVector tols, IntegerVector counts, IntegerVector whichind")
    } else {
        Rcpp_code[1] <- paste0(compType, "simFunction(NumericVector gdata, double tstart, double tstop, IntegerVector u")
    }
    
    ## add additional variables to parser
    if(!is.null(addVars)) {
        Rcpp_code[1] <- paste(Rcpp_code[1], addVars)
    }
    Rcpp_code[1] <- paste0(Rcpp_code[1], ") { ")
    
    ## extract rate markers
    ratelines <- sort(c(grep("RATELINES", Rcpp_code), grep("MATCHCRIT", Rcpp_code), grep("TSPAN", Rcpp_code)))
    stopifnot(length(ratelines) == 7)
    
    ## number of transitions
    nrates <- length(transitions)
    currline <- ratelines[1]
    lines <- paste0("    NumericVector rates(", nrates, ");")
    Rcpp_code <- c(Rcpp_code[1:(currline - 1)], lines, Rcpp_code[(currline + 1):length(Rcpp_code)])
    currline <- currline + length(lines) - 1
    ratelines <- ratelines[-1] + length(lines) - 1
    
    ## update rates
    upRates <- c("", "    // update rates")
    for(i in 1:length(transitions)) {
        temp <- transitions[[i]]$propensity
        temp <- paste0("    rates[", i - 1, "] = ", temp, ";")
        upRates <- c(upRates, temp)
    }
    upRates <- c(upRates, "    totrate = sum(rates);")
    Rcpp_code <- c(Rcpp_code[1:currline], upRates, Rcpp_code[(currline + 1):length(Rcpp_code)])
    ratelines <- ratelines + length(upRates)
    
    ## update output size
    if(is.null(matchCrit)){
        ## initialise tspan outputs
        if(!is.null(tspan)){
            outsize <- paste0(paste(rep(" ", 4), collapse = ""), "NumericMatrix out(tspan.size() + 1, u.size() + 2);")
        } else {
            outsize <- paste0(paste(rep(" ", 4), collapse = ""), "NumericVector out(u.size() + 2);")
        }
    } else {
        outsize <- paste0(paste(rep(" ", 4), collapse = ""), "IntegerVector out(u.size() + 1);")
    }
    currline <- ratelines[1]
    Rcpp_code <- c(Rcpp_code[1:(currline - 1)], outsize, Rcpp_code[(currline + 1):length(Rcpp_code)])
    ratelines <- ratelines[-1]
    
    ## update tspan
    if(!is.null(tspan)) {
        tempnSpace <- 4
        tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
        upTspan <- paste0(tempSpace, "while(tspan[k] < t && k < tspan.size()) {")
        tempSpace <- paste0(rep(" ", tempnSpace + 4), collapse = "")
        upTspan <- c(upTspan, paste0(tempSpace, "out(k, 0) = NA_REAL;"))
        upTspan <- c(upTspan, paste0(tempSpace, "out(k, 1) = tspan[k];"))
        upTspan <- c(upTspan, paste0(tempSpace, "out(k, Range(2, u.size() + 1)) = as<NumericVector>(u);"))
        upTspan <- c(upTspan, paste0(tempSpace, "k++;"))
        tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
        upTspan <- c(upTspan, paste0(tempSpace, "}"))
        
        currline <- ratelines[1]
        Rcpp_code <- c(
            Rcpp_code[1:(currline - 1)], 
            paste0(tempSpace, "// update tspan"), 
            paste0(tempSpace, "k = 0;"), 
            upTspan, 
            Rcpp_code[(currline + 1):length(Rcpp_code)]
        )
        ratelines <- ratelines[-1] + length(upTspan) + 1
    } else {
        Rcpp_code <- Rcpp_code[-ratelines[1]]
        ratelines <- ratelines[-1] - 1
    }
    
    ## update states
    tempnSpace <- 12
    tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
    upStates <- paste0(tempSpace, "cumrate = rates[0];")
    for(i in 1:length(transitions)) {
        if(i < length(transitions)) {
            if(i > 1) {
                upStates <- c(upStates, paste0(tempSpace, "cumrate += rates[", i - 1, "];"))
            }
            upStates <- c(upStates, paste0(tempSpace, "if(u_tmp < cumrate) {"))
            tempnSpace <- tempnSpace + 4
            tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
        }
        ## update negative states
        temp <- which(transitions[[i]]$S < 0)
        if(length(temp) > 1){
            stop("Cannot update more than one state for each transition")
        }
        if(length(temp) > 0) {
            upStates <- c(upStates, paste0(tempSpace, "u[", temp - 1, "]--;"))
        }
        ## update positive states
        temp <- which(transitions[[i]]$S > 0)
        if(length(temp) > 1){
            stop("Cannot update more than one state for each transition")
        }
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
    currline <- ratelines[1]
    ratelines <- ratelines[-1]
    Rcpp_code <- c(Rcpp_code[1:(currline - 1)], upStates, upRates, Rcpp_code[(currline + 1):length(Rcpp_code)])
    ratelines <- ratelines + length(upStates) + length(upRates) - 1
    
    ## update tspan
    if(!is.null(tspan)) {
        currline <- ratelines[1]
        upTspan <- paste(paste(rep(" ", 8), collapse = ""), upTspan)
        Rcpp_code <- c(
            Rcpp_code[1:(currline - 1)], 
            paste0(paste(rep(" ", 12), collapse = ""), "// update tspan"),
            upTspan, 
            Rcpp_code[(currline + 1):length(Rcpp_code)]
        )
        ratelines <- ratelines[-1] + length(upTspan)
    } else {
        Rcpp_code <- Rcpp_code[-ratelines[1]]
        ratelines <- ratelines[-1] - 1
    }
    
    ## add stopping criteria
    if(!is.null(stopCrit)) {
        currline <- ratelines[1]
        Rcpp_code <- c(Rcpp_code[1:(currline - 1)], "", stopCrit, Rcpp_code[(currline + 1):length(Rcpp_code)])
        ratelines <- ratelines[-1] + length(stopCrit)
    } else {
        Rcpp_code <- Rcpp_code[-ratelines[1]]
        ratelines <- ratelines[-1] - 1
    }
    if(is.null(matchCrit)) {
        tempSpace <- paste(rep(" ", 4), collapse = "")
        if(is.null(tspan)){
            output <- paste0(tempSpace, "out[0] = (totrate == 0.0 ? 1:0);\n",
                                tempSpace, "out[1] = t;\n",
                                tempSpace, "out[Range(2, u.size() + 1)] = as<NumericVector>(u);")
        } else {
            output <- paste0(tempSpace, "out(tspan.size(), 0) = (totrate == 0.0 ? 1:0);\n",
                                tempSpace, "out(tspan.size(), 1) = t;\n",
                                tempSpace, "out(tspan.size(), Range(2, u.size() + 1)) = as<NumericVector>(u);")
        }
    }
    currline <- ratelines[1]
    Rcpp_code <- c(Rcpp_code[1:(currline - 1)], output, Rcpp_code[(currline + 1):length(Rcpp_code)])
    Rcpp_code
}

