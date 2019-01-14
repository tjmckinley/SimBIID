#' @title Predicts future course of outbreak from \code{PMCMC} objects
#'
#' @description Predict method for \code{PMCMC} objects.
#'
#' @param object             A \code{PMCMC} object.
#' @param tspan         A vector of times over which to output predictions.
#' @param npart         The number of particles to use in the bootstrap filter.
#' @param ...           Not used here.
#'
#' @return A \code{SimBIID_runs} object.
#'         
#' @method predict PMCMC
#' @export 

predict.PMCMC <- function(object, tspan, npart = 50, ...) {
    
    ## check object
    if(class(object) != "PMCMC"){
        stop("'object' not a PMCMC object")
    }
    
    if(class(object$func) != "SimBIID_model"){
        stop("'object$func' not a SimBIID_model object")
    }
    
    if(missing(tspan)){
        stop("'tspan' can't be missing")
    }
    
    ## check stop
    checkInput(tspan, c("vector", "numeric"), gt = max(object$data$t))
    tspan <- sort(tspan)
    
    ## check npart
    checkInput(npart, c("vector", "numeric"), 1, int = T, gt = 0)
    
    ## extract parameters and remove extraneous columns
    pars <- as.matrix(object$pars) %>%
        as.data.frame() %>%
        select_(.dots = object$func$pars) %>%
        as.matrix()
    
    ## run bootstrap filter to get states at each time point of the dataset
    prevStates <- bootStates(object$data, object$func, pars, object$u, npart)
    
    ## extract final states
    iniStates <- prevStates %>%
        group_by(rep) %>%
        slice(n()) %>%
        ungroup() %>%
        select(one_of(object$func$compartments)) %>%
        as.matrix()
    
    ## generate model
    func <- object$func
    if(!func$tspan){
        cat("For predictions 'SimBIID_model' object will have 'tspan' set to T\n")
    }
    # if(!is.null(func$obsProcess[1])){
    #     cat("For predictions, 'obsProcess' will be removed from 'SimBIID_model' simulations\n")
    # }
    if(!is.null(func$addVars[1])){
        stop("'SimBIID_model' can't have non-NULL 'addVars'")
    }
    if(!is.null(func$stopCrit[1])){
        stop("'SimBIID_model' can't have non-NULL 'stopCrit'")
    }
    func <- mparseRcpp(
        transitions = func$transitions,
        compartments = func$compartments,
        pars = func$pars,
        obsProcess = func$obsProcess,
        addVars = NULL,
        stopCrit = NULL,
        tspan = T,
        afterTstar = NULL,
        PF = F,
        runFromR = T
    )
    compfunc <- compileRcpp(func)
    
    ## use run method to produce forward predictions
    out <- list()
    outsums <- list()
    for(i in 1:nrow(pars)) {
        out[[i]] <- compfunc(pars[i, ], max(object$data$t), max(tspan), iniStates[i, ], tspan)
        outsums[[i]] <- out[[i]][[1]]
        out[[i]] <- out[[i]][[2]]
        outsums[[i]] <- as.data.frame(matrix(outsums[[i]], nrow = 1))
        out[[i]] <- as.data.frame(out[[i]])
        tempnames <- c("completed", "t", colnames(iniStates))
        if(is.data.frame(func$obsProcess)) {
            tempnames <- c(tempnames, func$obsProcess$dataNames)
        }
        colnames(outsums[[i]]) <- tempnames
        tempnames <- c("t", colnames(iniStates))
        if(is.data.frame(func$obsProcess)) {
            tempnames <- c(tempnames, func$obsProcess$dataNames)
        }
        colnames(out[[i]]) <- tempnames
    }
    ## bind to prevStates
    out <- out %>%
        bind_rows(.id = "rep") %>%
        rbind(prevStates) %>%
        arrange(rep, t)
    
    ## extract summaries
    outsums <- outsums %>%
        bind_rows(.id = "rep")
    
    ## bind as output list and return
    out <- list(sums = outsums, runs = out, bootEnd = max(object$data$t))
    class(out) <- "SimBIID_runs"
    out
}
    
