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
    checkInput(tspan, c("vector", "numeric"), gt = max(object$data$time))
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
    
    ## convert to data frame
    prevStates <- prevStates %>%
        as.data.frame()
    
    ## extract final states
    iniStates <- prevStates %>%
        group_by(rep) %>%
        slice(n()) %>%
        ungroup() %>%
        select(-rep, -t) %>%
        as.matrix()
    
    ## generate model
    func <- object$func
    if(!func$tspan){
        warning("For predictions 'SimBIID_model' object will have 'tspan' set to T")
    }
    if(!is.null(func$obsProcess[1])){
        warning("For predictions, 'SimBIID_model' must have NULL 'obsProcess'")
    }
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
        obsProcess = NULL,
        addVars = NULL,
        stopCrit = NULL,
        tspan = T,
        afterTstar = NULL,
        runFromR = T
    )
    func <- compileRcpp(func)
    
    ## use run method to produce forward predictions
    out <- list()
    outsums <- list()
    for(i in 1:nrow(pars)) {
        out[[i]] <- func(pars[i, ], max(object$data$time), max(tspan), iniStates[i, ], tspan)
        outsums[[i]] <- out[[i]][[1]]
        out[[i]] <- out[[i]][[2]]
        outsums[[i]] <- as.data.frame(matrix(outsums[[i]], nrow = 1))
        out[[i]] <- as.data.frame(out[[i]])
        colnames(outsums[[i]]) <- c("completed", "t", colnames(iniStates))
        colnames(out[[i]]) <- c("t", colnames(iniStates))
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
    out <- list(sums = outsums, runs = out)
    class(out) <- "SimBIID_runs"
    out
}
    
