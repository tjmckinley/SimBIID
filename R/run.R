#' @title Runs \code{SimBIID_model} object
#'
#' @description Wrapper function that compiles (if necessary) and runs
#'              a \code{SimBIID_model} object. Returns results in a 
#'              user-friendly manner as a \code{SimBIID_run} object, 
#'              for which \code{summary()} and \code{plot()} generics 
#'              are provided.
#'
#' @export
#'
#' @param model: An object of class \code{SimBIID_model}.
#' @param pars: A named vector of parameters.
#' @param tstart: The time at which to start the simulation.
#' @param tstop: The time at which to stop the simulation.
#' @param u: A vector of initial states.
#' @param tspan: A numeric vector containing the times at which to 
#'               save the states of the system.
#' @param nrep: Specifies the number of simulations to run.
#' @param parallel  A \code{logical} determining whether to use parallel processing or not.
#' @param mc.cores  Number of cores to use if using parallel processing.
#'
#' @return An object of class \code{SimBIID_run}, essentially a list 
#'         containing elements:
#'         \itemize{
#'             \item{sums:}{ a \code{data.frame()} with summaries of the model runs. This
#'             includes columns \code{run}, \code{completed}, \code{t}, \code{u*} 
#'             (see help file for \code{SimBIID_model} for more details);}
#'             \item{runs:}{ a \code{data.frame()} object, containing columns: \code{run},
#'             \code{t}, \code{u*} (see help file for \code{SimBIID_model} for more details).
#'             These contain time series counts for the simulations. Note that this will
#'             only be returned if \code{tspan = T} in the original \code{SimBIID_model} object.}
#'         } 
#' 

run <- function(
    model,
    pars,
    tstart,
    tstop,
    u,
    tspan,
    nrep = 1,
    parallel = F,
    mc.cores = NA
) {
    ## check inputs
    if(missing(model)) {
        stop("'model' object missing")
    }
    if(class(model) != "SimBIID_model") {
        stop("'model' object not of class 'SimBIID_model'")
    }
    if(missing(pars)) {
        stop("'pars' object missing")
    }
    if(missing(tstart)) {
        stop("'tstart' object missing")
    }
    if(missing(tstop)) {
        stop("'tstop' object missing")
    }
    if(missing(u)) {
        stop("'u' object missing")
    }
    ## check inputs
    if(!model$runFromR) {
        stop("'model' must be specified with 'runFromR = T'")
    }
    checkInput(pars, c("vector", "numeric"), length(model$pars))
    parnames <- names(pars)
    checkInput(parnames, inSet = model$pars)
    pars <- pars[match(model$pars, parnames)]
    checkInput(tstart, c("vector", "numeric"), 1)
    checkInput(tstop, c("vector", "numeric"), 1, gt = tstart)
    checkInput(u, c("vector", "numeric"), length(model$compartments), gte = 0, int = T)
    unames <- names(u)
    checkInput(unames, c("vector", "character"), length(model$compartments))
    checkInput(unames, inSet = model$compartments)
    u <- u[match(model$compartments, unames)]
    if(missing(tspan) & model$tspan){
        stop("'SimBIID_model' requires that 'tspan' must be set")
    }
    if(!missing(tspan) & !model$tspan){
        stop("'SimBIID_model' does not allow a 'tspan' argument")
    }
    if(!missing(tspan)){
        checkInput(tspan, c("vector", "numeric"), int = T)
    }
    checkInput(nrep, c("numeric"), 1, int = T, gt = 0)
    checkInput(parallel, c("vector", "logical"), 1)
    if(parallel) {
        if(!require(parallel)) {
            stop("Must have 'parallel' package installed to use parallelisation")
        }
        nc <- detectCores()
        nc <- ifelse(is.na(nc), 1, nc)
        if(!is.na(mc.cores[1])) {
            checkInput(mc.cores, "numeric", 1, int = T)
            mc.cores <- min(nc, mc.cores)
        } else {
            mc.cores <- nc
        }
        parallel <- (mc.cores > 1)
        cat(paste0("Number of cores: ", mc.cores, "\n"))
    }
    
    ## compile model
    compModel <- compileRcpp(model)
    
    ## run simulations
    if(!parallel | nrep == 1) {
        if(!missing(tspan)) {
            sims <- lapply(1:nrep, function(i, model, pars, tstart, tstop, u, tspan){
                model(pars, tstart, tstop, u, tspan)
            }, model = compModel, pars = pars, tstart = tstart, tstop = tstop, u = u, tspan = tspan)
        } else {
            sims <- lapply(1:nrep, function(i, model, pars, tstart, tstop, u){
                model(pars, tstart, tstop, u)
            }, model = compModel, pars = pars, tstart = tstart, tstop = tstop, u = u)
        }
    } else {
        if(!missing(tspan)) {
            sims <- mclapply(1:nrep, function(i, model, pars, tstart, tstop, u, tspan){
                model(pars, tstart, tstop, u, tspan)
            }, model = compModel, pars = pars, tstart = tstart, tstop = tstop, u = u, tspan = tspan, mc.cores = mc.cores)
        } else {
            sims <- mclapply(1:nrep, function(i, model, pars, tstart, tstop, u){
                model(pars, tstart, tstop, u)
            }, model = compModel, pars = pars, tstart = tstart, tstop = tstop, u = u, mc.cores = mc.cores)
        }
    }
    ## extract simulations in useful format
    if(missing(tspan)) {
        sims <- do.call("rbind", sims)
        sims <- cbind(1:nrow(sims), sims)
        colnames(sims) <- c("rep", "completed", "t", model$compartments)
        sums <- as.data.frame(sims)
        runs <- NA
    } else {
        sums <- do.call("rbind", lapply(sims, function(x){
            x[nrow(x), , drop = F]
        }))
        sums <- rbind(1:nrow(sums), sums)
        colnames(sums) <- c("rep", "completed", "t", model$compartments)
        sums <- as.data.frame(sums)
        
        sims <- do.call("rbind", lapply(1:length(sims), function(i, x){
            x <- x[[i]]
            x <- x[-nrow(x), , drop = F]
            x <- cbind(rep(i, nrow(x)), x)
            x
        }, x = sims))
        sims <- do.call("rbind", sims)
        colnames(sims) <- c("rep", "t", model$compartments)
        runs <- as.data.frame(sims)
    }
    ## return list of simulated outputs
    out <- list(sums = sums, runs = runs)
    class(out) <- "SimBIID_runs"
    out
}

