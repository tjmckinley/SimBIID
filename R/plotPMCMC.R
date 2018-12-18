#' @title Plots \code{PMCMC} objects
#'
#' @description Plot method for \code{PMCMC} objects.
#'
#' @export
#'
#' @param x             A \code{PMCMC} object.
#' @param type          Takes the value \code{"post"} if you want to plot posterior distributions.
#'                      Takes the value \code{"trace"} if you want to plot the trace plots.
#' @param joint         A logical describing whether joint or marginal distributions are wanted.
#' @param transfunc     Is a \code{function} object where the arguments to the function must
#'                      match all or a subset of the parameters in the model. This function needs 
#'                      to return a \code{data.frame} object with columns containing the transformed
#'                      parameters.
#' @param ask           Should the user ask before moving onto next trace plot.
#'
#' @return A plot of the (approximate) posterior distributions from the particle MCMC algorithm,
#'         or corresponding trace plots.

plot.PMCMC <- function(x, type = c("post", "trace"), joint = F, transfunc = NA, ask = T) {
    
    ## check x
    if(class(x) != "PMCMC"){
        stop("'x' not a PMCMC object")
    }
    
    ## check type
    type <- type[1]
    checkInput(type, "character", inSet = c("post", "trace"))
    
    ## check joint
    checkInput(joint, c("vector", "logical"), 1)
    if(joint) {
        if(!require(GGally)) {
            stop("'GGally' package required for joint distribution plots")
        }
    }
    ## check ask
    checkInput(ask, c("vector", "logical"), 1)
    
    ## extract levels
    collev <- colnames(x$pars)
    
    if(type == "post") { 
        p <- as.matrix(x$pars) %>% as.data.frame()
                
        ## check for transformations if required
        checkInput(transfunc, length = 1, naAllow = T)
        if(is.function(transfunc)) {
        
            ## check function arguments
            fargs <- formals(transfunc)
            checkInput(names(fargs), inSet = colnames(p))
            
            ## perform transformations if required
            temppars <- p[, match(names(fargs), colnames(p))]
            temppars <- as.list(temppars)
            names(temppars) <- names(fargs)
            temp <- do.call("transfunc", temppars)
            checkInput(temp, "data.frame", nrow = nrow(p))
            
            ## bind to current posterior samples
            p <- cbind(p, temp)
            collev <- colnames(p)
        } else {
            if(!is.na(transfunc)) {
                stop("'transfunc' incorrectly specified")
            }
        }
        
        ## plot posteriors
        if(!joint) {
             p <- p %>%
                gather(Parameter, value) %>%
                mutate(Parameter = factor(Parameter, levels = collev)) %>%
                ggplot(aes(x = value)) +
                    geom_density() +
                    xlab("Parameter value") + ylab("Density") +
                    facet_wrap(~ Parameter, scales = "free")
         } else {
            p <- ggpairs(p,
                diag = list(continuous = wrap("densityDiag", alpha = 0.8)),
                lower = list(continuous = "density"),
                upper = list(continuous = "blank"))
         }
         print(p)
     } else {
        plot(x$pars, density = F, ask = ask)
     }
}
    
