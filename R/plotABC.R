#' @title Plots \code{ABCSMC} objects
#'
#' @description Plot method for \code{ABCSMC} objects.
#'
#' @export
#'
#' @param x     An \code{ABCSMC} object.
#' @param type  Takes the value \code{"post"} if you want to plot posterior distributions.
#'              Takes the value \code{"output"} if you want to plot the simulated outputs.
#' @param gen   A vector of generations to plot. If left missing then defaults to all generations.
#' @param joint A logical describing whether joint or marginal distributions are wanted.
#'
#' @return A plot of the ABC posterior distributions for different generations, or the distributions
#'         of the simulated summary measures for different generations.

plot.ABCSMC <- function(x, type = c("post", "output"), gen = NA, joint = F) {
    
    ## check x
    stopifnot(class(x) == "ABCSMC")
    
    ## check type
    type <- type[1]
    stopifnot(is.character(type))
    stopifnot(type %in% c("post", "output"))
    
    ## check gens
    if(is.na(gen[1])) {
        gen <- 1:length(x$pars)
    }
    checkInput(gen, c("vector", "numeric"), int = T)
    stopifnot(all(gen %in% 1:length(x$pars)))
    gen <- as.character(sort(gen))
    
    ## check joint
    checkInput(joint, c("vector", "logical"), 1)
    if(joint) {
        if(!require(GGally)) {
            stop("'GGally' package required for joint distribution plots")
        }
    }
    
    if(type == "post") {
        ## generate colorRamp
        getPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))
        fillCols <- getPalette(length(x$pars))
        fillCols <- fillCols[as.numeric(gen)]
        
        ## get parameter names
        pnames <- priors$parnames
        
        p <- x$pars %>%
                map(~{
                    as_tibble(.) %>%
                    set_names(pnames)
                }, pnames = pnames) %>%
                bind_rows(.id = "Generation") %>%
                filter(Generation %in% gen) %>%
                mutate(Generation = factor(Generation, levels = gen))
        
        ## plot posteriors
        if(!joint) {
             p <- p %>%
                gather(Parameter, value, -Generation) %>%
                ggplot(aes(x = value, fill = Generation)) +
                    geom_density(alpha = 0.8) +
                    xlab("Parameter value") + ylab("Density") +
                    scale_fill_manual(values = fillCols) +
                    facet_wrap(~ Parameter, scales = "free")
         } else {
            p <- ggpairs(p, mapping = aes(color = Generation), 
                columns = 2:ncol(p),
                diag = list(continuous = wrap("densityDiag", alpha = 0.8)),
                lower = list(continuous = "density"),
                upper = list(continuous = "blank"),
                legend = c(1, 1))
            ## amend colours
            for(i in 1:p$nrow) {
                for(j in 1:p$ncol){
                    p[i, j] <- p[i, j] + scale_fill_manual(values = fillCols)
                    p[i, j] <- p[i, j] + scale_colour_manual(values = fillCols)
                }
            }
         }
     } else {
        ## generate colorRamp
        getPalette <- colorRampPalette(brewer.pal(9, "YlGnBu"))
        fillCols <- getPalette(length(x$pars))
        fillCols <- fillCols[as.numeric(gen)]
        
        ## get output names
        onames <- colnames(x$output[[1]])
        
        ## plot outputs
        p <- x$output %>%
            map(~{
                as_tibble(.) %>%
                set_names(onames)
            }, onames = onames) %>%
            bind_rows(.id = "Generation") %>%
            filter(Generation %in% gen) %>%
            mutate(Generation = factor(Generation, levels = gen))
            
        if(!joint) {
            p <- p %>%
                gather(Output, value, -Generation) %>%
                ggplot(aes(x = value, fill = Generation)) +
                    geom_density(alpha = 0.8) +
                    xlab("Output") + ylab("Density") +
                    scale_fill_manual(values = fillCols) +
                    facet_wrap(~ Output, scales = "free")
        } else {
            p <- ggpairs(p, mapping = aes(color = Generation), 
                columns = 2:ncol(p),
                diag = list(continuous = wrap("densityDiag", alpha = 0.8)),
                lower = list(continuous = "density"),
                upper = list(continuous = "blank"),
                legend = c(1, 1))
            ## amend colours
            for(i in 1:p$nrow) {
                for(j in 1:p$ncol){
                    p[i, j] <- p[i, j] + scale_fill_manual(values = fillCols)
                    p[i, j] <- p[i, j] + scale_colour_manual(values = fillCols)
                }
            }
        }
            
     }
     print(p)
}
    
