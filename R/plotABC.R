#' @title Plots \code{ABCSMC} objects
#'
#' @description Plot method for \code{ABCSMC} objects.
#'
#' @param x     An \code{ABCSMC} object.
#' @param type  Takes the value \code{"post"} if you want to plot posterior distributions.
#'              Takes the value \code{"output"} if you want to plot the simulated outputs.
#' @param gen   A vector of generations to plot. If left missing then defaults to all generations.
#' @param joint A logical describing whether joint or marginal distributions are wanted.
#' @param transfunc Is a \code{function} object where the arguments to the function must
#'                  match all or a subset of the parameters in the model. This function needs 
#'                  to return a \code{data.frame} object with columns containing the transformed
#'                  parameters.
#' @param ... Not used here.
#'
#' @return A plot of the ABC posterior distributions for different generations, or the distributions
#'         of the simulated summary measures for different generations.
#' 
#' @method plot ABCSMC        
#' @export

plot.ABCSMC <- function(x, type = c("post", "output"), gen = NA, joint = F, transfunc = NA, ...) {
    
    ## check x
    if(class(x) != "ABCSMC"){
        stop("'x' is not a ABCSMC object")
    }
    
    ## check type
    type <- type[1]
    checkInput(type, "character", inSet = c("post", "output"))
    
    ## check gens
    if(is.na(gen[1])) {
        gen <- 1:length(x$pars)
    }
    checkInput(gen, c("vector", "numeric"), int = T, inSet = 1:length(x$pars))
    gen <- as.character(sort(gen))
    
    ## check joint
    checkInput(joint, c("vector", "logical"), 1)
    if(joint) {
        if(!requireNamespace("GGally", quietly = TRUE)) {
            stop("'GGally' package required for joint distribution plots")
        }
    }
    
    if(type == "post") {
        ## generate colorRamp
        getPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))
        fillCols <- getPalette(length(x$pars))
        fillCols <- fillCols[as.numeric(gen)]
        
        ## get parameter names
        pnames <- x$priors$parnames
        
        p <- x$pars %>%
                map(~{
                    as_tibble(.) %>%
                    set_names(pnames)
                }, pnames = pnames) %>%
                bind_rows(.id = "Generation") %>%
                dplyr::filter(Generation %in% gen) %>%
                mutate(Generation = factor(Generation, levels = gen))
                
        ## check for transformations if required
        if(length(transfunc) != 1){
            stop("'transfunc' not of length 1")
        }
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
        } else {
            if(!is.na(transfunc)) {
                stop("'transfunc' incorrectly specified")
            }
        }
        
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
            p <- GGally::ggpairs(p, mapping = aes(colour = Generation), 
                columns = 2:ncol(p),
                diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.8)),
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
            dplyr::filter(Generation %in% gen) %>%
            mutate(Generation = factor(Generation, levels = gen))
            
        if(length(transfunc) != 1){
            stop("'transfunc' not of length 1")
        }
        if(is.function(transfunc)) {
            stop("'transfunc' can't be used for output plots")
        } else {
            if(!is.na(transfunc)) {
                stop("'transfunc' incorrectly specified")
            }
        }
        
        ## extract data
        dat <- x$data %>% as.list() %>% data.frame() %>%
            gather(Output, value)
            
        if(!joint) {
            p <- p %>%
                gather(Output, value, -Generation) %>%
                ggplot(aes(x = value, fill = Generation)) +
                    geom_density(alpha = 0.8) +
                    xlab("Output") + ylab("Density") +
                    scale_fill_manual(values = fillCols) +
                    facet_wrap(~ Output, scales = "free") +
                    geom_vline(aes(xintercept = value), data = dat, linetype = 2, colour = "red", size = 1.5)
        } else {
            p <- GGally::ggpairs(p, mapping = aes(colour = Generation), 
                columns = 2:ncol(p),
                diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.8)),
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
    
