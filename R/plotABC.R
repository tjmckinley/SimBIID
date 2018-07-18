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
#'
#' @return A plot of the ABC posterior distributions for different generations, or the distributions
#'         of the simulated summary measures for different generations.

plot.ABCSMC <- function(x, type = c("post", "output"), gen = NA) {
    
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
    
    if(type == "post") {
        ## generate colorRamp
        getPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))
        fillCols <- getPalette(length(x$pars))
        fillCols <- fillCols[as.numeric(gen)]
        
        ## plot posteriors
        p <- x$pars %>%
            map(~{
                as_tibble(.) %>%
                set_names(paste("par", 1:ncol(.)))
            }) %>%
            bind_rows(.id = "Generation") %>%
            filter(Generation %in% gen) %>%
            mutate(Generation = factor(Generation, levels = gen)) %>%
            gather(Parameter, value, -Generation) %>%
            ggplot(aes(x = value, fill = Generation)) +
                geom_density(alpha = 0.8) +
                xlab("Parameter value") + ylab("Density") +
                scale_fill_manual(values = fillCols) +
                facet_wrap(~ Parameter, scales = "free")
     } else {
        ## generate colorRamp
        getPalette <- colorRampPalette(brewer.pal(9, "YlGnBu"))
        fillCols <- getPalette(length(x$pars))
        fillCols <- fillCols[as.numeric(gen)]
        
        ## plot outputs
        p <- x$output %>%
            map(~{
                as_tibble(.) %>%
                set_names(paste("output", 1:ncol(.)))
            }) %>%
            bind_rows(.id = "Generation") %>%
            filter(Generation %in% gen) %>%
            mutate(Generation = factor(Generation, levels = gen)) %>%
            gather(Output, value, -Generation) %>%
            ggplot(aes(x = value, fill = Generation)) +
                geom_density(alpha = 0.8) +
                xlab("Output") + ylab("Density") +
                scale_fill_manual(values = fillCols) +
                facet_wrap(~ Output, scales = "free")
     }
     print(p)
}
    
