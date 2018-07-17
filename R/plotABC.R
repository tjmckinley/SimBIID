#' @title Plots \code{ABCSMC} objects
#'
#' @description Plot method for \code{ABCSMC} objects.
#'
#' @export
#'
#' @param x     An \code{ABCSMC} object.
#'
#' @return A plot of the ABC posterior distributions for different generations.
#'

plot.ABCSMC <- function(x, ...) {
    p <- x$pars %>%
        map(~{
            as_tibble(.) %>%
            set_names(paste("par", 1:ncol(.)))
        }) %>%
        bind_rows(.id = "Generation") %>%
        gather(Parameter, value, -Generation) %>%
        ggplot(aes(x = value, fill = Generation)) +
            geom_density(alpha = 0.8) +
            xlab("Parameter value") + ylab("Density") +
            scale_fill_brewer(palette = "YlOrRd") +
            facet_wrap(~ Parameter, scales = "free")
     print(p)
}
    
