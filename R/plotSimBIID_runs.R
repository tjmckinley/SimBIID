#' @title Plots \code{SimBIID_runs} objects
#'
#' @description Plot method for \code{SimBIID_runs} objects.
#'
#' @export
#'
#' @param x     An \code{SimBIID_runs} object.
#' @param which A character vector of states to plot. Can be \code{"all"} to plot all
#'              states (and final event times), or \code{"t"} to plot final event times.
#' @param type Character stating whether to plot full simulations over time (\code{"runs"}) or
#'             summaries (\code{"sums"}).
#' @param rep An integer vector of simulation runs to plot.
#' @param quant A vector of quantiles (> 0.5) to plot if \code{type == "runs"}.
#'
#' @return A plot of individual simulations and/or summaries of repeated simulations 
#'         extracted from \code{SimBIID_runs} object.

plot.SimBIID_runs <- function(x, which = c("all", "t"), type = c("runs", "sums"), 
                              rep = NA, quant = seq(0.55, 0.95, by = 0.05)) {
    ## check x
    if(class(x) != "SimBIID_runs"){
        stop("'x' is not a SimBIID_runs object")
    }
    ## check which
    whichset <- c("t", colnames(x$sums)[-c(1:3)])
    if(which[1] != "all") {
        checkInput(which, c("vector", "character"), inSet = whichset)
    } else {
        which <- whichset
    }
    ## check type
    checkInput(type, c("vector", "character"))
    type <- type[1]
    checkInput(type, inSet = c("runs", "sums"))
    ## check rep
    if(!is.na(rep[1])){
        checkInput(rep, c("vector", "numeric"), inSet = x$sums$rep, int = T)
    }
    ## check quant
    checkInput(quant, c("vector", "numeric"), inSet = seq(0.55, 0.95, by = 0.05))
    quant <- sort(quant)
    quant <- cbind(rev(1 - quant), rev(quant))
    quant1 <- sort(as.vector(quant))
    
    ## plot final times and final epidemic sizes
    if(!is.data.frame(x$runs) | type == "sums"){
        if(nrow(x$sums) < 20){
            stop("Can't plot summary of final sizes and times for n < 20 replicates")
        }
        p <- x$sums %>%
            select(!!which) %>%
            gather(output, value) %>%
            mutate(output = factor(output, levels = which)) %>%
            ggplot(aes(x = value)) +
                geom_histogram() +
                facet_wrap(~ output, scales = "free") +
                xlab("")
        if(!is.na(rep[1])){
            rep1 <- rep
            repSums <- x$sums %>%
                slice(rep1) %>%
                select(!!which) %>%
                gather(output, value) %>%
                mutate(output = factor(output, levels = which))
            p <- p + geom_point(aes(x = value), data = repSums, y = 0, colour = "red", shape = 16)
        }
    } else {
        ## produce plot
        which <- unique(c("rep", "t", which))
        if(nrow(x$sums) < 20) {
            ## produce plot
            p <- x$runs %>%
                gather(output, value, -rep, -t) %>%
                mutate(output = factor(output, levels = which[-match(c("rep", "t"), which)])) %>%
                ggplot(aes(x = t, y = value)) +
                xlab("Time") + ylab("Counts") +
                facet_wrap(~ output)
            ## add individual simulations if required
            rep1 <- unique(x$runs$rep)
            repSums <- x$runs %>%
                select(!!which) %>%
                gather(output, value, -rep, -t) %>%
                mutate(output = factor(output, levels = which[-match(c("rep", "t"), which)]))
            for(i in rep1) {
                temp <- filter(repSums, rep == i)
                p <- p + geom_line(aes(x = t, y = value), data = temp)
            }
        } else {
            p <- x$runs %>%
                select(!!which) %>%
                gather(output, value, -rep, -t) %>%
                group_by(output, t) %>%
                summarise(mean = mean(value), value = list(enframe(quantile(value, probs = quant1)))) %>%
                unnest() 
            p1 <- quant %>%
                as.data.frame() %>%
                rename(lci = V1, uci = V2) %>%
                mutate(pair = 1:n()) %>%
                gather(output, name, -pair) %>%
                mutate(name = 100 * name) %>%
                mutate(name = paste0(name, "%")) %>%
                inner_join(p, by = "name") %>%
                select(-name) %>%
                spread(output.x, value) %>%
                mutate(output.y = factor(output.y, levels = which[-match(c("rep", "t"), which)])) %>%
                mutate(pair = as.character(quant[pair, 2]))
           
           ## produce plot
           p <- p1 %>%
               ggplot(aes(x = t)) +
                   xlab("Time") + ylab("Counts")
           for(i in unique(p1$pair)){
                temp <- filter(p1, pair == i)
                p <- p + geom_ribbon(aes(ymin = lci, ymax = uci), 
                     data = temp, alpha = 0.2)
           }
           
           ## add individual simulations if required
           if(!is.na(rep[1])){
               rep1 <- rep
               repSums <- x$runs %>%
                   filter(rep %in% rep1) %>%
                   select(!!which) %>%
                   gather(output, value, -rep, -t) %>%
                   mutate(output.y = factor(output, levels = levels(p1$output.y)))
               for(i in rep1){
                   temp <- filter(repSums, rep == i)
                   p <- p + geom_line(aes(x = t, y = value), data = temp)
               }
           }
           p <- p + geom_line(aes(y = mean), colour = "red") +
               facet_wrap(~ output.y) +
               labs(title = paste0("Intervals = ", paste0(paste0(rev(quant[, 2]) * 100, "%"), collapse = ", ")))
        }
    }
    print(p)
}
    
