#' @title Prints \code{PMCMC} objects
#'
#' @description Print method for \code{PMCMC} objects.
#'
#' @param x    A \code{PMCMC} object.
#' @param ...       Not used here.
#'
#' @return Summary outputs printed to the screen.
#' @export

print.PMCMC <- function(x, ...) {
    ## check object is a PMCMC object
    if(class(x) != "PMCMC"){
        stop("'x' is not a PMCMC object")
    }
    
    cat("An object of class: 'PMCMC'\n")
    cat(paste0("Chain consists of ", nrow(x$pars), " iterations with ", nrow(x$priors), " parameters.\n"))
    
    ## print prior information
    temp <- x$priors %>%
        mutate(p1 = as.character(signif(p1, 2))) %>%
        mutate(p2 = as.character(signif(p2, 2))) %>%
        mutate(temp = ifelse(dist == "unif", paste0("U(lower = ", p1, ", upper = ", p2, ")"), NA)) %>%
        mutate(temp = ifelse(dist == "gamma", paste0("G(shape = ", p1, ", rate = ", p2, ")"), temp)) %>%
        mutate(temp = ifelse(dist == "norm", paste0("N(mean = ", p1, ", sd = ", p2, ")"), temp)) %>%
        mutate(temp = paste0(parnames, " ~ ", temp)) %>%
        select(temp)
    colnames(temp) <- ""
    cat("\nPriors:\n")
    print(temp, row.names = F, col.names = F, quote = F)
}
