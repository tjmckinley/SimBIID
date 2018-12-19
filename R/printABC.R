#' @title Prints \code{ABCSMC} objects
#'
#' @description Print method for \code{ABCSMC} objects.
#'
#' @param x    An \code{ABCSMC} object.
#' @param ...           Not used here.
#'
#' @return Summary outputs printed to the screen.
#' 
#' @export

print.ABCSMC <- function(x, ...) {
    ## check object is a ABCSMC object
    if(class(x) != "ABCSMC"){
        stop("'x' is not a ABCSMC object")
    }
    
    cat("An object of class: 'ABCSMC'\n")
    cat(paste0("Consists of ", nrow(x$tols), " generations with ", nrow(x$priors), " parameters.\n"))
    
    ## print data information
    cat("\nData:\n")
    print(x$data, row.names = F)
    
    ## print tolerance information
    temp <- x$tols %>%
        as.data.frame() %>%
        mutate(Generation = 1:n()) %>%
        select(Generation, everything()) %>%
        mutate(ESS = do.call("c", x$ESS))
    cat("\nTolerances:\n")
    print(temp, row.names = F)
    
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

