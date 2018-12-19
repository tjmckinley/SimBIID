#' @title Prints \code{SimBIID_runs} objects
#'
#' @description Print method for \code{SimBIID_runs} objects.
#'
#' @param x    A \code{SimBIID_runs} object.
#' @param ...           Not used here.
#'
#' @return Summary outputs printed to the screen.
#' 
#' @export

print.SimBIID_runs <- function(x, ...) {
    cat(paste("'SimBIID_runs' object with n =", nrow(x$sums), "replicates.\n"))
    if(nrow(x$sums) == 1){
        cat("\nOutput at final time point:\n")
        x$sums %>%
            select(-rep) %>%
            as.tibble() %>%
            print()
        if(is.data.frame(x$runs)) {
            cat("\nTime-series counts:\n")
            x$runs %>%
                select(-rep) %>%
                as.tibble() %>%
                print()
        }
    } else {
        cat("\nSummaries of outputs at final time point:\n")
        x$sums %>%
            select(-rep) %>%
            summary() %>%
            print()
        # if(is.data.frame(x$runs)) {
        #     cat("\nSummaries of time-series counts:\n")
        #     x$runs %>%
        #         select(-rep) %>%
        #         gather(output, value, -t) %>%
        #         group_by(t, output) %>%
        #         summarise(list(enframe(c(mean(value), quantile(value, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))))) %>%
        #         unnest() %>%
        #         spread(name, value) %>%
        #         print()
        # }
    }
}
