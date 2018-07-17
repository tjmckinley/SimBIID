# Checks input arguments for consistency
checkInput <- function(input, type, length = NA, 
        nrow = NA, ncol = NA, int = F, naAllow = F) {
    # type is set of characters denoting is.*() type functions to test input
    # length is numeric testing length of input (similarly for nrow and ncol)
    # int is logical testing for integer
    # naAllow determines whether NA values are allowed or not
     
    type <- paste0("is.", type)
    output <- sapply(type, function(f, x) {
        do.call(f, list(x))
    }, x = input)
    if(!is.na(length)) {
        output <- c(output, length(input) == length)
    }
    if(!is.na(nrow)) {
        output <- c(output, nrow(input) == nrow)
    }
    if(!is.na(ncol)) {
        output <- c(output, ncol(input) == ncol)
    }
    if(int) {
        output <- c(output, all(input %% 1 == 0))
    }
    if(!naAllow){
        if(!is.function(input)) {
            output <- c(output, all(!is.na(input)))
        }
    }
    all(output)
}
