## simple ABC-SMC algorithm of Toni et al. (2009)

library(mvtnorm)
library(GillespieSSA)

#' @priors  Matrix of upper and lower bounds for uniform priors
ABC_SMC <- function(npart, tols, priors, func, data) {
    
    ## set initial particles
    accrate <- rep(NA, nrow(tols))
    pars <- list(); out <- list(); 
    pars[[1]] <- matrix(NA, npart, nrow(priors))
    out[[1]] <- matrix(NA, npart, length(data))
    accrate[1] <- 0
    cat("Generation 1:\n")
    for(i in 1:npart) {
        valid <- 0
        while(valid == 0) {
            ## sample from prior
            for(j in 1:nrow(priors)) {
                pars[[1]][i, j] <- runif(1, priors[, 1], priors[, 2])
            }
            
            ## simulate from model
            out[[1]][i, ] <- func(pars[[1]][i, ])
            
            ## check matching## check matching
            if(all(!is.na(out[[1]][i, ]))) {
                if(all(abs(out[[1]][i, ] - data) < tols[1, ])) {
                    valid <- 1
                }
            }
            
            ## update counter
            accrate[1] <- accrate[1] + 1
        }
        if(i %% 10 == 0) {
            cat(paste0("i = ", i, ", accrate = ", signif(i / accrate[1], 2), "\n"))
        }
    }
    
    ## set accrate
    accrate[1] <- npart / accrate[1] 
    
    ## set weights
    weights <- rep(1 / npart, npart)
    
    ## set proposal covariance
    propCov <- cov(pars[[1]]) * 2
    
    ## run sequential algorithm
    for(t in 2:nrow(tols)) {
        cat(paste0("Generation ", t, ":\n"))
        pars[[t]] <- matrix(NA, npart, nrow(priors))
        out[[t]] <- matrix(NA, npart, length(data))
        accrate[t] <- 0
        weightsNew <- rep(NA, npart)
        for(i in 1:npart) {
            valid <- 0
            while(valid == 0) {
                ## sample from previous generation
                k <- sample(1:npart, 1, prob = weights)
                pars[[t]][i, ] <- rmvnorm(1, 
                    mean = pars[[t - 1]][k, ], 
                    sigma = propCov
                )
                
                if(all(apply(cbind(pars[[t]][i, ], priors), 1, function(x) {
                    x[1] > x[2] & x[1] < x[3]
                }))) {
                    ## simulate from model
                    out[[t]][i, ] <- func(pars[[t]][i, ])
                    
                    ## check matching
                    if(all(!is.na(out[[t]][i, ]))) {
                        if(all(abs(out[[t]][i, ] - data) < tols[t, ])) {
                            valid <- 1
                        }
                    }
                }
                
                ## update counter
                accrate[t] <- accrate[t] + 1
            }
            
            ## calculate unnormalised weight
            weightsNew[i] <- prod(apply(cbind(pars[[t]][i, ], priors), 1, function(x) {
                dunif(x[1], x[2], x[3])
            }))
            weightsNew[i] <- weightsNew[i] / sum(weights * apply(pars[[t - 1]], 1, function(x) {
                dmvnorm(pars[[t]][i, ], mean = x, sigma = propCov)
            }))
            if(i %% 10 == 0) {
                cat(paste0("i = ", i, ", accrate = ", signif(i / accrate[t], 2), "\n"))
            }
        }
        
        ## set accrate
        accrate[t] <- npart / accrate[t] 
        
        ## set weights
        weights <- weightsNew / sum(weightsNew)
        
        ## set proposal covariance
        propCov <- cov(pars[[t]]) * 2
    }
    
    ## output results
    list(pars = pars, out = out, accrate = accrate)
}

## read in the data and convert to vector
smallpox <- as.numeric(unlist(read.table("smallpox.txt", header = F)))
summary(smallpox)

## set up function to perform simulation
simSIR <- function(pars) {
    ## set up parameters in order to use Gillespie SSA
    parameters <- c(beta = pars[1], gamma = pars[2])
    
    ## set initial states (1 initial infection in population of 120)
    initial_state <- c(S = 119, I = 1, R = 0)
    
    ## set rate vector and transition matrix for ssa function
    rates <- c("beta*S*I", "gamma*I")
    events <- matrix(c(-1, 1, 0, 0, -1, 1), nrow = 3)
    
    ## run ssa function to produce simulation
    ## set large final time to limit runtime
    sims <- ssa(initial_state, rates, events, parameters, tf = 2000)
    
    ## extract simulation from output
    sims <- sims$data
    
    ## extract final epidemic size and date of final removal
    finalsize <- sims[nrow(sims), 4]
    finaltime <- sims[nrow(sims), 1]
    
    ## reset time to be relative to first removal time
    firstremtime <- try(sims[sims[, 4] > 0, ][1, 1])
    if(class(firstremtime) == "try-error") {
        return(rep(NA, 2))
    }
    finaltime <- finaltime - firstremtime
    
    # return summary statistics
    return(c(finalsize, finaltime))
}

output <- ABC_SMC(100, matrix(rep(seq(100, 20, by = -10), each = 2), ncol = 2, byrow = T), 
    matrix(c(0, 0.002, 0, 0.2), 2, 2, byrow = T), 
    simSIR, 
    c(length(smallpox), smallpox[length(smallpox)])
)

## plot function
library(tidyverse)

plotABC <- function(ABC) {
    p <- ABC$pars %>%
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

plotABC(output)
            
