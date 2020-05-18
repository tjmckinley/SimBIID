test_that("ABCSMC works", {
    
    ## set up SIR simulationmodel
    transitions <- c(
        "S -> beta * S * I -> I",
        "I -> gamma * I -> R"
    )
    compartments <- c("S", "I", "R")
    pars <- c("beta", "gamma")
    model <- mparseRcpp(
        transitions = transitions,
        compartments = compartments,
        pars = pars
    )
    model <- compileRcpp(model)

    ## generate function to run simulators
    ## and return summary statistics
    simSIR <- function(pars, data, tols, u, model) {

        ## run model
        sims <- model(pars, 0, data[2] + tols[2], u)

        ## this returns a vector of the form:
        ## completed (1/0), t, S, I, R (here)
        if(sims[1] == 0) {
            ## if simulation rejected
            return(NA)
        } else {
            ## extract finaltime and finalsize
            finaltime <- sims[2]
            finalsize <- sims[5]
        }

        ## return vector if match, else return NA
        if(all(abs(c(finalsize, finaltime) - data) <= tols)){
            return(c(finalsize, finaltime))
        } else {
            return(NA)
        }
    }

    ## set priors
    priors <- data.frame(
        parnames = c("beta", "gamma"),
        dist = rep("gamma", 2),
        stringsAsFactors = FALSE
    )
    priors$p1 <- c(10, 10)
    priors$p2 <- c(10^4, 10^2)

    ## define the targeted summary statistics
    data <- c(
        finalsize = 30,
        finaltime = 76
    )

    ## set initial states (1 initial infection
    ## in population of 120)
    iniStates <- c(S = 119, I = 1, R = 0)

    ## set initial tolerances
    tols <- c(
        finalsize = 50,
        finaltime = 50
    )
    
    set.seed(50)

    ## run 2 generations of ABC-SMC
    ## setting tolerance to be 50th
    ## percentile of the accepted
    ## tolerances at each generation
    post <- ABCSMC(
        x = data,
        priors = priors,
        func = simSIR,
        u = iniStates,
        tols = tols,
        ptol = 0.2,
        ngen = 2,
        npart = 50,
        model = model
    )
    ## file to save results
    tmp <- "ABCSMC"
    ## the first run always succeeds, but warns
    expect_known_output(post, tmp, print = TRUE)

    ## run one further generation
    post <- ABCSMC(post, ptols = 0.5, ngen = 1)
    ## file to save results
    tmp <- "ABCSMC1"
    ## the first run always succeeds, but warns
    expect_known_output(post, tmp, print = TRUE)
    
    ## file to save results
    tmp <- "ABCSMCsum"
    ## the first run always succeeds, but warns
    expect_known_output(summary(post), tmp, print = TRUE)
    
    ## try to run 2 generations of ABC-SMC using fixed tolerances
    tols <- matrix(c(50, 50, 50, 50), 2, 2, byrow = TRUE)
    colnames(tols) <- c("finalsize", "finaltime")
    expect_error(ABCSMC(
        x = data,
        priors = priors,
        func = simSIR,
        u = iniStates,
        tols = tols,
        ngen = 2,
        npart = 50,
        model = model
    ))
    tols <- matrix(c(50, 50, 51, 20), 2, 2, byrow = TRUE)
    colnames(tols) <- c("finalsize", "finaltime")
    expect_error(ABCSMC(
        x = data,
        priors = priors,
        func = simSIR,
        u = iniStates,
        tols = tols,
        ngen = 2,
        npart = 50,
        model = model
    ))
    
    ## run 2 generations of ABC-SMC using fixed tolerances
    tols <- matrix(c(50, 50, 50, 20), 2, 2, byrow = TRUE)
    colnames(tols) <- c("finalsize", "finaltime")
    post <- ABCSMC(
        x = data,
        priors = priors,
        func = simSIR,
        u = iniStates,
        tols = tols,
        ngen = 2,
        npart = 50,
        model = model
    )
    ## file to save results
    tmp <- "ABCSMCtol"
    ## the first run always succeeds, but warns
    expect_known_output(post, tmp, print = TRUE)
    
    ## run one further generation
    newtols <- c(finalsize = 20, finaltime = 15)
    post <- ABCSMC(post, tols = newtols)
    ## file to save results
    tmp <- "ABCSMCtol1"
    ## the first run always succeeds, but warns
    expect_known_output(post, tmp, print = TRUE)
    
    ## run one further generation
    expect_error(ABCSMC(post, tols = newtols))
    newtols <- c(finalsize = 21, finaltime = 15)
    expect_error(ABCSMC(post, tols = newtols))
    
    ## run further generation using ptols
    expect_error(ABCSMC(post, tols = newtols, ptols = 0.8, ngen = 1))
    post <- ABCSMC(post, ptols = 0.8, ngen = 1)
    ## file to save results
    tmp <- "ABCSMCtol2"
    ## the first run always succeeds, but warns
    expect_known_output(post, tmp, print = TRUE)
    
    ## run 2 generations of ABC-SMC using fixed tolerances
    tols <- matrix(c(50, 50), 1, 2, byrow = TRUE)
    colnames(tols) <- c("finalsize", "finaltime")
    post <- ABCSMC(
        x = data,
        priors = priors,
        func = simSIR,
        u = iniStates,
        tols = tols,
        npart = 50,
        model = model
    )
    
    ## run two further generations
    newtols <- matrix(c(50, 50, 40, 40), 2, 2, byrow = TRUE)
    colnames(newtols) <- c("finalsize", "finaltime")
    expect_error(ABCSMC(post, tols = newtols))
    
    ## run two further generations
    newtols <- matrix(c(45, 45, 45, 45), 2, 2, byrow = TRUE)
    colnames(newtols) <- c("finalsize", "finaltime")
    expect_error(ABCSMC(post, tols = newtols))
    
    ## run two further generations
    newtols <- matrix(c(45, 45, 40, 40), 2, 2, byrow = TRUE)
    colnames(newtols) <- c("finalsize", "finaltime")
    post <- ABCSMC(post, tols = newtols)
    
})
