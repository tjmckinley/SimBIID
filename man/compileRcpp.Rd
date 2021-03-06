% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compileRcpp.R
\name{compileRcpp}
\alias{compileRcpp}
\title{Compiles \code{SimBIID_model} object}
\usage{
compileRcpp(model)
}
\arguments{
\item{model}{An object of class \code{SimBIID_model}.}
}
\value{
An object of class \code{XPtr} that points to the compiled function, or
        an R \code{function} object for calling directly from R.
}
\description{
Compiles an object of class \code{SimBIID_model} into an
             \code{XPtr} object for use in Rcpp functions, or an
             object of class \code{function} for calling directly from R.
}
\examples{
\donttest{
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

## compile model to be run directly
model <- compileRcpp(model)
model

## set initial states (1 initial infection 
## in population of 120)
iniStates <- c(S = 119, I = 1, R = 0)

## set parameters
pars <- c(beta = 0.001, gamma = 0.1)

## run compiled model
model(pars, 0, 100, iniStates)
}

}
\seealso{
\code{\link{mparseRcpp}}
}
