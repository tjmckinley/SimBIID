#' @title Simulation-based inference for infectious disease models
#'
#' @description Package implements various simulation-based inference routines for infectious
#' disease models.
#'
#' @details Package provides some code to run simulations of state-space models, and then
#' use these in the ABC-SMC algorithm of Toni et al. (2009)
#' and the bootstrap particle filter based particle MCMC algorithm (Andrieu et al., 2010). 
#' Also provides functions to plot and summarise the outputs.
#'
#' @docType package
#' @name SimBIID-package
#' @author Trevelyan J. McKinley <t.mckinley@@exeter.ac.uk>
#' @import dplyr
#' @import purrr
#' @import tibble
#' @import ggplot2
#' @import tidyr
#' @import Rcpp 
#' @importFrom Rcpp evalCpp
#' @importFrom graphics plot
#' @importFrom stats window
#' @useDynLib SimBIID, .registration = TRUE
NULL

