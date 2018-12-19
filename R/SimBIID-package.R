#' Simulation-based inference for infectious disease models
#'
#' Package implements various simulation-based inference routines for infectious
#' disease models.
#'
#' Package provides some code to run simulations of state-space models, and then
#' use these in the ABC-SMC algorithm of Toni et al. (2009)
#' and the bootstrap particle filter based particle MCMC algorithm (Andrieu et al., 2010). 
#' Also provides functions to plot and summarise the outputs.
#'
#' @docType package
#' @name SimBIID-package
#' @author TJ McKinley <t.mckinley@@exeter.ac.uk>
#' @import stats
#' @import dplyr
#' @import purrr
#' @import tibble
#' @import ggplot2
#' @import tidyr
#' @import mvtnorm
#' @import grDevices
#' @import RColorBrewer
#' @import coda 
#' @import Rcpp 
#' @import RcppArmadillo 
#' @import RcppXPtrUtils
#' @import SimInf
#' @useDynLib SimBIID
NULL

