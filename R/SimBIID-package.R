#' Simulation-based inference for infectious disease models
#'
#' Package implements various simulation-based inference routines for infectious
#' disease models.
#'
#' Package provides some code to run the ABC-SMC algorithm of Toni et al. (2009)
#' and the Alive Particle Filter based particle MCMC algorithm (Jasra et al., 2013;
#' Drovandi et al., 2016). Also provides functions to plot the outputs.
#'
#' @docType package
#' @name SimBIID-package
#' @author TJ McKinley <t.mckinley@@exeter.ac.uk>
#' @import dplyr
#' @import purrr
#' @import tibble
#' @import ggplot2
#' @import tidyr
#' @import mvtnorm
#' @import RColorBrewer
#' @import coda 
#' @import Rcpp 
#' @import RcppArmadillo 
#' @import RcppXPtrUtils
#' @import SimInf
#' @useDynLib SimBIID
NULL

