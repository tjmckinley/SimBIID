#' Sequential Monte Carlo based ABC
#'
#' Package implements various inference routines
#'
#' Package provides some simple code to run the ABCSMC algorithm of Toni et al. (2009)
#' and the Alive Particle Filter based particle MCMC algorithm (Jasra et al., 2013;
#' Drovandi et al., 2016). Also provides functions to plot the outputs.
#'
#' @docType package
#' @name ABCSMC-package
#' @author TJ McKinley <t.mckinley@@exeter.ac.uk>
#' @import dplyr
#' @import purrr
#' @import ggplot2
#' @import tidyr
#' @import mvtnorm
#' @import RColorBrewer
#' @import coda 
#' @import Rcpp 
#' @import RcppArmadillo 
#' @import RcppXPtrUtils
#' @useDynLib ABCSMC
NULL

