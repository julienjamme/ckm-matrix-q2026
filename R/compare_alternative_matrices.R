library(ptable)
library(ckm)
library(purrr)
library(dplyr)

# Parameters tested ----------------

source("R/parameters.R")

# Build transition matrices -----------------------  

ckm_matrices <- pmap(ckm_parameters, ptable::create_cnt_ptable, optim = 5, .progress=TRUE)

# Assess risk with alternative biased version -------------------------

ckm_matrices[[1]]@empResults
