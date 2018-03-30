#-------------------------------------------------------------------------
#  This file is part of the TVAR forecasting application https://github.com/maximegodin/tvar_forecasting
#
#  This application is governed by the MIT license. 
#  You can  use, modify and/ or redistribute this code under the terms
#  of the MIT license:  https://opensource.org/licenses/MIT
#
#  Maxime Godin, Amina Bouchafaa, Ecole Polytechnique
#  March 29th, 2018
#-------------------------------------------------------------------------

library(Rcpp)
library(reshape2)
library(dplyr)
library(ggplot2)
library(rmutil)
library(extraDistr)
library(xtable)
library(ggforce)
library(forecast)
library(gridExtra)

options(xtable.floating = FALSE)
options(xtable.timestamp = "")

sourceCpp("TVAR_tools.cpp")
sourceCpp("local_yule_walker.cpp")
sourceCpp("NLMS.cpp")
sourceCpp("exp_agg_NLMS.cpp")
source("graphical_outputs.R")
source("TVAR_tools.R")