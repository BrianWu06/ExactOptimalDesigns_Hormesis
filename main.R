# This is an example code for finding optimal exact designs for different models under D, h, and tau criterion.
# For the detials of these criterions, please refer to the shiny app main page or our paper 
# "Exact Optimal Designs for Small Studies in Toxicology with Applications to Hormesis via Metaheuristics"
# The main functions are in the source code "HROD_PSO.R".
# For the globpso package, please refer to https://github.com/PingYangChen/globpso for more details.
library(globpso)
library(dplyr)
library(ggplot2)
library(MASS)
source("EODH.R")


# Set the number of particles and iterations for PSO
psoinfo <- psoinfo_setting(nSwarms = 64, Iters = 1000)

# Run the hormesis_pso for searching the optimal exact designs
# model options: "HuntBowman", "ExpLog", "Logistic", "qlogistic", "clogistic"
# criterion options: "D" for all 5 models, "tau" and "h" for HuntBowman and ExpLog
# nPoints: The run size / number of observations of experiment
# parms: The parameter vector for the model
# psoinfo: The number of particles and iterations set by the psoinfo_setting function
# upper: Upper bound of the design space
# lower: Lower bound of the design space
# nRep: Number of reruns for PSO to search the optimal exact design

# An example of searching the D-optimal exact design for the Hunt-Bowman model for N = 10.
hb_par <- hb_parms(c1 = 170, tau = 0.04, b0 = 1.46, b1 = 40)
hb_res <- hormesis_pso(model = "HuntBowman", criterion = "D", nPoints = 10, parms = hb_par, 
                         psoinfo = psoinfo, upper = 0.15, lower = 0, nRep = 5)
hb_res

# Another example of searching the h-optimal exact design for the Exp-Log model for N = 10.
el_par <- el_parms(c0 = 0.15, c1 = 89, b0 = 3.2, b1 = 41)
el_res <- hormesis_pso(model = "ExpLog", criterion = "h", nPoints = 10, parms = el_par, 
                             psoinfo = psoinfo, upper = 0.15, lower = 0, nRep = 5)
el_res

# Examples of searching the D-optimal exact design for the logistic models for N = 10.
# Logistic model
log_par <- logistic_params(alpha = 2, beta = 1)
log_res <- hormesis_pso(model = "logistic", criterion = "D", nPoints = 10, parms = log_par, 
                        psoinfo = psoinfo, upper = 5, lower = -5, nRep = 5)
log_res

# Quadratic logistic model
qlog_par <- qlogistic_params(alpha = 3, beta1 = 0, beta2 = -1)
qlog_res <- hormesis_pso(model = "qlogistic", criterion = "D", nPoints = 10, parms = qlog_par, 
                         psoinfo = psoinfo, upper = 5, lower = -5, nRep = 5)
qlog_res

# Cubic logistic model
clog_par <- clogistic_params(alpha = 1, beta1 = 3, beta2 = 2, beta3 = -1)
clog_res <- hormesis_pso(model = "clogistic", criterion = "D", nPoints = 10, parms = clog_par, 
                         psoinfo = psoinfo, upper = 5, lower = -5, nRep = 5)
clog_res