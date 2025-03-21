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
psoinfo_approx <- psoinfo_setting(nSwarms = 64, Iters = 1000)
psoinfo_exact <- psoinfo_setting(nSwarms = 32, Iters = 500)

# Run the hormesis_pso for searching the optimal exact designs
# model options: "HuntBowman" and "ExpLog"
# criterion options: "D", "tau" and "h"
# nPoints: The run size / number of observations of experiment
# parms: The parameter vector for the model
# psoinfo: The number of particles and iterations set by the psoinfo_setting function
# upper: Upper bound of the design space
# lower: Lower bound of the design space
# nRep: Number of reruns for PSO to search the optimal exact design

# An example of searching the D-optimal exact design for the Hunt-Bowman model for N = 10.
hb_par <- hb_parms(c1 = 170, tau = 0.04, b0 = 1.46, b1 = 40)
hb_res <- hormesis_pso(model = "HuntBowman", criterion = "h", nPoints = 10, parms = hb_par, 
                       psoinfo_exact = psoinfo_exact, psoinfo_approx = psoinfo_approx, 
                       upper = 0.15, lower = 0, nRep = 1)
hb_res

# Another example of searching the h-optimal exact design for the Exp-Log model for N = 10.
el_par <- el_parms(c0 = 0.15, c1 = 89, b0 = 3.2, b1 = 41)
el_res <- hormesis_pso(model = "ExpLog", criterion = "h", nPoints = 10, parms = el_par, 
                       psoinfo_exact = psoinfo_exact, psoinfo_approx = psoinfo_approx,  
                       upper = 0.15, lower = 0, nRep = 1)
el_res