# Load the functions in "HEOD_PSO.R" for the functions of the following examples.
# For the globpso package, please refer to https://github.com/PingYangChen/globpso for more details.
library(globpso)
source("HEOD_PSO.R")


# Hunt-Bowman d-optimal

# Set Hunt-Bowman model parameters and PSO settings.
hb_par <- hb_parms(c1 = 170, tau = 0.04, b0 = 1.46, b1 = 40) 
psoinfo_hb <- psoinfo_setting(nSwarms = 64, Iters = 1000)

# Replicate Hunt-Bowman model 4-Point D-Optimal Exact Design 10 times.
hb_res <- hb_pso_rep(nRep = 10, nPoints = 4, parms = hb_par, psoinfo = psoinfo_hb)

# Results of Hunt-Bowman model 4-Point D-Optimal Exact Design
hb_res$result$best_val
hb_res$result$best_points


# exp-log Model

# Set exp-log model parameters and PSO settings.
el_par <- exp_log_params(0.15, 89, 3.2, 41) 
psoinfo_el <- psoinfo_setting(nSwarms = 128, Iters = 2000)


# Replicate Hunt-Bowman model 4-Point h-Optimal Exact Design 10 times.
el_res <- exp_log_h_pso_rep(nRep = 10, nPoints = 4, parms = el_par, psoinfo = psoinfo_el)

# Results of Hunt-Bowman model 4-Point h-Optimal Exact Design
el_h_res$result$best_val
el_h_res$result$best_points


# Replicate Hunt-Bowman model 2-Point tau-Optimal Exact Design 10 times.
psoinfo_el <- psoinfo_setting(nSwarms = 128, Iters = 1000)
el_tau_res <- exp_log_tau_pso_rep(nRep = 10, nPoints = 2, parms = el_par, psoinfo = psoinfo_el)

# Results of Hunt-Bowman model 2-Point tau-Optimal Exact Design
el_tau_res$result$best_val
el_tau_res$result$best_points


# Simple logistic model
log_par <- logistic_params(alpha = 2, beta = 1)
psoinfo_log <- psoinfo_setting(nSwarms = 128, Iters = 1000)

log_res <- logistic_d_pso_rep(nRep = 10, nPoints = 2, parms = log_par, psoinfo = psoinfo_log)

log_res$result$best_val
log_res$result$best_points

#Quadratic logistic model
qlog_par <- qlogistic_params(alpha = -3, beta1 = 0, beta2 = -1)
psoinfo_qlog <- psoinfo_setting(nSwarms = 128, Iters = 1000)

qlog_res <- qlogistic_pso_rep(nRep = 10, nPoints = 4, parms = qlog_par, psoinfo = psoinfo_qlog)

qlog_res$result$best_val
qlog_res$result$best_points
