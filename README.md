# Optimal experimental design strategies for detecting hormesis
This repository contains the sources code of the app [Optimal experimental design strategies for detecting hormesis](https://brianwu.shinyapps.io/hormeis_ed_pso/).

The source code is written in the R programming language. Users are required to download the package [globpso](https://github.com/willgertsch/SingleObjApp/blob/main/app.R).
The package can be downloaded by the following code.
    install.packages("devtools")
    devtools::install_github("PingYangChen/globpso")

In our app, we only allow users to search for up to 15 number of design points for all 5 optimal exact design problems mentioned in our paper to avoid exhausting the avaiable computing resource on the server. 
Hence, we suggest useres to adjust and run the codes on their own device for larger number of design points.

For more details of the background knowledge of our works, please refer to the user manual page of our [app](https://brianwu.shinyapps.io/hormeis_ed_pso/). 

The following R code (main.R) is an example of how to use the source code. 
    ##Filename: main.R
    # Load the required packages and the main functions
    library(globpso)
    source("HEOD_PSO.R")
    
    # Hunt-Bowman d-optimal
    
    # Set Hunt-Bowman model parameters and PSO settings.
    hb_par <- hb_parms(c1 = 170, tau = 0.04, b0 = 1.46, b1 = 40) 
    psoinfo_hb <- psoinfo_setting(nSwarms = 64, Iters = 1000)
    
    # Replicate 4-Point D-Optimal Exact Design for Hunt-Bowman model 10 times.
    hb_res <- hb_pso_rep(nRep = 10, nPoints = 4, parms = hb_par, psoinfo = psoinfo_hb)
    
    # Results of 4-Point D-Optimal Exact Design for Hunt-Bowman model
    hb_res$result$best_val
    hb_res$result$best_points
    
    
    # exp-log Model
    
    # Set exp-log model parameters and PSO settings.
    el_par <- exp_log_params(0.15, 89, 3.2, 41) 
    psoinfo_el <- psoinfo_setting(nSwarms = 128, Iters = 2000)
    
    
    # Replicate 4-Point h-Optimal Exact Design for exp-log model 10 times.
    el_res <- exp_log_h_pso_rep(nRep = 10, nPoints = 4, parms = el_par, psoinfo = psoinfo_el)
    
    # Results of 4-Point h-Optimal Exact Design for exp-log model
    el_h_res$result$best_val
    el_h_res$result$best_points
    
    
    # Replicate 2-Point tau-Optimal Exact Design for exp-log model 10 times.
    psoinfo_el <- psoinfo_setting(nSwarms = 128, Iters = 1000)
    el_tau_res <- exp_log_tau_pso_rep(nRep = 10, nPoints = 2, parms = el_par, psoinfo = psoinfo_el)
    
    # Results of 2-Point tau-Optimal Exact Design for exp-log model
    el_tau_res$result$best_val
    el_tau_res$result$best_points
    
    
    # Simple logistic model

    # Set simple logistic model parameters and PSO settings.
    log_par <- logistic_params(alpha = 2, beta = 1)
    psoinfo_log <- psoinfo_setting(nSwarms = 128, Iters = 1000)

    # Replicate 2-point D-optimal exact design for simple logistic model 10 times.
    log_res <- logistic_d_pso_rep(nRep = 10, nPoints = 2, parms = log_par, psoinfo = psoinfo_log)

    # Results of 2-point D-optimal exact design for simple logistic model
    log_res$result$best_val
    log_res$result$best_points
    
    #Quadratic logistic model

    # Set quadratic logistic model parameters and PSO settings.
    qlog_par <- qlogistic_params(alpha = -3, beta1 = 0, beta2 = -1)
    psoinfo_qlog <- psoinfo_setting(nSwarms = 128, Iters = 1000)

    # Replicate 4-point D-optimal exact design for quadratic logistic model 10 times.
    qlog_res <- qlogistic_pso_rep(nRep = 10, nPoints = 4, parms = qlog_par, psoinfo = psoinfo_qlog)

    # Results of 4-point D-optimal exact design for quadratic logistic model
    qlog_res$result$best_val
    qlog_res$result$best_points
