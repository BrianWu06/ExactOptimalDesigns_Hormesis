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
    
    # This is an example code for finding optimal exact designs for different models under D, h, and tau criterion.
    # For the detials of these criterions, please refer to the shiny app main page or our paper 
    # "Exact Optimal Designs for Small Studies in Toxicology with Applications to Hormesis via Metaheuristics"
    # The main functions are in the source code "HROD_PSO.R".
    # For the globpso package, please refer to https://github.com/PingYangChen/globpso for more details.
    library(globpso)
    source("HEOD_PSO.R")
    
    
    ## Hunt-Bowman models
    
    # Set Hunt-Bowman model parameters and PSO settings.
    hb_par <- hb_parms(c1 = 170, tau = 0.04, b0 = 1.86, b1 = 40) 
    psoinfo_hb <- psoinfo_setting(nSwarms = 128, Iters = 2000)
    
    # Replicate 4-Point D-Optimal Exact Design for Hunt-Bowman model 5 times.
    hb_doptiaml_res <- hb_doptimal_pso_rep(nRep = 5, nPoints = 4, parms = hb_par, psoinfo = psoinfo_hb)
    hb_doptiaml_res$result
    hb_doptiaml_res$approximate_design
    
    # Replicate 4-Point h-Optimal Exact Design for Hunt-Bowman model 5 times.
    hb_hoptiaml_res <- hb_hoptimal_pso_rep(nRep = 5, nPoints = 4, parms = hb_par, psoinfo = psoinfo_hb)
    hb_hoptiaml_res$result
    hb_hoptiaml_res$approximate_design
    
    # Replicate 2-Point tau-Optimal Exact Design for Hunt-Bowman model 5 times.
    hb_tauoptiaml_res <- hb_tauoptimal_pso_rep(nRep = 5, nPoints = 2, parms = hb_par, psoinfo = psoinfo_hb)
    hb_tauoptiaml_res$result
    hb_tauoptiaml_res$approximate_design
    
    
    ## exp-log models
    
    # Set exp-log model parameters and PSO settings.
    el_par <- exp_log_params(0.15, 89, 3.2, 41) 
    psoinfo_el <- psoinfo_setting(nSwarms = 128, Iters = 2000)
    
    # Replicate 4-Point D-Optimal Exact Design for exp-log model 5 times.
    el_doptimal_res <- exp_log_doptimal_pso_rep(nRep = 5, nPoints = 4, parms = el_par, psoinfo = psoinfo_el)
    el_doptimal_res$result
    el_doptimal_res$approximate_design
    
    # Replicate 4-Point h-Optimal Exact Design for exp-log model 5 times.
    el_hoptimal_res <- exp_log_hoptimal_pso_rep(nRep = 5, nPoints = 4, parms = el_par, psoinfo = psoinfo_el)
    el_hoptimal_res$result
    el_hoptimal_res$approximate_design
    
    # Replicate 2-Point tau-Optimal Exact Design for exp-log model 5 times.
    el_tauoptimal_res <- exp_log_tauoptimal_pso_rep(nRep = 5, nPoints = 2, parms = el_par, psoinfo = psoinfo_el)
    el_tauoptimal_res$result
    el_tauoptimal_res$approximate_design
    
    
    ## Logistic models
    
    # Set simple logistic model parameters and PSO settings.
    log_par <- logistic_params(alpha = 2, beta = 1)
    psoinfo_log <- psoinfo_setting(nSwarms = 128, Iters = 1000)
    
    # Replicate 2-Point D-Optimal Exact Design for simple logistic model 5 times.
    log_res <- logistic_doptimal_pso_rep(nRep = 5, nPoints = 2, parms = log_par, psoinfo = psoinfo_log)
    log_res$result
    log_res$approximate_design
    
    # Quadratic logistic model
    # Set quadratic logistic model parameters and PSO settings.
    qlog_par <- qlogistic_params(alpha = 3, beta1 = 0, beta2 = -1)
    psoinfo_qlog <- psoinfo_setting(nSwarms = 128, Iters = 1000)
    
    # Replicate 4-Point D-Optimal Exact Design for quadratic logistic model 5 times.
    qlog_res <- qlogistic_pso_rep(nRep = 1, nPoints = 4, parms = qlog_par, psoinfo = psoinfo_qlog)
    qlog_res$result
    qlog_res$approximate_design
    
    # Cubic logistic model
    # Set cubic logistic model parameters and PSO settings.
    clog_par <- clogistic_params(alpha = 1, beta1 = 3, beta2 = 2, beta3 = -1)
    psoinfo_clog <- psoinfo_setting(nSwarms = 128, Iters = 1000)
    
# Replicate 5-Point D-Optimal Exact Design for cubic logistic model 5 times.
clog_res <- clogistic_pso_rep(nRep = 5, nPoints = 5, parms = clog_par, psoinfo = psoinfo_clog)
clog_res$result
clog_res$approximate_design

