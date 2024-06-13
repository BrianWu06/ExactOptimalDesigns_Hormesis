# Optimal experimental design strategies for detecting hormesis
This repository contains the sources code of the app [Optimal experimental design strategies for detecting hormesis](https://brianwu.shinyapps.io/hormeis_ed_pso/).

The source code is written in the R programming language. Users are required to download the package [globpso](https://github.com/willgertsch/SingleObjApp/blob/main/app.R).
The package can be downloaded by the following code. 

    install.packages("devtools")
    devtools::install_github("PingYangChen/globpso")

In our app, we only allow users to search for up to 15 number of design points for all 5 optimal exact design problems mentioned in our paper to avoid exhausting the avaiable computing resource on the server. 
Hence, we suggest useres to adjust and run the codes on their own device for larger number of design points.

For more details of the background knowledge of our works, please refer to the user manual page of our [app](https://brianwu.shinyapps.io/hormeis_ed_pso/). 

The following R code (main.R) is an example of how to use the source code for the desired optimal exact designs. 

    # This is an example code for finding optimal exact designs for different models under D, h, and tau criterion.
    # For the detials of these criterions, please refer to the shiny app main page or our paper 
    # "Exact Optimal Designs for Small Studies in Toxicology with Applications to Hormesis via Metaheuristics"
    # The main functions are in the source code "HROD_PSO.R".
    # For the globpso package, please refer to https://github.com/PingYangChen/globpso for more details.
    library(globpso)
    library(dplyr)
    library(ggplot2)
    library(MASS)
    library(AlgDesign)
    source("EODH.R")
    
    
    # Run the hormesis_pso for searching the optimal exact designs
    # model options: "HuntBowman", "ExpLog", "Logistic", "qlogistic", "clogistic"
    # criterion options: "D" for all 5 models, "tau" and "h" for HuntBowman and ExpLog
    # nPoints: The run size / number of observations of experiment
    # parms: The parameter vector for the model
    # psoinfo: The number of particles and iterations set by the psoinfo_setting function
    # upper: Upper bound of the design space
    # lower: Lower bound of the design space
    # nRep: Number of reruns for PSO to search the optimal exact design
    
    # Set the number of particles and iterations for PSO
    psoinfo_exact <- psoinfo_setting(32, 500)
    psoinfo_approx <- psoinfo_setting(64, 1000)
    
    # An example of searching the D-optimal exact design by PSO for the Hunt-Bowman model for N = 10.
    hb_par <- hb_parms(c1 = 170, tau = 0.04, b0 = 1.46, b1 = 40)
    hb_res <- hormesis_pso(model = "HuntBowman", criterion = "D", nPoints = 10, parms = hb_par, 
                           psoinfo_exact = psoinfo_exact, psoinfo_approx = psoinfo_approx, 
                           upper = 0.15, lower = 0, nRep = 5)
    hb_res
    
    # Another example of searching the h-optimal exact design for the Exp-Log model for N = 10.
    el_par <- el_parms(c0 = 0.15, c1 = 89, b0 = 3.2, b1 = 41)
    el_res <- hormesis_pso(model = "ExpLog", criterion = "h", nPoints = 10, parms = el_par, 
                           psoinfo_exact = psoinfo_exact, psoinfo_approx = psoinfo_approx,
                           upper = 0.15, lower = 0, nRep = 5)
    el_res
    
    # Examples of searching the D-optimal exact design for the logistic models for N = 10.
    # Logistic model
    log_par <- logistic_parms(alpha = 2, beta = 1)
    log_res <- hormesis_pso(model = "logistic", criterion = "D", nPoints = 10, parms = log_par, 
                            psoinfo_exact = psoinfo_exact, psoinfo_approx = psoinfo_approx,
                            upper = 5, lower = -5, nRep = 5)
    log_res
    
    # Quadratic logistic model
    qlog_par <- qlogistic_parms(alpha = 3, beta1 = 0, beta2 = -1)
    qlog_res <- hormesis_pso(model = "qlogistic", criterion = "D", nPoints = 10, parms = qlog_par, 
                             psoinfo_exact = psoinfo_exact, psoinfo_approx = psoinfo_approx,
                             upper = 5, lower = -5, nRep = 5)
    qlog_res
    
    # Cubic logistic model
    clog_par <- clogistic_parms(alpha = 1, beta1 = 3, beta2 = 2, beta3 = -1)
    clog_res <- hormesis_pso(model = "clogistic", criterion = "D", nPoints = 10, parms = clog_par, 
                             psoinfo_exact = psoinfo_exact, psoinfo_approx = psoinfo_approx,
                             upper = 5, lower = -5, nRep = 5)
    clog_res
    
    
    # An example of searching the D-optimal exact design by DE for the Hunt-Bowman model for N = 10. 
    deinfo_exact <- getDEInfo(nPop = 32, maxIter = 500, deType = "rand-to-best-1")
    deinfo_approx <- getDEInfo(nPop = 64, maxIter = 1000, deType = "rand-to-best-1")
    
    hb_de_res <- hormesis_de(model = "HuntBowman", criterion = "D", nPoints = 10, parms = hb_par, 
                             deinfo_exact = deinfo_exact, deinfo_approx = deinfo_approx, 
                             upper = 0.15, lower = 0, nRep = 5)
        
