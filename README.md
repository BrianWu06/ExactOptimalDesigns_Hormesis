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
        
    # An example of searching the D-optimal exact design by DE for the Hunt-Bowman model for N = 10. 
    deinfo_exact <- getDEInfo(nPop = 32, maxIter = 500, deType = "rand-to-best-1")
    deinfo_approx <- getDEInfo(nPop = 64, maxIter = 1000, deType = "rand-to-best-1")
        
    hb_de_res <- hormesis_de(model = "HuntBowman", criterion = "D", nPoints = 10, parms = hb_par, 
                            deinfo_exact = deinfo_exact, deinfo_approx = deinfo_approx, 
                            upper = 0.15, lower = 0, nRep = 5)
            
