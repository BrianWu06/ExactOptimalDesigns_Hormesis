#devtools::install_github("PingYangChen/globpso")

library(shiny)
library(globpso)
library(waiter)

# Define UI for application that draws a histogram
ui <- fluidPage(

  withMathJax(), 
  tabsetPanel(
    selected = "User Manual", 
    type = "tabs", 
    id = "mainpanel", 
    
    tabPanel(
      "User Manual",
      
      tags$h2("Usage"), 
      tags$p("This app is for finding the optimal exact designs via the particle swarm 
      optimization(PSO) algorithm presented in \"Exact Optimal Designs for Small Studies in 
      Toxicology with Applications to Hormesis\".
      We support 5 different optimal exact designs, the D-optimal exact design for Hunt-Bowman 
      model, the \\(\\tau\\) and h-optimal exact designs for exp-log model, and the D-optimal exact 
      design for quadratic and cubic logistic models."),
      tags$p("Users are allow to adjust the PSO parameters and the model parameters. To start
      the searching process, just click the \"Start Searching\" button, and the PSO algorithm 
      will be triggered and search for the optimal exact designs of the corresponding model. 
      As the searching process terminates."), 
      
      tags$h2("Parameters users are allowed to adjust."), 
      tags$div(
        tags$b("Number of Particles:"), 
        "How many particles are used for searching in the PSO algorithm. 
        The default is the suggested value in our paper for each design.",
        tags$br(), 
        tags$b("Number of Iterations:"), 
        "How many iterations does the PSO algorithim run. 
        The default is the suggested value in our paper for each design.", 
        tags$br(), 
        tags$b("Number of Replications:"), 
        "How many PSO should be run to find the optimal exact design. 
         Since PSO can sometimes be trapped in the local minima, to run through more times 
         of PSO and pick the one with the best result can usually provide a more reliable 
         result. The default is 1, users can increase the value up to 10.", 
        tags$br(), 
        tags$b("Number of Design Points:"), 
        "The number of design points for the current design. 
         The default is set to the minima design points needed for the design. The upper
         bound is set to 10. For more number of design points, please refer to the following
         github link and feel free to adjust the code.", 
        tags$br(), 
        tags$a(href="https://github.com/brian32609/Hormesis_Optimal_Exact_Design_PSO", 
               "https://github.com/brian32609/Hormesis_Optimal_Exact_Design_PSO")
        
      ),
      
      tags$h2("Particle Swarm Optimization"), 
      tags$p("Particle swarm optimization is a population-based metaheuristic algorithim that 
             aims to find the optimal solutions that minimize the target function. A swarm is 
             a population of n particles, indexed as \\(p_1,...,p_n\\). Each particle will
             search through the entire solution space iteratively with an initial velocity. 
             In each iteration, the particles will update their velocity \\(v_i\\) and 
             position \\(x_i\\) by the following formulas."),
      tags$p("$$
v^{t+1}_i = wv^t_i+c_1\\beta_1(p_i-x^t_i)+c_2\\beta_2(g-x^t_i)
$$"),
      tags$p("$$
x^{t+1}_i = x^t_i + v^{t+1}_i
$$"),
      tags$p("\\(v_i^{t}\\): The velocity of the \\(i^{th}\\) particle at the \\(t^{th}\\) iteration."),
      tags$p("\\(x_i^{t}\\): The position of the \\(i^{th}\\) particle at the \\(t^{th}\\) iteration."),
      tags$p("\\(p_i\\): The personal best solution of the \\(i^{th}\\) particle."), 
      tags$p("\\(\\beta_1\\) and \\(\\beta_2\\): Two random variables drawn independently
             and uniformly from [0,1]"),
      tags$p("\\(\\omega\\): The inertia wight which represtns how active the particle is. 
             Can be a decreasing function or a fixed constant. A uniformly decreasing from 
             0.9 to 0.4 in this app."),
      tags$p("\\(c_1\\) and \\(c_2\\): Two tuning parameters called cognitive coefficient 
             and social coefficient respectively. Can be viewed as how much information
             the particle used to update its velocity from the personal best and the global
             best. \\(c_1=c_2=2.05\\) in this app."),
      tags$p("Particle swarm optimization has been shown to solve exact optimal design 
             problems efficiently in some previous tasks. Hence we will use the algorithim
             to generate the optimal exact designs in this app."), 
      
      tags$h2("Hormesis"), 
      tags$p("Hormesis is a dose response relationship that has low-dose positive effect 
      and high-dose negative effect, and the threshold dose level is defined as the maximum 
      nonzero exposure level below which no adverse events above background response occur.
      Which can be described by the following equation."),
      tags$p("$$
\\tau=\\tau(\\theta)=max \\{ d\\in\\Omega:\\mu(d,\\theta)\\leq\\mu(d,\\theta) \\}
$$"),
      tags$p("\\(\\mu(d,\\theta)\\) is the mean response function that characterize the 
      overall dose-response relationship, and \\(\\Omega = [0, \\hat{d}]\\) is the prespecified 
      dose interval."),
      
      tags$h2("D-optimal Exact Design for Hunt-Bowman Model"), 
      tags$p("Hunt and Bowman have charaterized the dose-response model by a piecewise quadratic
             logistic model as following."),
      tags$p("$$
 \\mu(d)=\\begin{cases} 
  c_1d^2+c_2d+\\kappa, & 0 \\leq d\\leq\\tau \\\\
  \\frac{1}{1+e^{\\beta_0-\\beta_1(d-\\tau)}}, & \\tau<d 
  \\end{cases}
 $$"),
      tags$p("Due to the constraints of hormesis threshold, we have 
             \\(\\kappa=\\frac{1}{1+e^{\\beta_0}}\\) and 
             \\(c_2=-c_1 \\tau\\). Hence the parameter set of the Hunt-Bowman model is
             \\(\\theta=(c_1,\\tau,\\beta_0,\\beta_1)\\)."),
      tags$p("For estimating the model parameters of the mean function, the D-optimal design 
             is defined as the determinant of the information matrix, which aims to minimize
             the variance of the estimation of model parameters."),
      tags$p("$$
M(\\xi,\\theta) = \\sum_i f(d_i,\\theta)f^T(d_i,\\theta)\\omega_i
$$"),
      tags$p("$$
f(d,\\theta)=\\frac{\\partial}{\\partial\\theta}\\mu(d,\\theta)
$$"),
      tags$p("For exact designs, we will fix the weights \\(\\omega_i=\\frac{1}{N}\\) for 
             all design points. Where N is the numebr of design points."),
      
      tags$h2("\\(\\tau\\)-optimal Exact Design for exp-log Model"), 
      tags$p("Dette et al. proposed a smooth analytic model that dose not involve the 
             threshold dose level parameter \\(\\tau\\), which is a sum of an exponential 
             decay and curve and a sigmoidal curve."),
      tags$p("$$
 \\mu(d) = c_0e^{-c_1d}+\\frac{1}{1+e^{\\beta_0-\\beta_1d}}
 $$"),
      tags$p("For estimating the threshold dose level \\(\\tau\\), the \\(\\tau\\)-optimal
             design is aimed to minimize the variance of \\(\\tau\\) estimation, which can 
             be viewed as a special case of c-optimal criterion:"),
      tags$p("$$
b^T(\\theta)M^{-1}(\\xi,\\theta)b(\\theta)
$$"),
      tags$p("$$
b(\\theta)=\\frac{\\partial}{\\partial\\theta}\\tau(\\theta)
$$"),
      
      
      tags$h2("h-optimal Exact Design for exp-log Model"), 
      tags$p("For detecting the existence of hormesis, Dette et al. proposed the h-optimal 
             criterion, which can also be treated as a special case of c-optimal criterion:"),
      tags$p("$$
h^T(0,\\theta)M^-1(\\xi,\\theta)h(0,\\theta)
$$"),
      tags$p("$$
h(d,\\theta)=\\frac{\\partial f(d,\\theta)}{\\partial d}
$$"),
      
      tags$h2("D-optimal Exact Design for Simple Logistic Model"), 
      tags$p("Logistic models are widely used to model the probability of a binary type 
             response variable. That is, whether an event will occur or not. The simple 
             logistic model with one predictor variable is proposed as following."),
      tags$p("$$
E(y) = \\frac{e^{\\alpha + \\beta x}}{1 + e^{\\alpha + \\beta x}}
$$"),
      tags$p("The D-optimal exact design for simpel logistic model maximize the determinent of 
             the information matrix and provides the smallest volume of the asymptotic 
             confidence region of the parameter set \\(\\theta=(\\alpha,\\beta)\\)"),
      
      tags$h2("D-optimal Exact Design for Quadratic Logistic Model"), 
      tags$p("Similar as the simple logistic model, the quadratic logistic model only adds the 
             quadratic term of the predictor variabls."),
      tags$p("$$
E(y) = \\frac{e^{\\alpha + \\beta_1 x + \\beta_2 x^2}}{1 + e^{\\alpha + \\beta_1 x + \\beta_2 x^2}}
$$"),
      tags$p("Which gives the parameter set \\(\\alpha,\\beta_1,\\beta_2\\). The D-optimal 
             exact design also aims to mazimize the determinent of the paramter set."),
      
    ),
    
    tabPanel(
      "Find Exact Optimal Designs", 
      
      tabsetPanel(
        selected = "D-optimal Exact Design for Hunt-Bowman Model", 
        type = "pills", 
        id = "HBD", 
        tabPanel("D-optimal Exact Design for Hunt-Bowman Model", 
                 sidebarLayout(
                   sidebarPanel(
                     tags$h3("PSO parameters"), 
                     numericInput("hb_swarm", 
                                  "Number of Swarms", 
                                  value = 64, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1), 
                     numericInput("hb_iter", 
                                  "Number of Iterations", 
                                  value = 1000, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1
                     ), 
                     numericInput("hb_rep", 
                                  "Number of Replications", 
                                  value = 1, 
                                  min = 1, 
                                  max = 10, 
                                  step = 1
                     ),
                     tags$h3("Design Parameters"), 
                     numericInput("hb_dim", 
                                  "Number of Design Points", 
                                  value = 4, 
                                  min = 1, 
                                  max = 15, 
                                  step = 1
                     ),
                   ), 
                   mainPanel(
                     tags$h3("D-optimal Exact Design for Hunt-Bowman Model"), 
                     tags$p("The Hunt-Bowman model is presented as the following equation."),
                     tags$p("$$
 \\mu(d)=\\begin{cases} 
  c_1d^2+c_2d+\\kappa, & 0 \\leq d\\leq\\tau \\\\
  \\frac{1}{1+e^{\\beta_0-\\beta_1(d-\\tau)}}, & \\tau<d 
  \\end{cases}
 $$"),
                     tags$p("The D-optimal criterion aims to minimize the variance of the 
                            estimation of the model parameters. It is expected to run 7 seconds
                            under the default setting."),
                     fluidRow(
                       column(8, 
                              tags$h4("Hunt-Bowman Model Parameters"), 
                              numericInput("hb_c1", 
                                           "\\(c_1\\)", 
                                           value = 170, 
                                           min = 1, 
                                           max = Inf), 
                              numericInput("hb_tau", 
                                           "\\(\\tau\\)", 
                                           value = 0.04, 
                                           min = 0, 
                                           max = 0.15), 
                              numericInput("hb_b0", 
                                           "\\(\\beta_0\\)", 
                                           value = 1.46, 
                                           min = 0, 
                                           max = Inf), 
                              numericInput("hb_b1", 
                                           "\\(\\beta_1\\)", 
                                           value = 40, 
                                           min = 0, 
                                           max = Inf), 
                              
                       )
                     ), 
                     use_waiter(),
                     actionButton("hb_pso", "Start Searching"), 
                     tagAppendAttributes(verbatimTextOutput("hb_out"), style = "height:100px;")
                   )
                 )
        ), 
        
        tabPanel("\\(\\tau\\)-Optiaml Exact Design for exp-log Model", 
                 sidebarLayout(
                   sidebarPanel(
                     tags$h3("PSO parameters"), 
                     numericInput("eltau_swarm", 
                                  "Number of Swarms", 
                                  value = 128, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1), 
                     numericInput("eltau_iter", 
                                  "Number of Iterations", 
                                  value = 2000, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1
                     ), 
                     numericInput("eltau_rep", 
                                  "Number of Replications", 
                                  value = 1, 
                                  min = 1, 
                                  max = 10, 
                                  step = 1
                     ), 
                     tags$h3("Design Parameters"), 
                     numericInput("eltau_dim", 
                                  "Number of Design Points", 
                                  value = 2, 
                                  min = 1, 
                                  max = 15, 
                                  step = 1
                     ),
                   ), 
                   mainPanel(
                     tags$h3("\\(\\tau\\)-optimal Exact Design for exp-log Model"), 
                     tags$p("The exp-log model is presented as the following equation:"),
                     tags$p("$$
 \\mu(d) = c_0e^{-c_1d}+\\frac{1}{1+e^{\\beta_0-\\beta_1d}}
 $$"),
                     tags$p("The \\(\\tau\\)-optimal criterion aims to minimize the variance of the 
                            estimation of the threshold dose level \\(\\tau\\). It is expected 
                            to run 35 seconds under the default setting."),
                     fluidRow(
                       column(8, 
                              tags$h4("exp-log Model Parameters"), 
                              numericInput("eltau_c0", 
                                           "\\(c_0\\)", 
                                           value = 0.15, 
                                           min = 0, 
                                           max = Inf), 
                              numericInput("eltau_c1", 
                                           "\\(c_1\\)", 
                                           value = 89, 
                                           min = 0, 
                                           max = Inf), 
                              numericInput("eltau_b0", 
                                           "\\(\\beta_0\\)", 
                                           value = 3.2, 
                                           min = 0, 
                                           max = Inf), 
                              numericInput("eltau_b1", 
                                           "\\(\\beta_1\\)", 
                                           value = 41, 
                                           min = 0, 
                                           max = Inf), 
                              
                       )
                     ), 
                     use_waiter(),
                     actionButton("eltau_pso", "Start Searching"), 
                     tagAppendAttributes(verbatimTextOutput("eltau_out"), style = "height:100px;")
                   )
                 )
        ), 
        
        tabPanel("h-Optiaml Exact Design for exp-log Model.", 
                 sidebarLayout(
                   sidebarPanel(
                     tags$h3("PSO parameters"), 
                     numericInput("elh_swarm", 
                                  "Number of Swarms", 
                                  value = 128, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1), 
                     numericInput("elh_iter", 
                                  "Number of Iterations", 
                                  value = 2000, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1
                     ),
                     numericInput("elh_rep", 
                                  "Number of Replications", 
                                  value = 1, 
                                  min = 1, 
                                  max = 10, 
                                  step = 1
                     ), 
                     tags$h3("Design Parameters"), 
                     numericInput("elh_dim", 
                                  "Number of Design Points", 
                                  value = 4, 
                                  min = 1, 
                                  max = 15, 
                                  step = 1
                     ),
                   ), 
                   mainPanel(
                     tags$h3("exp-log Model h-optimal Exact Design"), 
                    tags$p("The exp-log model is presented as the following equation:"), 
                     tags$p("$$
 \\mu(d) = c_0e^{-c_1d}+\\frac{1}{1+e^{\\beta_0-\\beta_1d}}
 $$"),
                    tags$p("The h-optimal criterion aims to detect the existence of hormesis. 
                           It is expected to run 40 seconds under the default setting."), 
                     fluidRow(
                       column(8, 
                              tags$h4("exp-log Model Parameters"), 
                              numericInput("elh_c0", 
                                           "\\(c_0\\)", 
                                           value = 0.15, 
                                           min = 0, 
                                           max = Inf), 
                              numericInput("elh_c1", 
                                           "\\(c_1\\)", 
                                           value = 89, 
                                           min = 0, 
                                           max = Inf), 
                              numericInput("elh_b0", 
                                           "\\(\\beta_0\\)", 
                                           value = 3.2, 
                                           min = 0, 
                                           max = Inf), 
                              numericInput("elh_b1", 
                                           "\\(\\beta_1\\)", 
                                           value = 41, 
                                           min = 0, 
                                           max = Inf), 
                              
                       )
                     ), 
                     use_waiter(),
                     actionButton("elh_pso", "Start Searching"), 
                     tagAppendAttributes(verbatimTextOutput("elh_out"), style = "height:100px;")
                   )
                 )
        ), 
        
        tabPanel("D-Optiaml Exact Design for Simple Logistic Model", 
                 sidebarLayout(
                   sidebarPanel(
                     tags$h3("PSO parameters"), 
                     numericInput("logistic_swarm", 
                                  "Number of Swarms", 
                                  value = 64, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1), 
                     numericInput("logistic_iter", 
                                  "Number of Iterations", 
                                  value = 1000, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1
                     ),
                     numericInput("logistic_rep", 
                                  "Number of Replications", 
                                  value = 1, 
                                  min = 1, 
                                  max = 10, 
                                  step = 1
                     ), 
                     tags$h3("Design Parameters"), 
                     numericInput("logistic_dim", 
                                  "Number of Design Points", 
                                  value = 2, 
                                  min = 1, 
                                  max = 15, 
                                  step = 1
                     ),
                   ), 
                   mainPanel(
                     tags$h3("D-optimal Exact Design for Simple Logistic Model"), 
                     tags$p("The simple logistic model is presented as the following equation:"),
                     tags$p("$$
 E(y) = \\frac{e^{\\alpha + \\beta x}}{1 + e^{\\alpha + \\beta x}}
 $$"),
                     tags$p("The D-optimal criterion aims to minimize the variance of the 
                            estimation of the model parameters. It is expected to run 2 seconds
                            under the default setting."), 
                     fluidRow(
                       column(8, 
                              tags$h4("Simple Logistic Model Parameters"), 
                              numericInput("logistic_a", 
                                           "\\(\\alpha\\)", 
                                           value = 2, 
                                           min = -Inf, 
                                           max = Inf), 
                              numericInput("logistic_b", 
                                           "\\(\\beta\\)", 
                                           value = 1, 
                                           min = -Inf, 
                                           max = Inf), 
                              
                       )
                     ), 
                     use_waiter(),
                     actionButton("logistic_pso", "Start Searching"), 
                     tagAppendAttributes(verbatimTextOutput("logistic_out"), style = "height:100px;")
                   )
                 )
        ), 
        
        tabPanel("D-Optiaml Exact Design for Quadratic Logistic Model.", 
                 sidebarLayout(
                   sidebarPanel(
                     tags$h3("PSO parameters"), 
                     numericInput("qlogistic_swarm", 
                                  "Number of Swarms", 
                                  value = 64, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1), 
                     numericInput("qlogistic_iter", 
                                  "Number of Iterations", 
                                  value = 1000, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1
                     ),
                     numericInput("qlogistic_rep", 
                                  "Number of Replications", 
                                  value = 1, 
                                  min = 1, 
                                  max = 10, 
                                  step = 1
                     ), 
                     tags$h3("Design Parameters"), 
                     numericInput("qlogistic_dim", 
                                  "Number of Design Points", 
                                  value = 4, 
                                  min = 1, 
                                  max = 15, 
                                  step = 1
                     ),
                   ), 
                   mainPanel(
                     tags$h3("D-optimal Exact Design for Quadratic Logistic Model."), 
                     tags$p("The quadratic logistic model is presented as the following equation:"),
                     tags$p("$$
 E(y) = \\frac{e^{\\alpha + \\beta_1 x + \\beta_2 x^2}}{1 + e^{\\alpha + \\beta_1 x + \\beta_2 x^2}}
 $$"),
                     tags$p("The D-optimal criterion aims to minimize the variance of the 
                            estimation of the model parameters. It is expected to run 2 seconds
                            under the default setting."),
                     fluidRow(
                       column(8, 
                              tags$h4("Quadratic Logistic Model Parameters"), 
                              numericInput("qlogistic_a", 
                                           "\\(\\alpha\\)", 
                                           value = -3, 
                                           min = -Inf, 
                                           max = Inf), 
                              numericInput("qlogistic_b1", 
                                           "\\(\\beta_1\\)", 
                                           value = 0, 
                                           min = -Inf, 
                                           max = Inf), 
                              numericInput("qlogistic_b2", 
                                           "\\(\\beta_2\\)", 
                                           value = -1, 
                                           min = -Inf, 
                                           max = Inf)
                              
                       )
                     ), 
                     use_waiter(),
                     actionButton("qlogistic_pso", "Start Searching"), 
                     tagAppendAttributes(verbatimTextOutput("qlogistic_out"), style = "height:100px;")
                   )
                 )
        )
        
      )
      
    )
  )
  
    
  )


# Define server logic required to draw a histogram
server <- function(input, output) {
  values <- reactiveValues()
  values$hbd <- list(eff = numeric(), points = c())
  values$eltau <- list(eff = numeric(), points = c())
  values$elh <- list(eff = numeric(), points = c())
  values$logistic <- list(val = c(), points = c())
  values$qlogistic <- list(val = c(), points = c())
  
  observeEvent(input$hb_pso, {
    waiter1 <- waiter::Waiter$new(
      id = "hb_out", 
      html = spin_circle(), 
      color = "lightgrey"
    )
    waiter1$show()
    on.exit(waiter1$hide())
    
    c1 = input$hb_c1
    tau = input$hb_tau
    b0 = input$hb_b0
    b1 = input$hb_b1
    nswarm = input$hb_swarm
    iter = input$hb_iter
    ndp = input$hb_dim
    nrep = input$hb_rep
    
    hb_par <- hb_parms(c1 = c1, tau = tau, b0 = b0, b1 = b1)
    psoinfo_hb <- psoinfo_setting(nSwarms = nswarm, Iters = iter)
    hb_res <- hb_pso_rep(nRep = nrep, nPoints = ndp, parms = hb_par, psoinfo = psoinfo_hb)
    
    values$hbd$eff = hb_res$result$best_eff
    values$hbd$points = hb_res$result$best_points
  }
               )
  
  output$hb_out <- renderPrint({
    cat("Design Points:", values$hbd$points, "\n")
  })
  
  
  observeEvent(input$eltau_pso, {
    waiter2 <- waiter::Waiter$new(
      id = "eltau_out", 
      html = spin_circle(), 
      color = "lightgrey"
    )
    waiter2$show()
    on.exit(waiter2$hide())
    
    c0 = input$eltau_c0
    c1 = input$eltau_c1
    b0 = input$eltau_b0
    b1 = input$eltau_b1
    nswarm = input$eltau_swarm
    iter = input$eltau_iter
    ndp = input$eltau_dim
    nrep = input$eltau_rep
    
    el_par <- exp_log_params(c0 = c0, c1 = c1, b0 = b0, b1 = b1)
    psoinfo_el <- psoinfo_setting(nSwarms = nswarm, Iters = iter)
    eltau_res <- exp_log_tau_pso_rep(nRep = nrep, nPoints = ndp, parms = el_par, psoinfo = psoinfo_el)
    
    values$eltau$eff = eltau_res$result$best_eff
    values$eltau$points = eltau_res$result$best_points
  }
  )
  
  output$eltau_out <- renderPrint({
    cat("Design Points:", values$eltau$points, "\n")
  })
  
  observeEvent(input$elh_pso, {
    waiter3 <- waiter::Waiter$new(
      id = "elh_out", 
      html = spin_circle(), 
      color = "lightgrey"
    )
    waiter3$show()
    on.exit(waiter3$hide())
    
    c0 = input$elh_c0
    c1 = input$elh_c1
    b0 = input$elh_b0
    b1 = input$elh_b1
    nswarm = input$elh_swarm
    iter = input$elh_iter
    ndp = input$elh_dim
    nrep = input$elh_rep
    
    el_par <- exp_log_params(c0 = c0, c1 = c1, b0 = b0, b1 = b1)
    psoinfo_el <- psoinfo_setting(nSwarms = nswarm, Iters = iter)
    elh_res <- exp_log_h_pso_rep(nRep = nrep, nPoints = ndp, parms = el_par, psoinfo = psoinfo_el)
    
    values$elh$eff = elh_res$result$best_eff
    values$elh$points = elh_res$result$best_points
  }
  )
  
  output$elh_out <- renderPrint({
    cat("Design Points:", values$elh$points, "\n")
  })
  
  observeEvent(input$logistic_pso, {
    waiter4 <- waiter::Waiter$new(
      id = "logistic_out", 
      html = spin_circle(), 
      color = "lightgrey"
    )
    waiter4$show()
    on.exit(waiter4$hide())
    
    alpha = input$logistic_a
    beta = input$logistic_b
    nswarm = input$logistic_swarm
    iter = input$logistic_iter
    ndp = input$logistic_dim
    nrep = input$logistic_rep
    
    logistic_par <- c(alpha, beta)
    psoinfo_logistic <- psoinfo_setting(nSwarms = nswarm, Iters = iter)
    logistic_res <- logistic_pso_rep(nRep = nrep, nPoints = ndp, parms = logistic_par, psoinfo = psoinfo_logistic)
    
    values$logistic$val = logistic_res$result$best_val
    print(values$logistic$val)
    values$logistic$points = logistic_res$result$best_points
  }
  )
  
  output$logistic_out <- renderPrint({
    cat("Design Points:", values$logistic$points, "\n")
  })
  
  observeEvent(input$qlogistic_pso, {
    waiter5 <- waiter::Waiter$new(
      id = "qlogistic_out", 
      html = spin_circle(), 
      color = "lightgrey"
    )
    waiter5$show()
    on.exit(waiter5$hide())
    
    alpha = input$qlogistic_a
    beta1 = input$qlogistic_b1
    beta2 = input$qlogistic_b2
    nswarm = input$qlogistic_swarm
    iter = input$qlogistic_iter
    ndp = input$qlogistic_dim
    nrep = input$qlogistic_rep
    
    ql_par <- c(alpha, beta1, beta2)
    psoinfo_ql <- psoinfo_setting(nSwarms = nswarm, Iters = iter)
    ql_res <- qlogistic_pso_rep(nRep = nrep, nPoints = ndp, parms = ql_par, psoinfo = psoinfo_ql)

    values$qlogistic$val = ql_res$result$best_val
    values$qlogistic$points = ql_res$result$best_points
  }
  )
  
  output$qlogistic_out <- renderPrint({
    cat("Design Points:", values$qlogistic$points, "\n")
  })
}



## Initalize PSO settings

psoinfo_setting <- function(nSwarms = 64, Iters = 1000){
  getPSOInfo(nSwarm = nSwarms, maxIter = Iters, w0 = 0.9, w1 = 0.4, w_var = 1)
}

## Hunt-Bowman model

# Hunt-Bowman parameters
hb_parms <- function(c1, tau, b0, b1){
  c(c1, tau, b0, b1)
}

# Hunt-Bowman information matrix of 1 point
hb_mat <- function(d, c1, tau, b0, b1){
  
  if (d <= tau){
    f1 <- d^2 - tau * d
    f2 <- -c1 * d
    f3 <- -exp(b0) / (1 + exp(b0))^2
    f4 <- 0
    f <- matrix(c(f1, f2, f3, f4))
    mat <- f %*% t(f)
  } 
  else if (d > tau){
    f1 <- 0
    f2 <- -b1 * exp(b0 - b1 * (d - tau)) / (1 + exp(b0 - b1 * (d - tau)))^2
    f3 <- -exp(b0 - b1 * (d - tau)) / (1 + exp(b0 - b1 * (d - tau)))^2
    f4 <- (d - tau) * exp(b0 - b1 * (d - tau)) / (1 + exp(b0 - b1 * (d - tau)))^2
    f <- matrix(c(f1, f2, f3, f4))
    mat <- f %*% t(f)
  }
  
  mat
}

# Hunt-Bowman d-optimality criterion value efficiency
hb_doptimal <- function(d, loc){
  
  # Hunt-Bowman parameters
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  eff_base <- 9.938331e-14 #Criterion value of 4-point uniform design
  
  # Number of experiment points
  d <- c(0, d)
  w <- length(d)
  
  # Evaluate d-optimality criterion value
  mat_list <- lapply(d, function(x) 1/w * hb_mat(x, c1, tau, b0, b1))
  inf_mat <- Reduce("+", mat_list)
  - (det(inf_mat) / eff_base)
}

# Find the exact design points that minimizes the d-optimality efficiency
# of hunt-bowman model by PSO
hb_pso <- function(nPoints, parms, psoinfo){
  
  # set the lower and upper bounds for PSO
  lb <- rep(0, nPoints-1)
  ub <- rep(0.15, nPoints-1)
  
  pso_res <- globpso(objFunc = hb_doptimal, 
                     lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, 
                     loc = parms)
  
  pso_res$val <- (-1 * pso_res$val)^(1/nPoints) |> round(4)
  pso_res$par <- c(0, pso_res$par) |> sort() |> round(4)
  pso_res$history <- pso_res$history * -1
  pso_res
}

# Replicate m PSO results of Hunt-Bowman model d-optimality
hb_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo){
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- hb_pso(nPoints, parms, psoinfo))
  
  hb_list <- list(efficiency = c(), design_points = c(), history = c(), 
                  result = list(best_eff = 0, best_points = c()))
  
  # Efficiency of each replication.
  hb_list$efficiency <- sapply(1:nRep, function(x) pso_results[[x]]$val) 
  # Design points found in each replication.
  hb_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par) 
  # Best result of each iteration in each replication. (For finding the appropriate nRep.)
  hb_list$history <- sapply(1:nRep, function(x) pso_results[[x]]$history) 
  
  best_idx <- which.max(hb_list$efficiency)
  # Best efficiency among all replications.
  hb_list$result$best_eff <- hb_list$efficiency[best_idx]
  # Design points correspond to the best efficiency.
  hb_list$result$best_points <- hb_list$design_points[,best_idx]
  
  hb_list
}

## exp-log model

# exp-log model parameters
exp_log_params <- function(c0, c1, b0, b1){
  c(c0, c1, b0, b1)
}

# exp-log model information matrix
exp_log_f <- function(d, c0, c1, b0, b1){
  f1 <- exp(-c1 * d)
  f2 <- -c0 * d * exp(-c1 * d)
  f3 <- -exp(b0 - b1 * d) / (1 + exp(b0 - b1 * d))^2
  f4 <- d * exp(b0 - b1 * d) / (1 + exp(b0 - b1 * d))^2
  f <- matrix(c(f1, f2, f3, f4))
  f %*% t(f)
}

exp_log_mat <- function(d, c0, c1, b0, b1){
  #d <- c(0, d)
  n <- length(d)
  mat_list <- lapply(d, function(x) 1/n * exp_log_f(x, c0, c1, b0, b1))
  
  #print(mat_list)
  M <- Reduce("+", mat_list)
  diag(M) <- diag(M) + 1e-10
  M
}

# exp-log model h-optimality criterion value
exp_log_hoptimal <- function(d, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  pen = 999999999 # Penalty value
  eff_base = 243738.4 # Criterion value of 4-point locally optimal design
  
  n <- length(d)
  h1 <- -c1
  h2 <- -c0
  h3 <- b1 * exp(b0) * (1 - exp(b0)) / (exp(b0) + 1)^3
  h4 <- exp(b0) / (exp(b0) + 1)^2 
  h <- matrix(c(h1, h2, h3, h4))
  
  M <- exp_log_mat(d, c0, c1, b0, b1)
  if (rcond(M) < 2.220446e-16){
    res = pen
  }
  else {
    M_inv <- solve(M)
    res = t(h) %*% M_inv %*% h
  }
  
  res / eff_base
}

# Find the exact design points that minimizes the h-optimality efficiency
# of exp-log model by PSO
exp_log_pso <- function(nPoints, parms, psoinfo){
  
  lb <- rep(0, nPoints)
  ub <- rep(0.15, nPoints)
  
  pso_res <- globpso(objFunc = exp_log_hoptimal, 
                     lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, 
                     loc = parms)
  
  pso_res$val <- (1 / pso_res$val)^(1/nPoints)  |> round(4)
  pso_res$par <- c(pso_res$par) |> sort() |> round(4)
  pso_res$history <- 1 / pso_res$history
  pso_res
}

# Replicate m PSO results of exp-log model h-optimality
exp_log_h_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo){
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- exp_log_pso(nPoints, parms, psoinfo))
  
  exp_log_list <- list(efficiency = c(), design_points = c(), history = c(), 
                       result = list(best_eff = 0, best_points = c()))
  
  # Efficiency of each replication.
  exp_log_list$efficiency <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  # Design points found in each replication.
  exp_log_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par)
  # Best result of each iteration in each replication. (For finding the appropriate nRep.)
  exp_log_list$history <- sapply(1:nRep, function(x) pso_results[[x]]$history)
  
  best_idx <- which.max(exp_log_list$efficiency)
  # Best efficiency among all replications.
  exp_log_list$result$best_eff <- exp_log_list$efficiency[best_idx]
  # Design points correspond to the best efficiency.
  exp_log_list$result$best_points <- exp_log_list$design_points[,best_idx]
  
  exp_log_list
}


## tau-optimal 

tau_func <- function(d, c0, c1, b0, b1){
  c0 * exp(-c1*d) + 1 / (1 + exp(b0 - b1*d)) - (c0 + 1 / (1 + exp(b0)))
}

exp_log_b <- function(d, c0, c1, b0, b1){
  h <- ((b1 * exp(b1*d + b0)) / (exp(b1*d) + exp(b0))^2) - c0 * c1 * exp(-c1*d)
  h1 <- -1*(exp(-c1*d) - 1) / h
  h2 <- -1*(-c0 * d * exp(-c1*d)) / h
  h3 <- -1*(exp(b0) * (1 / (exp(b0) + 1)^2 - exp(b1*d) / (exp(b1*d) + exp(b0))^2 )) / h
  h4 <- -1*(d * exp(b1*d + b0) / (exp(b1 * d) + exp(b0))^2) / h
  
  b <- matrix(c(h1, h2, h3, h4))
}

exp_log_loc_mat <- function(d, w, c0, c1, b0, b1){
  n <- length(d)
  mat_list <- lapply(1:n, function(x) w[x] * exp_log_f(d[x], c0, c1, b0, b1))
  
  #print(mat_list)
  M <- Reduce("+", mat_list)
  diag(M) <- diag(M) + 1e-10
  M
}

exp_log_loc_tauoptimal <- function(x, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  tau = loc[5]
  pen = 999999999
  
  n <- round(length(x) / 2)
  d <- x[1:n]
  w <- x[-(1:n)]
  w[n] = 1 - sum(w)
  
  #tau <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, c0 = c0, c1 = c1, b0 = b0, b1 = b1)$root
  #print(tau)
  b <- exp_log_b(tau, c0, c1, b0, b1)
  
  M <- exp_log_loc_mat(d, w, c0, c1, b0, b1)
  if (rcond(M) < 2.220446e-16){
    res = pen
  } 
  else if (w[n] <= 0){
    res = pen
  }
  else {
    M_inv <- solve(M)
    #print(M_inv)
    #print(b)
    res = t(b) %*% M_inv %*% b
  }
  
  res
}

exp_log_loc_tau_pso <- function(nPoints, parms, psoinfo){
  
  lb <- c(rep(0, nPoints), rep(0, nPoints-1))
  ub <- c(rep(0.15, nPoints), rep(1, nPoints-1))
  tau <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, 
                 c0 = parms[1], c1 = parms[2], b0 = parms[3], b1 = parms[4])$root
  
  pso_res <- globpso(objFunc = exp_log_loc_tauoptimal, 
                     lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, 
                     loc = c(parms, tau))
  
  pso_res$val <- pso_res$val
  pso_res$points <- pso_res$par[1:nPoints] |> round(4)
  pso_res$weight <- c(pso_res$par[-c(1:nPoints)], 1 - sum(pso_res$par[-c(1:nPoints)])) |> round(4)
  pso_res
}

exp_log_loc_tau_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo){
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- exp_log_loc_tau_pso(nPoints, parms, psoinfo))
  
  exp_log_list <- list(design_points = c(), weight = c(), val = c())
  
  # Efficiency of each replication.
  exp_log_list$weight <- sapply(1:nRep, function(x) pso_results[[x]]$weight)
  # Design points found in each replication.
  exp_log_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$points)
  exp_log_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  
  exp_log_list
}

exp_log_tauoptimal <- function(d, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  tau = loc[5]
  pen = 999999999
  
  n <- length(d)
  w <- rep(1/n, n)
  
  b <- exp_log_b(tau, c0, c1, b0, b1)
  
  M <- exp_log_loc_mat(d, w, c0, c1, b0, b1)
  if (rcond(M) < 2.220446e-16){
    res = pen
  }
  else {
    M_inv <- solve(M)
    #print(M_inv)
    #print(b)
    res = t(b) %*% M_inv %*% b
  }
  
  #print(tau)
  res / 0.1157493
}

exp_log_tau_pso <- function(nPoints, parms, psoinfo){
  
  lb <- c(rep(0, nPoints))
  ub <- c(rep(0.15, nPoints))
  tau <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, 
                 c0 = parms[1], c1 = parms[2], b0 = parms[3], b1 = parms[4])$root
  
  pso_res <- globpso(objFunc = exp_log_tauoptimal, 
                     lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, 
                     loc = c(parms, tau))
  
  pso_res$val <- (1 / pso_res$val)^(1/nPoints) |> round(4)
  pso_res$par <- c(pso_res$par) |> sort() |> round(4)
  pso_res$history <- 1 / pso_res$history
  pso_res
}

exp_log_tau_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo){
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- exp_log_tau_pso(nPoints, parms, psoinfo))
  
  exp_log_list <- list(design_points = c(), weight = c(), val = c())
  
  # Efficiency of each replication.
  exp_log_list$efficiency <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  # Design points found in each replication.
  exp_log_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par)
  # Best result of each iteration in each replication. (For finding the appropriate nRep.)
  exp_log_list$history <- sapply(1:nRep, function(x) pso_results[[x]]$history)
  
  best_idx <- which.max(exp_log_list$efficiency)
  # Best efficiency among all replications.
  exp_log_list$result$best_eff <- exp_log_list$efficiency[best_idx]
  # Design points correspond to the best efficiency.
  exp_log_list$result$best_points <- exp_log_list$design_points[,best_idx]
  
  exp_log_list
}

logistic_doptimal <- function(d, loc){
  a <- loc[1]
  b <- loc[2]
  n <- length(d)
  
  theta <- a + b * d
  omega <- exp(theta) / (1 + exp(theta))^2
  
  inf_mat <- matrix(c(sum(omega), sum(d * omega), sum(d * omega), sum(d^2 * omega)), nrow = 2)
  -det(inf_mat) / n^2
}

logistic_pso <- function(npoints, parms, psoinfo){
  lb <- rep(-10, npoints)
  ub <- rep(10, npoints)
  
  pso_res <- globpso(objFunc = logistic_doptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$val <- -1 * pso_res$val
  pso_res
}

logistic_pso_rep <- function(nRep, nPoints = 2, parms, psoinfo){
  
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- logistic_pso(nPoints, parms, psoinfo))
  
  logistic_list <- list(design_points = c(), val = c())
  logistic_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par)
  logistic_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  
  best_idx <- which.max(logistic_list$val)
  logistic_list$result$best_val <- logistic_list$val[best_idx]
  logistic_list$result$best_points <- logistic_list$design_points[, best_idx] |> sort() |> round(4)
  
  logistic_list
}

qlogistic_doptimal <- function(d, loc){
  a <- loc[1]
  b1 <- loc[2]
  b2 <- loc[3]
  n <- length(d)
  
  theta <- a + b1 * d + b2 * d^2
  omega <- exp(theta) / (1 + exp(theta))^2
  
  inf_mat <- matrix(c(sum(omega), sum(d * omega), sum(d^2 * omega), 
                      sum(d * omega), sum(d^2 * omega), sum(d^3 * omega), 
                      sum(d^2 * omega), sum(d^3 * omega), sum(d^4 * omega)), 
                    nrow = 3)
  #print(inf_mat)
  
  -det(inf_mat) / n^3
}

qlogistic_pso <- function(npoints, parms, psoinfo){
  lb <- rep(-10, npoints)
  ub <- rep(10, npoints)
  
  pso_res <- globpso(objFunc = qlogistic_doptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$val <- -1 * pso_res$val
  pso_res
}

qlogistic_pso_rep <- function(nRep, nPoints = 2, parms, psoinfo){
  
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- qlogistic_pso(nPoints, parms, psoinfo))
  
  qlogistic_list <- list(design_points = c(), val = c())
  qlogistic_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par)
  qlogistic_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  
  best_idx <- which.max(qlogistic_list$val)
  qlogistic_list$result$best_val <- qlogistic_list$val[best_idx]
  qlogistic_list$result$best_points <- qlogistic_list$design_points[, best_idx] |> sort() |> round(4)
  
  qlogistic_list
}

# Run the application 
shinyApp(ui = ui, server = server)
