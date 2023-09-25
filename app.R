library(shiny)
library(globpso)
library(waiter)
library(dplyr)

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
      We support various kinds of optimal exact designs, the D, h and \\(\\tau\\)-optimal exact design for 
      Hunt-Bowman model and the exp-log model, 
      and the D-optimal exact design for simple, quadratic and cubic logistic models."),
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
        
        tags$b("Criterion:"), 
        "For Hunt-Bowman and exp-log models, users are allowed to choose the criterion of the optimal exact design. 
        There are three chocies: D, tau and h. The details can be found in the following paragraphs. ", 
        tags$br(), 
        
        tags$a(href="https://github.com/BrianWu06/ExactOptimalDesigns_Hormesis", 
               "https://github.com/BrianWu06/ExactOptimalDesigns_Hormesis")
        
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
             problems efficiently in some previous tasks [1]. Hence we will use the algorithim
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
      
      
      tags$h2("D-optimal Exact Designs"), 
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
      
      
      tags$h2("\\(\\tau\\)-optimal Exact Designs"), 
      tags$p("For estimating the threshold dose level \\(\\tau\\), the \\(\\tau\\)-optimal
             design is aimed to minimize the variance of \\(\\tau\\) estimation, which can 
             be viewed as a special case of c-optimal criterion:"),
      tags$p("$$
b^T(\\theta)M^{-1}(\\xi,\\theta)b(\\theta)
$$"),
      tags$p("$$
b(\\theta)=\\frac{\\partial}{\\partial\\theta}\\tau(\\theta)
$$"),
      
      
      tags$h2("h-optimal Exact Designs"), 
      tags$p("For detecting the existence of hormesis, Dette et al. [3] proposed the h-optimal 
             criterion, which can also be treated as a special case of c-optimal criterion:"),
      tags$p("$$
h^T(0,\\theta)M^-1(\\xi,\\theta)h(0,\\theta)
$$"),
      tags$p("$$
h(d,\\theta)=\\frac{\\partial f(d,\\theta)}{\\partial d}
$$"),
      
      
      tags$h2("Hunt-Bowman Model"), 
      tags$p("Hunt and Bowman [2] have charaterized the dose-response model by a piecewise quadratic
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
      
      
      tags$h2("exp-log Model"), 
      tags$p("Dette et al. [3] proposed a smooth analytic model that dose not involve the 
             threshold dose level parameter \\(\\tau\\), which is a sum of an exponential 
             decay and curve and a sigmoidal curve."),
      tags$p("$$
 \\mu(d) = c_0e^{-c_1d}+\\frac{1}{1+e^{\\beta_0-\\beta_1d}}
 $$"),
      tags$p("The parameter set of the exp-log model is 
             \\(\\theta=(c_0,c_1,\\beta_0,\\beta_1)\\)."),
      
      
      tags$h2("D-optimal Exact Design for Simple Logistic Model"), 
      tags$p("Logistic models are widely used to model the probability of a binary type 
             response variable. That is, whether an event will occur or not. The simple 
             logistic model with one predictor variable is proposed as following."),
      tags$p("$$
E(y) = \\frac{e^{\\alpha + \\beta x}}{1 + e^{\\alpha + \\beta x}}
$$"),
      tags$p("The D-optimal exact design for simpel logistic model [4] maximize the determinent of 
             the information matrix and provides the smallest volume of the asymptotic 
             confidence region of the parameter set \\(\\theta=(\\alpha,\\beta)\\)"),
      
      
      tags$h2("D-optimal Exact Design for Quadratic Logistic Model"), 
      tags$p("Similar as the simple logistic model, the quadratic logistic model only adds the 
             quadratic term of the predictor variabls."),
      tags$p("$$
E(y) = \\frac{e^{\\alpha + \\beta_1 x + \\beta_2 x^2}}{1 + e^{\\alpha + \\beta_1 x + \\beta_2 x^2}}
$$"),
      tags$p("Which gives the parameter set \\(\\alpha,\\beta_1,\\beta_2\\). The D-optimal 
             exact design also aims to mazimize the determinent of the paramter set. [5]"),
      
      
      tags$h2("D-optimal Exact Design for Cubic Logistic Model"), 
      tags$p("The cubic logistic models include the cubic term of the predictor in the 
             quadratic logistic modelss."),
      tags$p("$$
E(y) = \\frac{e^{\\alpha + \\beta_1 x + \\beta_2 x^2 + \\beta_3 x^3}}
{1 + e^{\\alpha + \\beta_1 x + \\beta_2 x^2 + \\beta_3 x^3}}
$$"),
      tags$p("Which gives the parameter set \\(\\alpha,\\beta_1,\\beta_2,\\beta_3\\)."),
      
      
      tags$h2("Reference"),
      tags$div(
        "[1] Ping-Yang Chen, Ray-Bing Chen, and Weng Kee Wong. ", 
        tags$b("Particle swarm optimization for searching efficient experimental designs: a review. "), 
        tags$i("WIREs Computational Statistics, "), 
        "14, 2022.", 
        tags$br(),
        
        "[2]  L. Daniel Hunt and Bowman Dale. ", 
        tags$b(" A parametric model for detecting hormetic effects in developmental toxicity studies. "), 
        tags$i("Risk analysis, "), 
        "24(1):65–72, 2004", 
        tags$br(),
        
        "[3] Holger Dette, Andrey Pepelyshev, and Weng Kee Wong. ", 
        tags$b("Optimal experimental design strategies for detecting hormesis. "), 
        tags$i("Risk Analysis: An International Journal, "), 
        "31(12):1949–1960, 2011.", 
        tags$br(), 
        
        "[4] Salomon Minkin. ", 
        tags$b("Optimal designs for binary data. "), 
        tags$i("Journal of the American Statistical Association, "), 
        "1987.", 
        tags$br(),
        
        "[5] Ellinor Fackle Fornius. ", 
        tags$b("Optimal Design of Experiments for the Quadratic Logistic Model. "), 
        tags$i("PhD thesis, Department of Statistics, Stockholm University, "), 
        "2008.", 
        tags$br(),
      )
      
    ),
    
    tabPanel(
      "Find Optimal Exact Designs", 
      
      tabsetPanel(
        selected = "Optimal Exact Designs for Hunt-Bowman Model", 
        type = "pills", 
        id = "HBD", 
        tabPanel("Optimal Exact Designs for Hunt-Bowman Model", 
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
                                  value = 5, 
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
                     selectInput("hb_criterion", 
                                 "Criterion:", 
                                 c("D" = "hbd", 
                                   "tau" = "hbtau", 
                                   "h" = "hbh")
                       
                     )
                   ), 
                   mainPanel(
                     tags$h3("Exact Designs for Hunt-Bowman Model"), 
                     tags$p("The Hunt-Bowman model is presented as the following equation."),
                     tags$p("$$
 \\mu(d)=\\begin{cases} 
  c_1d^2+c_2d+\\kappa, & 0 \\leq d\\leq\\tau \\\\
  \\frac{1}{1+e^{\\beta_0-\\beta_1(d-\\tau)}}, & \\tau<d 
  \\end{cases}
 $$"),
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
                     tagAppendAttributes(verbatimTextOutput("hb_out"), style = "height:200px;")
                   )
                 )
        ), 
        
        tabPanel("Optiaml Exact Design for exp-log Model", 
                 sidebarLayout(
                   sidebarPanel(
                     tags$h3("PSO parameters"), 
                     numericInput("el_swarm", 
                                  "Number of Swarms", 
                                  value = 128, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1), 
                     numericInput("el_iter", 
                                  "Number of Iterations", 
                                  value = 2000, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1
                     ), 
                     numericInput("el_rep", 
                                  "Number of Replications", 
                                  value = 5, 
                                  min = 1, 
                                  max = 10, 
                                  step = 1
                     ), 
                     tags$h3("Design Parameters"), 
                     numericInput("el_dim", 
                                  "Number of Design Points", 
                                  value = 4, 
                                  min = 1, 
                                  max = 15, 
                                  step = 1
                     ),
                     selectInput("el_criterion", 
                                 "Criterion:", 
                                 c("D" = "eld", 
                                   "tau" = "eltau", 
                                   "h" = "elh")
                                 
                     )
                   ), 
                   mainPanel(
                     tags$h3("Optimal Exact Designs for exp-log Model"), 
                     tags$p("The exp-log model is presented as the following equation:"),
                     tags$p("$$
 \\mu(d) = c_0e^{-c_1d}+\\frac{1}{1+e^{\\beta_0-\\beta_1d}}
 $$"),
                     fluidRow(
                       column(8, 
                              tags$h4("exp-log Model Parameters"), 
                              numericInput("el_c0", 
                                           "\\(c_0\\)", 
                                           value = 0.15, 
                                           min = 0, 
                                           max = Inf), 
                              numericInput("el_c1", 
                                           "\\(c_1\\)", 
                                           value = 89, 
                                           min = 0, 
                                           max = Inf), 
                              numericInput("el_b0", 
                                           "\\(\\beta_0\\)", 
                                           value = 3.2, 
                                           min = 0, 
                                           max = Inf), 
                              numericInput("el_b1", 
                                           "\\(\\beta_1\\)", 
                                           value = 41, 
                                           min = 0, 
                                           max = Inf), 
                              
                       )
                     ), 
                     use_waiter(),
                     actionButton("el_pso", "Start Searching"), 
                     tagAppendAttributes(verbatimTextOutput("el_out"), style = "height:200px;")
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
                                  value = 5, 
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
                     tagAppendAttributes(verbatimTextOutput("logistic_out"), style = "height:200px;")
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
                                  value = 5, 
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
                     tagAppendAttributes(verbatimTextOutput("qlogistic_out"), style = "height:200px;")
                   )
                 )
        ), 
        
        tabPanel("D-Optiaml Exact Design for Cubic Logistic Model.", 
                 sidebarLayout(
                   sidebarPanel(
                     tags$h3("PSO parameters"), 
                     numericInput("clogistic_swarm", 
                                  "Number of Swarms", 
                                  value = 64, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1), 
                     numericInput("clogistic_iter", 
                                  "Number of Iterations", 
                                  value = 1000, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1
                     ),
                     numericInput("clogistic_rep", 
                                  "Number of Replications", 
                                  value = 5, 
                                  min = 1, 
                                  max = 10, 
                                  step = 1
                     ), 
                     tags$h3("Design Parameters"), 
                     numericInput("clogistic_dim", 
                                  "Number of Design Points", 
                                  value = 4, 
                                  min = 1, 
                                  max = 15, 
                                  step = 1
                     ),
                   ), 
                   mainPanel(
                     tags$h3("D-optimal Exact Design for Cubic Logistic Model."), 
                     tags$p("The cubic logistic model is presented as the following equation:"),
                     tags$p("$$
 E(y) = \\frac{e^{\\alpha + \\beta_1 x + \\beta_2 x^2 + \\beta_3 x^3}}
 {1 + e^{\\alpha + \\beta_1 x + \\beta_2 x^2 + \\beta_3 x^3}}
 $$"),
                     fluidRow(
                       column(8, 
                              tags$h4("Cubic Logistic Model Parameters"), 
                              numericInput("clogistic_a", 
                                           "\\(\\alpha\\)", 
                                           value = -3, 
                                           min = -Inf, 
                                           max = Inf), 
                              numericInput("clogistic_b1", 
                                           "\\(\\beta_1\\)", 
                                           value = 0, 
                                           min = -Inf, 
                                           max = Inf), 
                              numericInput("clogistic_b2", 
                                           "\\(\\beta_2\\)", 
                                           value = 0, 
                                           min = -Inf, 
                                           max = Inf), 
                              numericInput("clogistic_b3", 
                                           "\\(\\beta_3\\)", 
                                           value = -1, 
                                           min = -Inf, 
                                           max = Inf)
                              
                       )
                     ), 
                     use_waiter(),
                     actionButton("clogistic_pso", "Start Searching"), 
                     tagAppendAttributes(verbatimTextOutput("clogistic_out"), style = "height:200px;")
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
  values$hb <- list(val = numeric(), eff = numeric(), points = c(), 
                    approximate_design = data.frame(suport_points = c(0), weight = c(0)))
  values$el <- list(val = numeric(), eff = numeric(), points = c(), 
                    approximate_design = data.frame(suport_points = c(0), weight = c(0)))
  values$logistic <- list(val = numeric(), eff = numeric(), points = c(), 
                          approximate_design = data.frame(suport_points = c(0), weight = c(0)))
  values$qlogistic <- list(val = numeric(), eff = numeric(), points = c(), 
                           approximate_design = data.frame(suport_points = c(0), weight = c(0)))
  values$clogistic <- list(val = numeric(), eff = numeric(), points = c(), 
                           approximate_design = data.frame(suport_points = c(0), weight = c(0)))
  
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
    criterion = input$hb_criterion
    
    hb_par <- hb_parms(c1 = c1, tau = tau, b0 = b0, b1 = b1)
    psoinfo_hb <- psoinfo_setting(nSwarms = nswarm, Iters = iter)
    if (criterion == "hbd"){
      hb_res <- hb_doptimal_pso_rep(nRep = nrep, nPoints = ndp, parms = hb_par, psoinfo = psoinfo_hb)
    } else if (criterion == "hbtau"){
      hb_res <- hb_tauoptimal_pso_rep(nRep = nrep, nPoints = ndp, parms = hb_par, psoinfo = psoinfo_hb)
    } else if (criterion == "hbh"){
      hb_res <- hb_hoptimal_pso_rep(nRep = nrep, nPoints = ndp, parms = hb_par, psoinfo = psoinfo_hb)
    }
    
    values$hb$val = hb_res$result$best_val
    values$hb$points = hb_res$result$design_points
    values$hb$eff = hb_res$result$efficiency
    values$hb$approx = hb_res$approximate_design
  }
               )
  
  output$hb_out <- renderPrint({
    cat("Design Points:", values$hb$points, "\n", "\n")
    cat("Efficiency:", values$hb$eff, "\n", "\n")
    cat("Approximate Design:", "\n")
    print(values$hb$approx)
  })
  
  
  observeEvent(input$el_pso, {
    waiter2 <- waiter::Waiter$new(
      id = "el_out", 
      html = spin_circle(), 
      color = "lightgrey"
    )
    waiter2$show()
    on.exit(waiter2$hide())
    
    c0 = input$el_c0
    c1 = input$el_c1
    b0 = input$el_b0
    b1 = input$el_b1
    nswarm = input$el_swarm
    iter = input$el_iter
    ndp = input$el_dim
    nrep = input$el_rep
    criterion = input$el_criterion
    
    el_par <- exp_log_params(c0 = c0, c1 = c1, b0 = b0, b1 = b1)
    psoinfo_el <- psoinfo_setting(nSwarms = nswarm, Iters = iter)
    if (criterion == "eld"){
      el_res <- exp_log_doptimal_pso_rep(nRep = nrep, nPoints = ndp, parms = el_par, psoinfo = psoinfo_el)
    } else if (criterion == "eltau"){
      el_res <- exp_log_tauoptimal_pso_rep(nRep = nrep, nPoints = ndp, parms = el_par, psoinfo = psoinfo_el)
    } else if (criterion == "elh"){
      el_res <- exp_log_hoptimal_pso_rep(nRep = nrep, nPoints = ndp, parms = el_par, psoinfo = psoinfo_el)
    }
    
    values$el$val = el_res$result$best_val
    values$el$points = el_res$result$design_points
    values$el$eff = el_res$result$efficiency
    values$el$approx = el_res$approximate_design
  }
  )
  
  output$el_out <- renderPrint({
    cat("Design Points:", values$el$points, "\n", "\n")
    cat("Efficiency:", values$el$eff, "\n", "\n")
    cat("Approximate Design:", "\n")
    print(values$el$approx)
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
    
    logistic_par <- logistic_params(alpha, beta)
    psoinfo_logistic <- psoinfo_setting(nSwarms = nswarm, Iters = iter)
    logistic_res <- logistic_pso_rep(nRep = nrep, nPoints = ndp, parms = logistic_par, psoinfo = psoinfo_logistic)
    
    values$logistic$val = logistic_res$result$best_val
    values$logistic$points = logistic_res$result$design_points
    values$logistic$eff = logistic_res$result$efficiency
    values$logistic$approx = logistic_res$approximate_design
  }
  )
  
  output$logistic_out <- renderPrint({
    cat("Design Points:", values$logistic$points, "\n", "\n")
    cat("Efficiency:", values$logistic$eff, "\n", "\n")
    cat("Approximate Design:", "\n")
    print(values$logistic$approx)
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
    values$qlogistic$points = ql_res$result$design_points
    values$qlogistic$eff = ql_res$result$efficiency
    values$qlogistic$approx = ql_res$approximate_design
  }
  )
  
  output$qlogistic_out <- renderPrint({
    cat("Design Points:", values$qlogistic$points, "\n", "\n")
    cat("Efficiency:", values$qlogistic$eff, "\n", "\n")
    cat("Approximate Design:", "\n")
    print(values$qlogistic$approx)
  })
  
  observeEvent(input$clogistic_pso, {
    waiter6 <- waiter::Waiter$new(
      id = "clogistic_out", 
      html = spin_circle(), 
      color = "lightgrey"
    )
    waiter6$show()
    on.exit(waiter6$hide())
    
    alpha = input$clogistic_a
    beta1 = input$clogistic_b1
    beta2 = input$clogistic_b2
    beta3 = input$clogistic_b3
    nswarm = input$clogistic_swarm
    iter = input$clogistic_iter
    ndp = input$clogistic_dim
    nrep = input$clogistic_rep
    
    cl_par <- c(alpha, beta1, beta2, beta3)
    psoinfo_cl <- psoinfo_setting(nSwarms = nswarm, Iters = iter)
    cl_res <- clogistic_pso_rep(nRep = nrep, nPoints = ndp, parms = cl_par, psoinfo = psoinfo_cl)
    
    values$clogistic$val = cl_res$result$best_val
    values$clogistic$points = cl_res$result$design_points
    values$clogistic$eff = cl_res$result$efficiency
    values$clogistic$approx = cl_res$approximate_design
  }
  )
  
  output$clogistic_out <- renderPrint({
    cat("Design Points:", values$clogistic$points, "\n", "\n")
    cat("Efficiency:", values$clogistic$eff, "\n", "\n")
    cat("Approximate Design:", "\n")
    print(values$clogistic$approx)
  })
}



## Initalize PSO settings

psoinfo_setting <- function(nSwarms = 64, Iters = 1000){
  getPSOInfo(nSwarm = nSwarms, maxIter = Iters, w0 = 0.9, w1 = 0.4, w_var = 1)
}

### Hunt-Bowman model

## Hunt-Bowman parameters
hb_parms <- function(c1, tau, b0, b1){
  c(c1, tau, b0, b1)
}

# Hunt-Bowman information matrix
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

## Hunt-Bowman D-optimal

hb_doptimal <- function(d, loc){
  
  # Hunt-Bowman parameters
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  
  # Number of experiment points
  d <- c(0, d)
  n <- length(d)
  
  # Evaluate d-optimality criterion value
  mat_list <- lapply(d, function(x) 1/n * hb_mat(x, c1, tau, b0, b1))
  inf_mat <- Reduce("+", mat_list)
  -det(inf_mat)
}

# Find the D-optimal exact design for Hunt-Bowman model
hb_doptimal_pso <- function(nPoints, parms, psoinfo){
  
  # set the lower and upper bounds for PSO
  lb <- rep(0, nPoints-1)
  ub <- rep(0.15, nPoints-1)
  
  #run PSO
  pso_res <- globpso(objFunc = hb_doptimal, 
                     lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, 
                     loc = parms)
  
  pso_res$val <- -1 * pso_res$val
  pso_res$par <- c(0, pso_res$par) |> sort() |> round(4)
  pso_res$history <- pso_res$history * -1
  pso_res
}

# Find D-optimal approximate design for the Hunt-Bowman model
hb_doptimal_approx <- function(parms){
  psoinfo <- psoinfo_setting()
  approx_design <- hb_doptimal_pso(nPoints = 4, parms, psoinfo)
  approx_design
}

# Replicate m PSO results of the D-optimal exact design for Hunt-Bowman model
hb_doptimal_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo){
  
  # The D-optimal approximate design for the Hunt-Bowman model under certain parameter set.
  approx_design <- hb_doptimal_approx(parms)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- hb_doptimal_pso(nPoints, parms, psoinfo))
  
  hb_list <- list(val = c(), design_points = c(),  
                  result = list(best_val = 0, efficiency = 0), 
                  approximate_design = data.frame(support_points=approx_design$par, weight = rep(0.25, 4)))
  hb_list$approximate_design <- hb_list$approximate_design %>% arrange(support_points)
  
  # Criterion Value of each replication.
  hb_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val) 
  # Design points found in each replication.
  hb_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par) 
  
  best_idx <- which.max(hb_list$val)
  # Best criterion value among all replications.
  hb_list$result$best_val <- hb_list$val[best_idx]
  # Design points correspond to the best criterion value.
  hb_list$result$design_points <- hb_list$design_points[,best_idx]
  # Efficiency correspond to the best criterion value.
  hb_list$result$efficiency <- (hb_list$val[best_idx] / approx_design$val) ^ (1/nPoints)
  
  hb_list
}

## Hunt-Bowman tau-optimal

hb_tauoptimal <- function(d, loc){
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  n <- length(d)
  pen = 999999999
  
  mat_list <- lapply(d, function(x) 1/n * hb_mat(x, c1, tau, b0, b1))
  M <- Reduce("+", mat_list)
  diag(M) <- diag(M) + 1e-10
  b <- matrix(c(0, 1, 0, 0))
  
  if (rcond(M) < 2.220446e-16){
    res = pen
  } 
  else{
    M_inv <- solve(M)
    res <- t(b) %*% M_inv %*% b
  }
  
  res
}

# Find tau-optimal exact design for the Hunt-Bowman model
hb_tauoptimal_pso <- function(nPoints, parms, psoinfo){
  lb <- rep(0, nPoints)
  ub <- rep(0.15, nPoints)
  
  pso_res <- globpso(objFunc = hb_tauoptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$par <- pso_res$par |> round(4) |> sort()
  pso_res
}

# Find tau-optimal approximate design for the Hunt-Bowman model
hb_tauoptimal_approx <- function(parms){
  psoinfo <- psoinfo_setting()
  approx_design <- hb_tauoptimal_pso(nPoints = 2, parms, psoinfo)
  approx_design
}

# Replicate m PSO results of the tau-optimal exact design for Hunt-Bowman model
hb_tauoptimal_pso_rep <- function(nRep, nPoints = 2, parms, psoinfo){
  
  # The tau-optimal approximate design for the Hunt-Bowman model under certain parameter set.
  approx_design <- hb_tauoptimal_approx(parms)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- hb_tauoptimal_pso(nPoints, parms, psoinfo))
  
  hb_list <- list(val = c(), design_points = c(), history = c(), 
                  result = list(best_val = 0), 
                  approximate_design = data.frame(support_points=approx_design$par, weight = rep(0.5, 2)))
  hb_list$approximate_design <- hb_list$approximate_design %>% arrange(support_points)
  
  # Criterion value of each replication.
  hb_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val) 
  # Design points found in each replication.
  hb_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par) 
  
  best_idx <- which.min(hb_list$val)
  # Best criterion value among all replications.
  hb_list$result$best_val <- hb_list$val[best_idx]
  # Design points correspond to the best criterion value.
  hb_list$result$design_points <- hb_list$design_points[,best_idx]
  # Efficiency correspond to the best criterion value.
  hb_list$result$efficiency <- (approx_design$val / hb_list$val[best_idx])^(1/nPoints)
  
  hb_list
}

## Hunt-Bowman h-optimal

hb_hoptimal <- function(d, loc){
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  n <- length(d)
  pen = 999999999
  
  mat_list <- lapply(d, function(x) 1/n * hb_mat(x, c1, tau, b0, b1))
  M <- Reduce("+", mat_list)
  h <- matrix(c(-tau, -c1, 0, 0))
  
  if (rcond(M) < 2.220446e-16){
    res = pen
  } 
  else{
    M_inv <- solve(M)
    res <- t(h) %*% M_inv %*% h
  }
  
  res
}

hb_hoptimal_approx <- function(input, loc){
  l <- length(input)
  n <- ceiling(l/2)
  d <- input[1:n]
  w <- input[(n+1):l]
  w[n] <- 1 - sum(w)
  
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  pen = 999999999
  
  mat_list <- lapply(1:n, function(x) w[x] * hb_mat(d[x], c1, tau, b0, b1))
  M <- Reduce("+", mat_list)
  h <- matrix(c(-tau, -c1, 0, 0))
  
  if (rcond(M) < 2.220446e-16 || w[n] < 0){
    res = pen
  } 
  else{
    M_inv <- solve(M)
    res <- t(h) %*% M_inv %*% h
  }
  
  res
}

# Find h-optimal exact design for the Hunt-Bowman model
hb_hoptimal_pso <- function(nPoints, parms, psoinfo){
  lb <- rep(0, nPoints)
  ub <- rep(0.15, nPoints)
  
  pso_res <- globpso(objFunc = hb_hoptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$par <- pso_res$par |> round(4) |> sort()
  pso_res
}

# Find h-optimal approximate design for the Hunt-Bowman model
hb_hoptimal_approx_pso <- function(nPoints, parms){
  psoinfo <- psoinfo_setting(nSwarms = 256, Iters = 3000)
  lb <- rep(0, 2 * nPoints - 1)
  ub <- c(rep(0.15, nPoints), rep(1, nPoints - 1))
  
  pso_res <- globpso(objFunc = hb_hoptimal_approx, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$par <- c(pso_res$par, 1 - sum(pso_res$par[(nPoints+1):(2*nPoints-1)]))
  pso_res$par <- pso_res$par |> round(4)
  pso_res
}

# Replicate m PSO results of the h-optimal exact design for Hunt-Bowman model
hb_hoptimal_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo){
  
  # The h-optimal approximate design for the Hunt-Bowman model under certain parameter set.
  approx_design <- hb_hoptimal_approx_pso(nPoints = 4, parms = parms)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- hb_hoptimal_pso(nPoints, parms, psoinfo))
  
  hb_list <- list(val = c(), design_points = c(), history = c(), 
                  result = list(best_val = 0), 
                  approximate_design = data.frame(support_points = approx_design$par[1:4], 
                                                  weight = approx_design$par[5:8]))
  hb_list$approximate_design <- hb_list$approximate_design %>% arrange(support_points)
  idx0 <- print(match(0, hb_list$approximate_design$weight, -1))
  if (idx0 != -1) hb_list$approximate_design <- hb_list$approximate_design[-idx0, ]
  
  # Criterion value of each replication.
  hb_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val) 
  # Design points found in each replication.
  hb_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par) 
  
  best_idx <- which.min(hb_list$val)
  # Best criterion value among all replications.
  hb_list$result$best_val <- hb_list$val[best_idx]
  # Design points correspond to the best criterion value.
  hb_list$result$design_points <- hb_list$design_points[,best_idx]
  # Efficiency correspond to the best criterion value.
  hb_list$result$efficiency <- (approx_design$val / hb_list$val[best_idx])^(1/nPoints)
  
  hb_list
}


### exp-log model

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

## exp-log D-optimal
exp_log_doptimal <- function(d, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  
  M <- exp_log_mat(d, c0, c1, b0, b1)
  
  -det(M)
}

# Find D-optimal exact design for the exp-log model
exp_log_doptimal_pso <- function(nPoints, parms, psoinfo){
  lb <- rep(0, nPoints)
  ub <- rep(0.15, nPoints)
  
  pso_res <- globpso(objFunc = exp_log_doptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$val <- pso_res$val * -1
  pso_res$par <- pso_res$par |> round(4) |> sort()
  pso_res
}

# Find D-optimal approximate design for the exp-log model
exp_log_doptimal_approx <- function(parms){
  psoinfo <- psoinfo_setting()
  approx_design <- exp_log_doptimal_pso(nPoints = 4, parms, psoinfo)
  approx_design
}

# Replicate m PSO results of the D-optimal exact design for exp-log model
exp_log_doptimal_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo){
  
  # The D-optimal approximate design for the exp-log model under certain parameter set.
  approx_design <- exp_log_doptimal_approx(parms)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- exp_log_doptimal_pso(nPoints, parms, psoinfo))
  
  exp_log_list <- list(val = c(), design_points = c(), history = c(), 
                       result = list(best_val = 0), 
                       approximate_design = data.frame(support_points=approx_design$par, weight = rep(0.25, 4)))
  exp_log_list$approximate_design <- exp_log_list$approximate_design %>% arrange(support_points)
  
  # Criterion value of each replication.
  exp_log_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  # Design points found in each replication.
  exp_log_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par)
  
  best_idx <- which.max(exp_log_list$val)
  # Best criterion value among all replications.
  exp_log_list$result$best_val <- exp_log_list$val[best_idx]
  # Design points correspond to the best criterion value.
  exp_log_list$result$design_points <- exp_log_list$design_points[,best_idx]
  # Efficiency correspond to the best criterion value.
  exp_log_list$result$efficiency <- (exp_log_list$val[best_idx] / approx_design$val)^(1/nPoints)
  
  exp_log_list
}


## exp-log model h-optimal
exp_log_hoptimal <- function(d, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  pen = 999999999 # Penalty value
  
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
  
  res
}

exp_log_hoptimal_approx <- function(input, loc){
  l <- length(input)
  n <- ceiling(l/2)
  d <- input[1:n]
  w <- input[(n+1):l]
  w[n] <- 1 - sum(w)
  
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  pen = 999999999
  
  h1 <- -c1
  h2 <- -c0
  h3 <- b1 * exp(b0) * (1 - exp(b0)) / (exp(b0) + 1)^3
  h4 <- exp(b0) / (exp(b0) + 1)^2 
  h <- matrix(c(h1, h2, h3, h4))
  
  mat_list <- lapply(1:n, function(x) w[x] * exp_log_f(d[x], c0, c1, b0, b1))
  M <- Reduce("+", mat_list)
  
  if (rcond(M) < 2.220446e-16 || w[n] < 0){
    res = pen
  } 
  else{
    M_inv <- solve(M)
    res <- t(h) %*% M_inv %*% h
  }
  
  res
}

# Find h-optimal exact design for the exp-log model
exp_log_hoptimal_pso <- function(nPoints, parms, psoinfo){
  
  lb <- rep(0, nPoints)
  ub <- rep(0.15, nPoints)
  
  pso_res <- globpso(objFunc = exp_log_hoptimal, 
                     lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, 
                     loc = parms)
  
  pso_res$val <- pso_res$val
  pso_res$par <- c(pso_res$par) |> sort() |> round(4)
  pso_res$history <- pso_res$history
  pso_res
}

# Find h-optimal approximate design for the exp-log model
exp_log_hoptimal_approx_pso <- function(nPoints, parms){
  psoinfo <- psoinfo_setting(nSwarms = 256, Iters = 2000)
  lb <- rep(0, 2 * nPoints - 1)
  ub <- c(rep(0.15, nPoints), rep(1, nPoints - 1))
  
  pso_res <- globpso(objFunc = exp_log_hoptimal_approx, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$par <- c(pso_res$par, 1 - sum(pso_res$par[(nPoints+1):(2*nPoints-1)]))
  pso_res$par <- pso_res$par |> round(4)
  pso_res
}

# Replicate m PSO results of the h-optimal exact design for exp-log model
exp_log_hoptimal_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo){
  
  # The h-optimal approximate design for the exp-log model under certain parameter set.
  approx_design <- exp_log_hoptimal_approx_pso(4, parms)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- exp_log_hoptimal_pso(nPoints, parms, psoinfo))
  
  exp_log_list <- list(val = c(), design_points = c(), history = c(), 
                       result = list(best_val = 0, design_points = c()), 
                       approximate_design = data.frame(support_points = approx_design$par[1:4], 
                                                       weight = approx_design$par[5:8]))
  exp_log_list$approximate_design <- exp_log_list$approximate_design %>% arrange(support_points)
  print(exp_log_list$approximate_design)
  
  # Criterion value of each replication.
  exp_log_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  # Design points found in each replication.
  exp_log_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par)
  
  best_idx <- which.min(exp_log_list$val)
  # Best criterion value among all replications.
  exp_log_list$result$best_val <- exp_log_list$val[best_idx]
  # Design points correspond to the best criterion value.
  exp_log_list$result$design_points <- exp_log_list$design_points[,best_idx]
  # Efficiency correspond to the best criterion value.
  exp_log_list$result$efficiency <- (approx_design$val / exp_log_list$val[best_idx])^(1/nPoints)
  
  exp_log_list
}


## exp-log model tau-optimal 

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
    res = t(b) %*% M_inv %*% b
  }
  
  res
}

# Find tau-optimal exact design for the exp-log model
exp_log_tauoptimal_pso <- function(nPoints, parms, psoinfo){
  
  lb <- c(rep(0, nPoints))
  ub <- c(rep(0.15, nPoints))
  # Evaluate the tau value
  tau <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, 
                 c0 = parms[1], c1 = parms[2], b0 = parms[3], b1 = parms[4])$root
  
  pso_res <- globpso(objFunc = exp_log_tauoptimal, 
                     lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, 
                     loc = c(parms, tau))
  
  pso_res$val <- pso_res$val
  pso_res$par <- pso_res$par |> sort() |> round(4)
  pso_res
}

# Find tau-optimal approximate design for the exp-log model
exp_log_tauoptimal_approx <- function(parms){
  psoinfo <- psoinfo_setting()
  approx_design <- exp_log_tauoptimal_pso(nPoints = 2, parms, psoinfo)
  approx_design
}

# Replicate m PSO results of the tau-optimal exact design for exp-log model
exp_log_tauoptimal_pso_rep <- function(nRep, nPoints = 2, parms, psoinfo){
  
  # The tau-optimal approximate design for the exp-log model under certain parameter set.
  approx_design <- exp_log_tauoptimal_approx(parms)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- exp_log_tauoptimal_pso(nPoints, parms, psoinfo))
  
  exp_log_list <- list(val = c(), design_points = c(), history = c(), 
                       result = list(best_val = 0, design_points = c()), 
                       approximate_design = data.frame(support_points=approx_design$par, weight = rep(0.5, 2)))
  
  # Criterion value of each replication.
  exp_log_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  # Design points found in each replication.
  exp_log_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par)
  
  best_idx <- which.min(exp_log_list$val)
  # Best criterion value among all replications.
  exp_log_list$result$best_val <- exp_log_list$val[best_idx]
  # Design points correspond to the best criterion value.
  exp_log_list$result$design_points <- exp_log_list$design_points[,best_idx]
  # Efficiency correspond to the best criterion value.
  exp_log_list$result$efficiency <- (approx_design$val / exp_log_list$val[best_idx])^(1/nPoints)
  
  exp_log_list
}

## Logistic models

## Simple logistic model D-optimal
logistic_params <- function(alpha, beta){
  c(alpha, beta)
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

# Find D-optimal exact design for the simple logistic model
logistic_doptimal_pso <- function(nPoints, parms, psoinfo){
  lb <- rep(-10, nPoints)
  ub <- rep(10, nPoints)
  
  pso_res <- globpso(objFunc = logistic_doptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$val <- -1 * pso_res$val
  pso_res$par <- pso_res$par |> round(4) |> sort()
  pso_res
}

# Find D-optimal approximate design for the simple logistic model
logistic_doptimal_approx <- function(parms){
  psoinfo <- psoinfo_setting()
  approx_design <- logistic_doptimal_pso(nPoints = 2, parms, psoinfo)
  approx_design
}

# Replicate m PSO results of the D-optimal exact design for simple logistic model
logistic_pso_rep <- function(nRep, nPoints = 2, parms, psoinfo){
  
  # The D-optimal approximate design for the simple logistic model under certain parameter set.
  approx_design <- logistic_doptimal_approx(parms)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- logistic_d_pso(nPoints, parms, psoinfo))
  
  logistic_list <- list(design_points = c(), val = c(), 
                        result = list(best_val = 0, design_points = c(), efficiency = 0), 
                        approximate_design = data.frame(support_points=approx_design$par, weight = rep(0.5, 2)))
  
  # Criterion value of each replication.
  logistic_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  # Design points found in each replication.
  logistic_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par)
  
  best_idx <- which.max(logistic_list$val)
  # Best criterion value among all replications.
  logistic_list$result$best_val <- logistic_list$val[best_idx]
  # Design points correspond to the best criterion value.
  logistic_list$result$design_points <- logistic_list$design_points[, best_idx] |> sort() |> round(4)
  # Efficiency correspond to the best criterion value.
  logistic_list$result$efficiency <- (logistic_list$val[best_idx] / approx_design$val)^(1/nPoints)
  
  logistic_list
}

## Quadratic logistic D-optimal
qlogistic_params <- function(alpha, beta1, beta2){
  c(alpha, beta1, beta2)
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

qlogistic_doptimal_approx <- function(input, loc){
  d <- input[1:4]
  w <- input[5:7]
  w[4] <- 1 - sum(w)
  
  a <- loc[1]
  b1 <- loc[2]
  b2 <- loc[3]
  n <- length(d)
  pen = 999999999
  
  theta <- a + b1 * d + b2 * d^2
  omega <- exp(theta) / (1 + exp(theta))^2
  
  if (w[4] < 0){
    res = pen
  } else{
    inf_mat <- matrix(c(sum(w * omega), sum(w * d * omega), sum(w * d^2 * omega), 
                        sum(w * d * omega), sum(w * d^2 * omega), sum(w * d^3 * omega), 
                        sum(w * d^2 * omega), sum(w * d^3 * omega), sum(w * d^4 * omega)), 
                      nrow = 3)
    res = -det(inf_mat)
  }
  
  res
}

# Find D-optimal exact design for the quadratic logistic model
qlogistic_pso <- function(nPoints, parms, psoinfo){
  lb <- rep(-10, nPoints)
  ub <- rep(10, nPoints)
  
  pso_res <- globpso(objFunc = qlogistic_doptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$val <- -1 * pso_res$val
  pso_res$par <- pso_res$par |> round(4) |> sort()
  pso_res
}

# Find D-optimal approximate design for the quadratic logistic model
qlogistic_approx_pso <- function(parms){
  lb <- c(rep(-10, 4), rep(0, 3))
  ub <- c(rep(10, 4), rep(1, 3))
  psoinfo <- psoinfo_setting()
  
  pso_res <- globpso(objFunc = qlogistic_doptimal_approx, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$val <- -1 * pso_res$val
  pso_res$par <- c(pso_res$par, 1 - sum(pso_res$par[5:7])) |> round(4)
  pso_res
}

# Replicate m PSO results of the D-optimal exact design for quadratic logistic model
qlogistic_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo){
  
  # The D-optimal approximate design for the quadratic logistic model under certain parameter set.
  approx_design <- qlogistic_approx_pso(parms)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- qlogistic_pso(nPoints, parms, psoinfo))
  
  qlogistic_list <- list(design_points = c(), val = c())
  # Criterion value of each replication.
  qlogistic_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  # Design points found in each replication.
  qlogistic_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par)
  
  
  qlogistic_list$approximate_design <- data.frame(support_points = approx_design$par[1:4], 
                                                  weight = approx_design$par[5:8])
  qlogistic_list$approximate_design <- qlogistic_list$approximate_design %>% arrange(support_points)
  idx0 <- match(0, qlogistic_list$approximate_design$weight, -1)
  if (idx0 != -1) qlogistic_list$approximate_design <- qlogistic_list$approximate_design[-idx0, ]
  
  best_idx <- which.max(qlogistic_list$val)
  # Best criterion value among all replications.
  qlogistic_list$result$best_val <- qlogistic_list$val[best_idx]
  # Design points correspond to the best criterion value.
  qlogistic_list$result$design_points <- qlogistic_list$design_points[, best_idx] |> sort() |> round(4)
  # Efficiency correspond to the best criterion value.
  qlogistic_list$result$efficiency <- (qlogistic_list$val[best_idx] / approx_design$val)^(1/nPoints)
  
  qlogistic_list
}


## Cubic logistic models
clogistic_params <- function(alpha, beta1, beta2, beta3){
  c(alpha, beta1, beta2, beta3)
}

clogistic_doptimal <- function(d, loc){
  a <- loc[1]
  b1 <- loc[2]
  b2 <- loc[3]
  b3 <- loc[4]
  n <- length(d)
  
  theta <- a + b1 * d + b2 * d^2 + b3 * d^3
  omega <- exp(theta) / (1 + exp(theta))^2
  
  inf_mat <- matrix(c(sum(omega), sum(d * omega), sum(d^2 * omega), sum(d^3 * omega), 
                      sum(d * omega), sum(d^2 * omega), sum(d^3 * omega),  sum(d^4 * omega), 
                      sum(d^2 * omega), sum(d^3 * omega), sum(d^4 * omega), sum(d^5 * omega), 
                      sum(d^3 * omega), sum(d^4 * omega), sum(d^5 * omega), sum(d^6 * omega)), 
                    nrow = 4)
  
  -det(inf_mat) / n^4
}

clogistic_doptimal_approx <- function(input, loc){
  a <- loc[1]
  b1 <- loc[2]
  b2 <- loc[3]
  b3 <- loc[4]
  
  d <- input[1:5]
  weight <- input[6:9]
  weight[5] <- 1- sum(weight)
  
  
  theta <- a + b1 * d + b2 * d^2 + b3 * d^3
  omega <- exp(theta) / (1 + exp(theta))^2
  
  
  inf_mat <- matrix(c(sum(weight * omega), sum(weight * d * omega), 
                      sum(weight * d^2 * omega), sum(weight * d^3 * omega), 
                      sum(weight * d * omega), sum(weight * d^2 * omega), 
                      sum(weight * d^3 * omega),  sum(weight * d^4 * omega), 
                      sum(weight * d^2 * omega), sum(weight * d^3 * omega), 
                      sum(weight * d^4 * omega), sum(weight * d^5 * omega), 
                      sum(weight * d^3 * omega), sum(weight * d^4 * omega), 
                      sum(weight * d^5 * omega), sum(weight * d^6 * omega)), 
                    nrow = 4)
  
  if (weight[5] < 0){
    res = 9999999
  } else {
    res = -det(inf_mat)
  }
  
  res
}

# Find D-optimal exact design for the cubic logistic model
clogistic_pso <- function(npoints, parms, psoinfo){
  lb <- rep(-5, npoints)
  ub <- rep(5, npoints)
  
  pso_res <- globpso(objFunc = clogistic_doptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$val <- pso_res$val * -1
  pso_res$par <- pso_res$par |> round(4) |> sort()
  pso_res
}

# Find D-optimal approximate design for the cubic logistic model
clogistic_approx_pso <- function(parms){
  lb <- c(rep(-5, 5), rep(0, 4))
  ub <- c(rep(5, 5), rep(1, 4))
  psoinfo <- psoinfo_setting()
  
  pso_res <- globpso(objFunc = clogistic_doptimal_approx, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$val <- pso_res$val * -1
  pso_res$par <- c(pso_res$par, 1 - sum(pso_res$par[6:9])) |> round(4)
  pso_res
}

# Replicate m PSO results of the D-optimal exact design for cubic logistic model
clogistic_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo){
  
  # The D-optimal approximate design for the cubic logistic model under certain parameter set.
  approx_design <- clogistic_approx_pso(parms)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- clogistic_pso(nPoints, parms, psoinfo))
  
  clogistic_list <- list(design_points = c(), val = c())
  # Criterion value of each replication.
  clogistic_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  # Design points found in each replication.
  clogistic_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par)
  
  clogistic_list$approximate_design <- data.frame(support_points = approx_design$par[1:5], 
                                                  weight = approx_design$par[6:10])
  clogistic_list$approximate_design <- clogistic_list$approximate_design %>% arrange(support_points)
  idx0 <- match(0, clogistic_list$approximate_design$weight, -1)
  if (idx0 != -1) clogistic_list$approximate_design <- clogistic_list$approximate_design[-idx0, ]
  
  best_idx <- which.max(clogistic_list$val)
  # Best criterion value among all replications.
  clogistic_list$result$best_val <- clogistic_list$val[best_idx]
  # Design points correspond to the best criterion value.
  clogistic_list$result$design_points <- clogistic_list$design_points[, best_idx] |> sort() |> round(4)
  # Efficiency correspond to the best criterion value.
  clogistic_list$result$efficiency <- (clogistic_list$val[best_idx] / approx_design$val)^(1/nPoints)
  
  clogistic_list
}

# Run the application 
shinyApp(ui = ui, server = server)
