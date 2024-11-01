library(shiny)
library(shinyvalidate)
library(globpso)
library(waiter)
library(dplyr)
library(ggplot2)
library(AlgDesign)
library(MASS)
library(shinyjs)

# Define UI for application that draws a histogram
ui <- fluidPage(
  shinyjs::useShinyjs(),
  withMathJax(), 
  tabsetPanel(
    selected = "User Manual", 
    type = "tabs", 
    id = "mainpanel", 
    
    tabPanel(
      "User Manual",
      
      tags$h2("Usage"), 
      tags$p("This app is for finding optimal designs for studying hormesis in toxicology experiments. 
      The key tool is a state-of-the-art nature-inspired metaheuristic algorithm called the particle swarm 
      optimization(PSO) briefly described in the paper \"Optimal Experimental Designs for Big and Small 
      Experiments in Toxicology with Applications to Studying Hormesis via Metaheuristics\".
      This app finds various kinds of optimal designs for studying hormesis in large and 
      small experiments.  Since theory is available for finding large sample designs or approximate 
      designs, an emphasis here is for finding optimal exact designs for small experiments, which is becoming 
      common in toxicology studies. 
      We consider two statistical models for studying hormesis, namely the 
      Hunt-Bowman model and the Exp-Log model and find optimal designs for estimating model 
      parameters (D-optimality), detecting existence of hormesis (h-optimality) and estimating the threshold dose 
      (\\(\\tau\\)-optimality).  The algorithm requires the user to input the swarm size (i.e. number of particles 
      used for the search) and the number of iterations allowed for the search."),
      tags$p("Users are allowed to adjust the PSO parameters and change the model parameters. To start
      the search process, please note the following. 
     "), 
      
      tags$h2("Parameters"), 
      tags$div(
        tags$b("Lower / Upper Bound of Design Space:"),
        "The design space for the optimal design is the dose range of interest; the lower and upper bounds can be any numbers, except that 
        the upper bound must be a number larger than the lower bound. A larger design space, i.e. a wider dose range of interest, may 
        require more computation time.",
        tags$br(), tags$br(), 

        tags$b("N:"), 
        "Total number of the observations in the study.", 
        tags$br(), tags$br(),
        
        tags$b("Criterion:"), 
        "Users are allowed to choose one of three design criteria for studying hormesis: D, \\(\\tau\\) and h-optimality.  More details below.", 
        tags$br(), tags$br(),
        
        tags$b("Number of Particles:"), 
        "This is the swarm size; i.e. how many particles are used to search in the PSO algorithm. 
        The default size is 64 swarms.",
        tags$br(), tags$br(),
        
        tags$b("Number of Iterations:"), 
        "How many iterations for the PSO algorithim. The default is 1000 iterations.", 
        tags$br(), tags$br(),
        
        "Users are allowed to adjust the PSO parameters and change the default values of the model parameters 
        in the boxes provided. The latter are nominal values that represent our best guesses for the true 
        values of the parameters."
      ),
      
      tags$h2("Output"), 
      
      tags$h3("Optimal Design for Model 1: "), 
      tags$h4(strong("Model 1 is our model of interest. Supply the nominal value for its parameters first.")),
      
      tags$b("Approximate Design: "),
      "After running PSO, the app will show the optimal approximate design with the dose levels and
      their corresponding weights. The optimality of the D-optimal designs can be verified by an equivalence 
      theorem graphically, and the graph is shown on the rights side in the output page. Optimal approximate designs 
      are especially useful when the total sample size  is relatively large, say for example, N is at least 20.",
      tags$br(), tags$br(),
      
      tags$b("Efficient Rounding Method (ERM)-generated Design: "), 
      "The ERM-generated design is useful when N is relatively small. It is found from the optimal approximate design 
      using an efficient rounding rule [1] as follows:",
      tags$br(), 
      
      "Let \\(r_i\\) be the ceiling of \\(w_i (N - p/2)\\) for the number of replicates on each 
      support. If the sum of all \\(r_i\\)'s is smaller than N, the \\(r_i\\) with the lowest 
      \\(r_i / w_i\\) are added with one in sequence until all \\(r_i\\)'s sum up to one. And if 
      the sum of all \\(r_i\\)'s is greater than N, the \\(r_i\\) with the highest \\( (r_i-1) / w_i\\) 
      is subtracted with one until the sum of all \\(r_i\\)'s achives N. If there are multiple \\(r_i\\)'s
      having the same \\(r_i / w_i \\) or \\( (r_i - 1) / w_i \\), then the one with the largest \\(r_i\\) 
      will be subtracted or added with one first.",
      tags$br(), tags$br(),
      
      tags$b("Efficiency of the ERM-generated design relative to the approximate design:"),
      "The app displays the efficiency of the ERM-generated design.  This number measures the worth of a design and 
      is usually calculated from the ratio of its criterion value relative to the optimal approximate design 
      and is a number between 0 and 1.  If its efficiency value is 0.5, this means that the ERM-generated design 
      must be replicated twice for it to do as well as the optimal approximate design in terms of the criterion 
      value. For some criteria, like D-optimality, the ratio is suitably modified to retain the same interpretation.",
      tags$br(), tags$br(),
      
      tags$h3("Optimal Design for Model 2"),
      "Users can enter another set of nominal parameter values for model 2. These values typically represent second 
      best guesses or alternative values for the model parameters. PSO searches for optimal design under Model 2 and 
      evaluates its efficiency relative to that found for Model 1. The tab 'Comparison of Optimal Designs' enables such 
      comparison, including how the ERM-generated design compares to the approximate optimal designs for Models 1 and 2.",
      
      tags$h3("Equivalence Theorem Plot: "),
      "The equivalence theorem plot allows the user to check whether the optimal approximate design found 
      by PSO is truly optimal, that is it is optimal among all approximate designs for the problem at hand. 
      We visually check whether the graph attains the same peak values with a value of 0 at the dose levels of the  
      PSO-generated approximate optimal design. If it does, the design is optimal; otherwise, it is not. In tha app, 
      we provide such plots for confirming whether a design is D-optimal or not.",
      
      tags$h3("User Specified Design: "),
      
      "If the user has a design of interest, this app also provides 
      a function to evaluate its efficiency relative to the optimal approximate design and ERM-generated designs 
      for both model 1 and model 2.\
      For example, if  N = 8 and the design of interest has 2 observations at each of doses at 0, 0.1, 0.2 and 0.3, 
      the app computes its efficiencies relative to the optimal approximate design and ERM-generated design.",
      tags$br(), tags$br(),
      
      tags$b("Doses: "), 
      "Input the doses of the user-specified design with each dose separated by a comma. For example, for the above 
      user-specified design, we input '0,0.1,0.2,0.3' and make sure the number of different doses should not exceed 
      the total number of observations N, and the doses should all be inside the design space specified in the 
      setting of optimal design.",
      tags$br(), tags$br(),
      
      tags$b("Number of Replicates: "),
      "The number of replicates at each dose of the user-specified design. Each number should also be
      separated by a comma, e.g. '(2,2,2,2)' and the sum of the number of replicates must sum 
      up to  N,. Each number should  correspond to one dose, 
      that is, the dose levels should  align with the number of replicates.",
      tags$br(), tags$br(),
      
      tags$b("Efficiency: "), 
      "The output includes the efficiencies of the user-specified design relative to the optimal approximate 
      design and the ERM-generated design for both model 1 and model 2.",
      tags$br(), tags$br(),
      
      tags$h2("Particle Swarm Optimization"), 
      tags$p("Particle swarm optimization is a swarm-based metaheuristic algorithim that can
             optimize a target function without requiring technical assumptions. A swarm is a population 
             of randomly generated particles, which are candidates
             of the global optimum. The number of particles is the user-selected swarm
             size, and the \\(i^{th}\\) particle will search through the entire solution space iteratively.
             At iteration \\(t\\), the \\(i^{th}\\) particle gets to its current position \\(x^t_i\\) with velocity 
             \\(v^t_i\\) and moves to its next position \\(x^{t+1}_i\\) using velocity \\(v^{t+1}_i\\). 
             In each iteration, the particles update their velocities and positions using the following equations:"),
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
             and uniformly from \\([0,1]\\)"),
      tags$p("\\(\\omega\\): The inertia weight which represtns how active the particle is. 
             Can be a decreasing function or a fixed constant. Uniformly decreasing from 
             \\(0.9\\) to \\(0.4\\) in this app."),
      tags$p("\\(c_1\\) and \\(c_2\\): Two tuning parameters called cognitive coefficient 
             and social coefficient respectively. Can be viewed as how much information
             the particle used to update its velocity from the personal best and the global
             best. \\(c_1=c_2=2.05\\) in this app."),
      
      tags$h2("Hormesis"), 
      tags$p("Hormesis is a dose response relationship that has low-dose positive effect 
      and high-dose negative effect, and the threshold dose level is defined as the maximum 
      nonzero exposure level below which no adverse events above background response occur.
      Mathematically, it can be described by the following equation."),
      tags$p("$$
\\tau=\\tau(\\theta)=max \\{ d\\in\\Omega:\\mu(d,\\theta)\\leq\\mu(0,\\theta) \\}
$$"),
      tags$p("\\(\\mu(d,\\theta)\\) is the mean response function that characterize the 
      overall dose-response relationship, and \\(\\Omega = [0, \\hat{d}]\\) is the prespecified 
      dose interval."),
      
      
      tags$h2("D-optimal Exact Designs"), 
      tags$p("For estimating the model parameters in the mean function, a D-optimal design 
             maximizes the determinant of the information matrix, which is equivalent to minimizing
             the generalized variance of the estimated model parameters."),
      tags$p("$$
M(\\xi,\\theta) = \\sum_i f(d_i,\\theta)f^T(d_i,\\theta)\\omega_i
$$"),
      tags$p("$$
f(d,\\theta)=\\frac{\\partial}{\\partial\\theta}\\mu(d,\\theta)
$$"),
      tags$p("For exact designs, we  fix the weights \\(\\omega_i=\\frac{1}{N}\\) for 
             all design points. Here \\(N\\) is the numebr of design points."),
      
      
      tags$h2("\\(\\tau\\)-optimal Exact Designs"), 
      tags$p("For estimating the threshold dose level \\(\\tau\\), a \\(\\tau\\)-optimal
             design   minimizes the variance of the estimated threshold dose \\(\\tau\\), which can 
             be viewed as a special case of c-optimal criterion:"),
      tags$p("$$
b^T(\\theta)M^{-1}(\\xi,\\theta)b(\\theta)
$$"),
      tags$p("$$
b(\\theta)=\\frac{\\partial}{\\partial\\theta}\\tau(\\theta)
$$"),
      
      
      tags$h2("h-optimal Exact Designs"), 
      tags$p("For detecting the existence of hormesis, Dette et al. [3] proposed h-optimality
            , which can also be treated as a special case of c-optimal criterion:"),
      tags$p("$$
h^T(0,\\theta)M^-1(\\xi,\\theta)h(0,\\theta)
$$"),
      tags$p("$$
h(d,\\theta)=\\frac{\\partial f(d,\\theta)}{\\partial d}
$$"),
      
      
      tags$h2("Hunt-Bowman Model"), 
      tags$p("Hunt and Bowman [2] suggest a dose-response model whose mean function is characterized by a piecewise quadratic
             logistic model in the following way:"),
      tags$p("$$
 \\mu(d)=\\begin{cases} 
  c_1d^2+c_2d+\\kappa, & 0 \\leq d\\leq\\tau \\\\
  \\frac{1}{1+e^{\\beta_0-\\beta_1(d-\\tau)}}, & \\tau<d 
  \\end{cases}
 $$"),
      tags$p("Due to the constraints on the hormesis threshold, we have 
             \\(\\kappa=\\frac{1}{1+e^{\\beta_0}}\\) and 
             \\(c_2=-c_1 \\tau\\). Hence the parameter set for the Hunt-Bowman model is
             \\(\\theta=(c_1,\\tau,\\beta_0,\\beta_1)^T\\)."),
      
      
      tags$h2("Exp-Log Model"), 
      tags$p("Dette et al. [3] proposed a smooth analytic model that dose not involve the 
             threshold dose level parameter \\(\\tau\\).  It models the mean reponse as a sum of an exponential 
             decay and curve and a sigmoidal curve."),
      tags$p("$$
 \\mu(d) = c_0e^{-c_1d}+\\frac{1}{1+e^{\\beta_0-\\beta_1d}}
 $$"),
      tags$p("The parameter set of the Exp-Log model is 
             \\(\\theta=(c_0,c_1,\\beta_0,\\beta_1)^T\\)."),
      
      tags$h2("Reference"),
      tags$div(
        "[1] Pukelsheim, Friedrich, and Sabine Rieder. ", 
        tags$b("Efficient Rounding of Approximate Designs. "), 
        tags$i("Biometrika, "), 
        "vol. 79, no. 4, pp. 763–70, 1992. ", 
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
      )
      
    ),
    
    tabPanel(
      "Find Optimal Designs", 
      
      tabsetPanel(
        selected = "Optimal Designs for Hunt-Bowman Model", 
        type = "pills", 
        id = "HBD", 
        tabPanel("Optimal Designs for Hunt-Bowman Model", 
                 tags$h2("Optimal Designs for Hunt-Bowman Model"), 
                 tags$p("$$
 \\mu(d)=\\begin{cases} 
  c_1d^2+c_2d+\\kappa, & 0 \\leq d\\leq\\tau \\\\
  \\frac{1}{1+e^{\\beta_0-\\beta_1(d-\\tau)}}, & \\tau<d 
  \\end{cases}
 $$"),
                 
                 use_waiter(),
                 fluidRow(
                   column(width = 5, 
                          actionButton("hb_plot_response", "Plot dose-response relationship"), 
                          plotOutput("hb_plot"),
                          tags$h6(strong("CAUTION: The model parameters, including the design space must be meaningfully provided, 
 - a good way to check them is to examine the shapes of the resulting mean response curves for the models."))
                   ),
                   column(width = 3, 
                          column(width = 6,
                                 wellPanel(
                                   tags$h3("Model 1 Parameters"),
                                   numericInput("hb1_c1", 
                                                "\\(c_1\\)", 
                                                value = 170, 
                                                min = -Inf, 
                                                max = Inf),
                                   numericInput("hb1_tau", 
                                                "\\(\\tau\\)", 
                                                value = 0.04, 
                                                min = -Inf, 
                                                max = Inf),
                                   numericInput("hb1_b0", 
                                                "\\(\\beta_0\\)", 
                                                value = 1.46, 
                                                min = -Inf, 
                                                max = Inf),
                                   numericInput("hb1_b1", 
                                                "\\(\\beta_1\\)", 
                                                value = 40, 
                                                min = -Inf, 
                                                max = Inf),
                                 )
                          ), 
                          column(width = 6,
                                 wellPanel(
                                   tags$h3("Model 2 Parameters"),
                                   numericInput("hb2_c1", 
                                                "\\(c_1\\)", 
                                                value = 170, 
                                                min = -Inf, 
                                                max = Inf),
                                   numericInput("hb2_tau", 
                                                "\\(\\tau\\)", 
                                                value = 0.03, 
                                                min = -Inf, 
                                                max = Inf),
                                   numericInput("hb2_b0", 
                                                "\\(\\beta_0\\)", 
                                                value = 1.46, 
                                                min = -Inf, 
                                                max = Inf),
                                   numericInput("hb2_b1", 
                                                "\\(\\beta_1\\)", 
                                                value = 40, 
                                                min = -Inf, 
                                                max = Inf),
                                 )
                          )
                          ),
                   column(width = 4, 
                          column(width = 6, 
                                 wellPanel(
                                   tags$h3("Design Parameters"),
                                   numericInput("hb_lb", 
                                                "Lower Bound of Design Space", 
                                                value = 0, 
                                                min = -Inf, 
                                                max = Inf, 
                                                step = 0.01
                                   ),
                                   numericInput("hb_ub", 
                                                "Upper Bound of Design Space", 
                                                value = 0.15, 
                                                min = -Inf, 
                                                max = Inf, 
                                                step = 0.01
                                   ),
                                   numericInput("hb_dim", 
                                                "Total Number of Observations for the Study N", 
                                                value = 10, 
                                                min = 1, 
                                                step = 1
                                   ),
                                   selectInput("hb_criterion", 
                                               "Criterion", 
                                               c("D" = "D", 
                                                 "tau" = "tau", 
                                                 "h" = "h")
                                   ),
                                 )
                                 ), 
                          column(width = 6, 
                                 wellPanel(
                                   tags$h3("PSO Parameters"),
                                   numericInput("hb_approx_swarm", 
                                                "Number of Particles for PSO", 
                                                value = 64, 
                                                min = 1, 
                                                max = Inf, 
                                                step = 1),
                                   numericInput("hb_approx_iter", 
                                                "Number of Iterations for PSO", 
                                                value = 1000, 
                                                min = 1, 
                                                max = Inf, 
                                                step = 1
                                   ) 
                                 ),
                                 )
                          )
                   ),
                 
                 fluidRow(
                   column(width = 4, 
                          actionButton("hb1_pso", "Optimal Design for model 1"), 
                          tagAppendAttributes(verbatimTextOutput("hb1_pso_out"), style = "height:400px;"),
                          ), 
                   column(width = 4,
                          actionButton("hb2_pso", "Optimal Design for model 2"), 
                          tagAppendAttributes(verbatimTextOutput("hb2_pso_out"), style = "height:400px;"),
                   ),
                   column(width = 4, 
                          selectInput("hb_get_md", 
                                      "Equivalence Theorem Plot", 
                                      c("Model 1" = "md1", 
                                        "Model 2" = "md2")
                          ),
                          plotOutput("hb_get_plot", height = "365px"),
                   ), 
                 ),
                 
                 fluidRow(
                   column(width = 6,
                          actionButton("hb_comp", "Comparison of Optimal Designs"),
                          tagAppendAttributes(verbatimTextOutput("hb_comp_out"), style = "height:400px;"),
                   ),
                   column(width = 6,
                          actionButton("hb_spec", "User Specified Design"),
                          fluidRow(
                            column(width = 6, 
                                   textInput("hb_spec_sup", "Doses", "0,0.0375,0.075,0.1125,0.15"),
                            ),
                            column(width = 6, 
                                   textInput("hb_spec_N", "Number of Replicates", "2,2,2,2,2"),
                            ),
                          ),
                          tagAppendAttributes(verbatimTextOutput("hb_spec_out"), style = "height:325px;"),
                   ),
                 ),
        ), 
        
        tabPanel("Optiaml Designs for Exp-Log Model", 
                 tags$h2("Optimal Designs for Exp-Log Model"), 
                 tags$p("$$
 \\mu(d) = c_0e^{-c_1d}+\\frac{1}{1+e^{\\beta_0-\\beta_1d}}
 $$"), 
                 use_waiter(),
                 fluidRow(
                   column(width = 5, 
                          actionButton("el_plot_response", "Plot dose-response relationship"), 
                          plotOutput("el_plot"), 
                          tags$h6(strong("CAUTION: The model parameters, including the design space must be meaningfully provided, 
 - a good way to check them is to examine the shapes of the resulting mean response curves for the models."))
                   ),
                   column(width = 3, 
                          column(width = 6,
                                 wellPanel(
                                   tags$h3("Model 1 Parameters"),
                                   numericInput("el1_c0", 
                                                "\\(c_0\\)", 
                                                value = 0.15, 
                                                min = -Inf, 
                                                max = Inf), 
                                   numericInput("el1_c1", 
                                                "\\(c_1\\)", 
                                                value = 89, 
                                                min = -Inf, 
                                                max = Inf), 
                                   numericInput("el1_b0", 
                                                "\\(\\beta_0\\)", 
                                                value = 3.2, 
                                                min = -Inf, 
                                                max = Inf), 
                                   numericInput("el1_b1", 
                                                "\\(\\beta_1\\)", 
                                                value = 41, 
                                                min = -Inf, 
                                                max = Inf),
                                 )
                          ), 
                          column(width = 6,
                                 wellPanel(
                                   tags$h3("Model 2 Parameters"),
                                   numericInput("el2_c0", 
                                                "\\(c_0\\)", 
                                                value = 0.15, 
                                                min = -Inf, 
                                                max = Inf), 
                                   numericInput("el2_c1", 
                                                "\\(c_1\\)", 
                                                value = 70, 
                                                min = -Inf, 
                                                max = Inf), 
                                   numericInput("el2_b0", 
                                                "\\(\\beta_0\\)", 
                                                value = 3.2, 
                                                min = -Inf, 
                                                max = Inf), 
                                   numericInput("el2_b1", 
                                                "\\(\\beta_1\\)", 
                                                value = 41, 
                                                min = -Inf, 
                                                max = Inf),
                                 )
                          ), 
                   ),
                   column(width = 4, 
                          column(width = 6, 
                                 wellPanel(
                                   tags$h3("Design Parameters"),
                                   numericInput("el_lb", 
                                                "Lower Bound of Design Space", 
                                                value = 0, 
                                                min = -Inf, 
                                                max = Inf, 
                                                step = 0.01
                                   ),
                                   numericInput("el_ub", 
                                                "Upper Bound of Design Space", 
                                                value = 0.15, 
                                                min = -Inf, 
                                                max = Inf, 
                                                step = 0.01
                                   ),
                                   numericInput("el_dim", 
                                                "Total Number of Observations for the Study N", 
                                                value = 10, 
                                                min = 1, 
                                                step = 1
                                   ),
                                   selectInput("el_criterion", 
                                               "Criterion", 
                                               c("D" = "D", 
                                                 "tau" = "tau", 
                                                 "h" = "h")
                                   ),
                                 )
                          ), 
                          column(width = 6, 
                                 wellPanel(
                                   tags$h3("PSO Parameters"),
                                   numericInput("el_approx_swarm", 
                                                "Number of Particles for PSO", 
                                                value = 64, 
                                                min = 1, 
                                                max = Inf, 
                                                step = 1),
                                   numericInput("el_approx_iter", 
                                                "Number of Iterations for PSO", 
                                                value = 1000, 
                                                min = 1, 
                                                max = Inf, 
                                                step = 1
                                   ) 
                                 ),
                          )
                   )
                 ),
                 
                 fluidRow(
                   column(width = 4,
                          actionButton("el1_pso", "Optimal Design for model 1"), 
                          tagAppendAttributes(verbatimTextOutput("el1_pso_out"), style = "height:400px;"),
                   ),
                   column(width = 4,
                          actionButton("el2_pso", "Optimal Design for model 2"), 
                          tagAppendAttributes(verbatimTextOutput("el2_pso_out"), style = "height:400px;"),
                   ),
                   column(width = 4, 
                          selectInput("el_get_md", 
                                      "Equivalence Theorem Plot", 
                                      c("Model 1" = "md1", 
                                        "Model 2" = "md2")
                          ),
                          plotOutput("el_get_plot", height = "365px"),
                   ), 
                 ),
                 
                 fluidRow(
                   column(width = 6,
                          actionButton("el_comp", "Comparison of Optimal Designs"),
                          tagAppendAttributes(verbatimTextOutput("el_comp_out"), style = "height:400px;"),
                   ),
                   column(width = 6,
                          actionButton("el_spec", "User Specified Design"),
                          fluidRow(
                            column(width = 6, 
                                   textInput("el_spec_sup", "Doses", "0,0.0375,0.075,0.1125,0.15"),
                            ),
                            column(width = 6, 
                                   textInput("el_spec_N", "Number of Replicates", "2,2,2,2,2"),
                            ),
                          ),
                          tagAppendAttributes(verbatimTextOutput("el_spec_out"), style = "height:325px;"),
                   ),
                 ),
                 )
        )
        
      )
      
    ))


# Define server logic required to draw a histogram
server <- function(input, output) {
  values <- reactiveValues()
  values$hb1 <- list(approx = list(design = data.frame("Support" = c(0), "Weight" = c(1)), val = numeric()), 
                     eff_round = list(design = data.frame("Support" = c(0), "N" = c(0)), efficiency = numeric()), 
                     spec = list(design = data.frame("Support" = c(0), "N" = c(0)), 
                                 efficiency = list(approx = numeric(), erm = numeric())), 
                     drplot = ggplot(), getplot = ggplot(), 
                     comp = list(approx = numeric(), erm = numeric(), erm_approx = numeric()))
  
  values$hb2 <- list(approx = list(design = data.frame("Support" = c(0), "Weight" = c(1)), val = numeric(), comp_eff = numeric()), 
                     eff_round = list(design = data.frame("Support" = c(0), "N" = c(0)), 
                                      efficiency = numeric(), comp_eff = numeric(), comp_eff_erm = numeric()), 
                     spec = list(design = data.frame("Support" = c(0), "N" = c(0)), 
                                 efficiency = list(approx1 = numeric(), erm1 = numeric(), approx2 = numeric(), erm2 = numeric())), 
                     drplot = ggplot(), getplot = ggplot())
  
  values$el1 <- list(approx = list(design = data.frame("Support" = c(0), "Weight" = c(1)), val = numeric()), 
                     eff_round = list(design = data.frame("Support" = c(0), "N" = c(0)), efficiency = numeric()), 
                     spec = list(design = data.frame("Support" = c(0), "N" = c(0)), 
                                 efficiency = list(approx = numeric(), erm = numeric())), 
                     drplot = ggplot(), getplot = ggplot(), 
                     comp = list(approx = numeric(), erm = numeric(), erm_approx = numeric()))
  
  values$el2 <- list(approx = list(design = data.frame("Support" = c(0), "Weight" = c(1)), val = numeric(), comp_eff = numeric()), 
                     eff_round = list(design = data.frame("Support" = c(0), "N" = c(0)), 
                                      efficiency = numeric(), comp_eff = numeric(), comp_eff_erm = numeric()), 
                     spec = list(design = data.frame("Support" = c(0), "N" = c(0)), 
                                 efficiency = list(approx1 = numeric(), erm1 = numeric(), approx2 = numeric(), erm2 = numeric())), 
                     drplot = ggplot(), getplot = ggplot())
  
  compute <- reactiveValues(hb1_approx = FALSE, hb1_effround = FALSE, hb_spec = FALSE, hb_comp = FALSE,
                            hb2_approx = FALSE, hb2_effround = FALSE, hb2_spec = FALSE, 
                            el1_approx = FALSE, el1_effround = FALSE, el_spec = FALSE, el_comp = FALSE,
                            el2_approx = FALSE, el2_effround = FALSE, el2_spec = FALSE)
  
  spec_error <- reactiveValues(hb_support = FALSE, hb_N = FALSE, hb_align = FALSE, 
                               hb_dop = FALSE, hb_Nrep = FALSE, hb_space = FALSE,
                               el_support = FALSE, el_N = FALSE, el_align = FALSE, 
                               el_dop = FALSE, el_Nrep = FALSE, el_space = FALSE)
  
  observeEvent(input$hb_plot_response, {
    c11 = input$hb1_c1
    tau1 = input$hb1_tau
    b01 = input$hb1_b0
    b11 = input$hb1_b1
    
    c12 = input$hb2_c1
    tau2 = input$hb2_tau
    b02 = input$hb2_b0
    b12 = input$hb2_b1
    
    lb = input$hb_lb
    ub = input$hb_ub
    
    fp <- seq(lb, ub, by = (ub - lb)/100)
    hb_df <- data.frame(dose = rep(fp, 2), 
                        response = c(sapply(fp, function(x) hunt_bowman(x, c11, tau1, b01, b11)), sapply(fp, function(x) hunt_bowman(x, c12, tau2, b02, b12))), 
                        model = c(rep("model1", 101), rep("model2", 101)))
    values$hb1$drplot <- ggplot(data = hb_df, aes(x = dose, y = response, group = model)) + 
      geom_line(aes(color = model))
  })
  
  output$hb_plot <- renderPlot({
    drplot = values$hb1$drplot
    drplot
  })
  
  observeEvent(input$hb1_pso, {
    waiter1 <- waiter::Waiter$new(
      id = "hb1_pso_out", 
      html = spin_circle(), 
      color = "lightgrey"
    )
    waiter1$show()
    on.exit(waiter1$hide())
    
    c1 = input$hb1_c1
    tau = input$hb1_tau
    b0 = input$hb1_b0
    b1 = input$hb1_b1
    approx_nswarm = input$hb_approx_swarm
    approx_iter = input$hb_approx_iter
    criterion = input$hb_criterion
    lb = input$hb_lb
    ub = input$hb_ub
    ndp = input$hb_dim
    
    hb_par <- hb_parms(c1 = c1, tau = tau, b0 = b0, b1 = b1)
    psoinfo_approx <- psoinfo_setting(nSwarms = approx_nswarm, Iters = approx_iter)
    approx_res <- pso_approx(model = "HuntBowman", criterion = criterion, parms = hb_par, 
                             upper = ub, lower = lb, psoinfo_approx = psoinfo_approx)
    
    values$hb1$approx = approx_res
    
    approx <- approx_res$design |> round_approx()
    approx_d <- approx$Support
    approx_w <- approx$Weight
    hb_n <- length(approx_d)
    hb_grid <- seq(lb, ub, by = 0.0001)
    hb_mat_list <- lapply(1:hb_n, function(x) approx_w[x] * hb_mat(approx_d[x], c1, tau, b0, b1))
    hb_M <- Reduce("+", hb_mat_list)
    diag(hb_M) <- diag(hb_M) + 1e-10
    hb_M_inv <- solve(hb_M)
    
    if (criterion == "D"){
      hbd_get_val <- sapply(hb_grid, function(x) hbd_get(x, hb_par, hb_M_inv))
      get_df1 <- data.frame(dose = hb_grid, val = hbd_get_val)
      get_df2 <- data.frame(dose = approx_d, val = hbd_get_val[(approx_d / 0.0001)+1])
      values$hb1$getplot <- ggplot(data = get_df1, aes(x = dose, y = val)) +
        geom_line(stat = "identity") +
        geom_point(data = get_df2, aes(x = dose, y = val), shape = 1, size = 10)
    } else if (criterion == "tau"){
      text = paste("This function does not support \n tau-optimal designs.")
      values$hb1$getplot <- ggplot() + 
        annotate("text", x = 4, y = 25, size = 8, label = text) + 
        theme_void()
    } else if (criterion == "h"){
      text = paste("This function does not support \n h-optimal designs.")
      values$hb1$getplot <- ggplot() + 
        annotate("text", x = 4, y = 25, size = 8, label = text) + 
        theme_void()
    }
    
    compute$hb1_approx = TRUE 

    eff_round_design <- eff_round(approx = approx, model = "HuntBowman", criterion = criterion, 
                                  parms = hb_par, nPoints = ndp)
    values$hb1$eff_round = eff_round_design
    
    compute$hb1_effround = TRUE
  })
  
  observeEvent(input$hb2_pso, {
    waiter3 <- waiter::Waiter$new(
      id = "hb2_pso_out", 
      html = spin_circle(), 
      color = "lightgrey"
    )
    waiter3$show()
    on.exit(waiter3$hide())
    
    c1 = input$hb2_c1
    tau = input$hb2_tau
    b0 = input$hb2_b0
    b1 = input$hb2_b1
    approx_nswarm = input$hb_approx_swarm
    approx_iter = input$hb_approx_iter
    criterion = input$hb_criterion
    lb = input$hb_lb
    ub = input$hb_ub
    ndp = input$hb_dim
    
    hb_par <- hb_parms(c1 = c1, tau = tau, b0 = b0, b1 = b1)
    psoinfo_approx <- psoinfo_setting(nSwarms = approx_nswarm, Iters = approx_iter)
    approx_res <- pso_approx(model = "HuntBowman", criterion = criterion, parms = hb_par, 
                             upper = ub, lower = lb, psoinfo_approx = psoinfo_approx)
    
    
    c11 = input$hb1_c1
    tau1 = input$hb1_tau
    b01 = input$hb1_b0
    b11 = input$hb1_b1
    hb_par1 <- c(c11, tau1, b01, b11)
    approx_d <- approx_res$design$Support
    approx_w <- approx_res$design$Weight
    approx_w1 <- approx_w[-length(approx_w)]
    approx_val <- values$hb1$approx$val
    
    if (criterion == "D"){
      comp_val <- hb_doptimal_approx(c(approx_d, approx_w), loc = hb_par1)
      approx_res$comp_eff <- (comp_val / approx_val)^(1/4) |> as.numeric()
    }
    if (criterion == "tau"){
      comp_val <- hb_tauoptimal_approx(c(approx_d, approx_w), loc = hb_par1)
      approx_res$comp_eff <- (approx_val / comp_val) |> as.numeric()
    }
    if (criterion == "h"){
      comp_val <- hb_hoptimal_approx(c(approx_d, approx_w), loc = hb_par1)
      approx_res$comp_eff <- (approx_val / comp_val) |> as.numeric()
    }
    
    values$hb2$approx = approx_res
    
    approx <- approx_res$design |> round_approx()
    approx_d <- approx$Support
    approx_w <- approx$Weight
    hb_n <- length(approx_d)
    hb_grid <- seq(lb, ub, by = 0.0001)
    hb_mat_list <- lapply(1:hb_n, function(x) approx_w[x] * hb_mat(approx_d[x], c1, tau, b0, b1))
    hb_M <- Reduce("+", hb_mat_list)
    diag(hb_M) <- diag(hb_M) + 1e-10
    hb_M_inv <- solve(hb_M)
    
    if (criterion == "D"){
      hbd_get_val <- sapply(hb_grid, function(x) hbd_get(x, hb_par, hb_M_inv))
      get_df1 <- data.frame(dose = hb_grid, val = hbd_get_val)
      get_df2 <- data.frame(dose = approx_d, val = hbd_get_val[(approx_d / 0.0001)+1])
      values$hb2$getplot <- ggplot(data = get_df1, aes(x = dose, y = val)) +
        geom_line(stat = "identity") +
        geom_point(data = get_df2, aes(x = dose, y = val), shape = 1, size = 10)
    } else if (criterion == "tau"){
      text = paste("This function does not support tau-optimal designs.")
      values$hb2$getplot <- ggplot() + 
        annotate("text", x = 4, y = 25, size = 8, label = text) + 
        theme_void()
    } else if (criterion == "h"){
      text = paste("This function does not support h-optimal designs.")
      values$hb2$getplot <- ggplot() + 
        annotate("text", x = 4, y = 25, size = 8, label = text) + 
        theme_void()
    }
    
    compute$hb2_approx = TRUE
    
    eff_round_design <- eff_round(approx = approx, model = "HuntBowman", criterion = criterion, 
                                  parms = hb_par, nPoints = ndp)
    effround_d <- lapply(1:nrow(eff_round_design$design), 
                         function(x) rep(eff_round_design$design$Support[x], eff_round_design$design$N[x])) |> unlist()
    approx_val <- values$hb1$approx$val
    
    if (criterion == "D"){
      comp_val <- hb_doptimal(effround_d, loc = hb_par1)
      eff_round_design$comp_eff <- (comp_val / approx_val)^(1/4) |> as.numeric()
    }
    if (criterion == "tau"){
      comp_val <- hb_tauoptimal(effround_d, loc = hb_par1)
      eff_round_design$comp_eff <- (approx_val / comp_val) |> as.numeric()
    }
    if (criterion == "h"){
      comp_val <- hb_hoptimal(effround_d, loc = hb_par1)
      eff_round_design$comp_eff <- (approx_val / comp_val) |> as.numeric()
    }
    
    eff_round_design$comp_eff_erm <- eff_round_design$comp_eff / values$hb1$eff_round$efficiency
    values$hb2$eff_round = eff_round_design
    
    compute$hb2_effround = TRUE
  })
  
  observeEvent(input$hb_comp, {
    if(compute$hb1_approx && compute$hb2_approx){
      hb1_c1 = input$hb1_c1
      hb1_tau = input$hb1_tau
      hb1_b0 = input$hb1_b0
      hb1_b1 = input$hb1_b1
      hb2_c1 = input$hb2_c1
      hb2_tau = input$hb2_tau
      hb2_b0 = input$hb2_b0
      hb2_b1 = input$hb2_b1
      approx_nswarm = input$hb_approx_swarm
      approx_iter = input$hb_approx_iter
      criterion = input$hb_criterion
      lb = input$hb_lb
      ub = input$hb_ub
      ndp = input$hb_dim
      
      hb_par1 <- c(hb1_c1, hb1_tau, hb1_b0, hb1_b1)
      
      approx_d <- values$hb2$approx$design$Support
      approx_w <- values$hb2$approx$design$Weight[-length(values$hb2$approx$design$Weight)]
      hb1_approx_val <- values$hb1$approx$val
      hb1_erm_eff <- values$hb1$eff_round$efficiency
      hb2_erm_eff <- values$hb2$eff_round$efficiency
      if (criterion == "D"){
        hb2_approx_val <- hb_doptimal_approx(c(approx_d, approx_w), loc = hb_par1)
        approx_comp_eff <- (hb2_approx_val / hb1_approx_val)^(1/4) |> as.numeric()
        erm_approx_comp_eff <- hb2_erm_eff * approx_comp_eff |> as.numeric()
        erm_comp_eff <- erm_approx_comp_eff / hb1_erm_eff |> as.numeric()
      }
      if (criterion == "tau"){
        hb2_approx_val <- hb_tauoptimal_approx(c(approx_d, approx_w), loc = hb_par1)
        approx_comp_eff <- (hb1_approx_val / hb2_approx_val) |> as.numeric()
        erm_approx_comp_eff <- hb2_erm_eff * approx_comp_eff |> as.numeric()
        erm_comp_eff <- erm_approx_comp_eff / hb1_erm_eff |> as.numeric()
      }
      if (criterion == "h"){
        hb2_approx_val <- hb_hoptimal_approx(c(approx_d, approx_w), loc = hb_par1)
        approx_comp_eff <- (hb1_approx_val / hb2_approx_val) |> as.numeric()
        erm_approx_comp_eff <- hb2_erm_eff * approx_comp_eff |> as.numeric()
        erm_comp_eff <- erm_approx_comp_eff / hb1_erm_eff |> as.numeric()
      }
      
      values$hb1$comp$approx = approx_comp_eff
      values$hb1$comp$erm_approx = erm_approx_comp_eff
      values$hb1$comp$erm = erm_comp_eff
      
      compute$hb_comp = TRUE
    }
  })
  
  observeEvent(input$hb_spec, {
    if(compute$hb1_approx && compute$hb2_approx){
      hb1_c1 = input$hb1_c1
      hb1_tau = input$hb1_tau
      hb1_b0 = input$hb1_b0
      hb1_b1 = input$hb1_b1
      hb2_c1 = input$hb2_c1
      hb2_tau = input$hb2_tau
      hb2_b0 = input$hb2_b0
      hb2_b1 = input$hb2_b1
      approx_nswarm = input$hb_approx_swarm
      approx_iter = input$hb_approx_iter
      criterion = input$hb_criterion
      lb = input$hb_lb
      ub = input$hb_ub
      ndp = input$hb_dim
      
      hb_par1 <- c(hb1_c1, hb1_tau, hb1_b0, hb1_b1)
      hb_par2 <- c(hb2_c1, hb2_tau, hb2_b0, hb2_b1)
      
      sup_error <- FALSE
      N_error <- FALSE
      tryCatch({support <- input$hb_spec_sup |> strsplit(",") |> unlist() |> as.numeric()}, 
               warning = function(w) {sup_error <<- TRUE}, 
               error = function(e) {sup_error <<- TRUE})
      if (sup_error == FALSE) spec_error$hb_support = FALSE
      else if (sup_error == TRUE) spec_error$hb_support = TRUE
      
      tryCatch({N <- input$hb_spec_N |> strsplit(",") |> unlist() |> as.numeric()}, 
               warning = function(w) {N_error <<- TRUE}, 
               error = function(e) {N_error <<- TRUE})
      if (N_error == FALSE) spec_error$hb_N = FALSE
      else if (N_error == TRUE) spec_error$hb_N = TRUE
      
      if(spec_error$hb_support == FALSE & spec_error$hb_N == FALSE){
        if (length(support) != length(N)){
          spec_error$hb_align <- TRUE
        } else {
          spec_error$hb_align <- FALSE
        }
        
        if (criterion == "D" & length(support) < 4){
          spec_error$hb_dop <- TRUE
        } else {
          spec_error$hb_dop <- FALSE
        }
        
        if(sum(N) != ndp){
          spec_error$hb_Nrep <- TRUE
        } else {
          spec_error$hb_Nrep <- FALSE
        }
        
        if(any(support < lb) || any(support > ub)){
          spec_error$hb_space = TRUE
        }
      } 
      if(!(TRUE %in% c(spec_error$hb_support, spec_error$hb_N, spec_error$hb_align))){
        hb_dose <- lapply(1:length(N), function(x) rep(support[x], N[x])) |> unlist()
        
        if(criterion == "D"){
          spec_val1 <- hb_doptimal(hb_dose, hb_par1)
          spec_val2 <- hb_doptimal(hb_dose, hb_par2)
          approx1 <- values$hb1$approx$val
          approx2 <- values$hb2$approx$val
          eff1 <- (spec_val1 / approx1)^(1/4)
          eff2 <- (spec_val2 / approx2)^(1/4)
        } else if (criterion == "tau"){
          spec_val1 <- hb_tauoptimal(hb_dose, hb_par1)
          spec_val2 <- hb_tauoptimal(hb_dose, hb_par2)
          approx1 <- values$hb1$approx$val
          approx2 <- values$hb2$approx$val
          eff1 <- approx1 / spec_val1
          eff2 <- approx2 / spec_val2
        } else if (criterion == "h"){
          spec_val1 <- hb_hoptimal(hb_dose, hb_par1)
          spec_val2 <- hb_hoptimal(hb_dose, hb_par2)
          approx1 <- values$hb1$approx$val
          approx2 <- values$hb2$approx$val
          eff1 <- approx1 / spec_val1
          eff2 <- approx2 / spec_val2
        }
        
        values$hb2$spec$design <- data.frame(Support = support, N = N)
        values$hb2$spec$efficiency$approx1 <- eff1
        values$hb2$spec$efficiency$erm1 <- eff1 /  values$hb2$eff_round$efficiency
        values$hb2$spec$efficiency$approx2 <- eff2
        values$hb2$spec$efficiency$erm2 <- eff2 /  values$hb2$eff_round$efficiency
        values$hb2$spec$design <- data.frame(Support = support, N = N)
        
        compute$hb_spec = TRUE
      }
    }
  })
  
  output$hb1_pso_out <- renderPrint({
    if(compute$hb1_approx == FALSE){
      cat("Press the button above to compute!")
    } else {
      cat("Approximate Design: ", "\n")
      approx_df <- values$hb1$approx$design |> round_approx()
      print(approx_df, row.names = FALSE)
      cat("\nEfficient Rounding Method (ERM)-generated Design: ", "\n")
      erm_df <- values$hb1$eff_round$design |> round(4)
      print(erm_df, row.names = FALSE)
      cat("\nEfficiency of the ERM-generated design relative \nto the approximate design for Model 1: ", 
          values$hb1$eff_round$efficiency |> round(4))
      cat("\n\nThe value of the design efficiency represents its worth; 
a design with 0.5 efficiency means that it has to be 
replicated twice to do as well as the optimum design. 
Similar interpretation applies when the design is compared 
to a non-optimal design.")
    }
    
  })
  
  output$hb2_pso_out <- renderPrint({
    if(compute$hb2_approx == FALSE) {
      cat("Press the button above to compute!")
    } else {
      cat("Approximate Design: ", "\n")
      approx_df <- values$hb2$approx$design |> round_approx()
      print(approx_df, row.names = FALSE)
      cat("\n\nEfficient Rounding Method (ERM)-generated Design: ", "\n")
      erm_df <- values$hb2$eff_round$design |> round(4)
      print(erm_df, row.names = FALSE)
      cat("\nEfficiency of the ERM-generated design relative \nto the approximate design for Model 2: ", 
          values$hb1$eff_round$efficiency |> round(4))
    }
    
  })
  
  output$hb_spec_out <- renderPrint({
    if(compute$hb1_approx == FALSE){
      cat("Please compute the optimal design for model 1 first!")
    } else if (compute$hb2_approx == FALSE){
      cat("Please compute the optimal design for model 2 first!")
    } else if (compute$hb_spec == FALSE){
      cat("Press the button above to compute!")
      cat("\nEach number must be seperated by a single\ncomma. e.g. 0,0.0375,0.075,0.1125,0.15")
    } else if (spec_error$hb_support){
      cat("Error in the format of Support.\nEach number must be seperated by a single\ncomma. e.g. 2,2,2,2,2")
    } else if (spec_error$hb_N){
      cat("Error in the format of N.\nEach number must be seperated by a single comma.\ne.g. 2,2,3")
    } else if (spec_error$hb_align){
      cat("The number of support points must be\nthe same with the number of replicates!")
    } else {
      if (spec_error$hb_dop) cat("Warning:\nFor D-optimal, at least 4 support is needed\nto have efficiency > 0.\n\n")
      if (spec_error$hb_Nrep) cat("Warning:\nThe total number of replicates you entered\ndoes not equal to the optimal approximate\ndesign.\n\n")
      if (spec_error$hb_space) cat("Warning:\nSome dose levels are outside the design space, please check!\n\n")
      cat("User Specified Design:\n")
      print(values$hb2$spec$design, row.names = FALSE)
      cat("\nEfficiency of the user specified design relative \nto the approximate design for Model 1: ", 
          values$hb2$spec$efficiency$approx1 |> round(4))
      cat("\n\nEfficiency of the user specified design relative \nto the ERM-generated design for Model 1: ", 
          values$hb2$spec$efficiency$erm1 |> round(4))
      cat("\n\nEfficiency of the user specified design relative \nto the approximate design for Model 2: ", 
          values$hb2$spec$efficiency$approx2 |> round(4))
      cat("\n\nEfficiency of the user specified design relative \nto the ERM-generated design for Model 2: ", 
          values$hb2$spec$efficiency$erm2 |> round(4))
    }
    
  })
  
  output$hb_comp_out <- renderPrint({
    if(compute$hb1_approx == FALSE){
      cat("Please compute the optimal design for model 1 first!")
    } else if (compute$hb2_approx == FALSE){
      cat("Please compute the optimal design for model 2 first!")
    } else if (compute$hb_comp == FALSE){
      cat("Press the button above to compute!")
    } else {
      cat("Efficiency of the approximate design for Model 2 relative \nto the approximate design for Model 1: ", 
          values$hb1$comp$approx |> round(4))
      cat("\n\nEfficiency of the ERM-generated design for Model 2 relative \nto the approximate design for Model 1: ", 
          values$hb1$comp$erm_approx |> round(4))
      cat("\n\nEfficiency of the ERM-generated design for Model 2 relative \nto the ERM-generated design for Model 1: ", 
          values$hb1$comp$erm |> round(4))
    }
    
  })
  
  output$hb_get_plot <- renderPlot({
    md = input$hb_get_md
    if (md == "md1") getplot = values$hb1$getplot
    else if (md == "md2") getplot = values$hb2$getplot
    getplot
  })
  
  observeEvent(input$el_plot_response, {
    c01 = input$el1_c0
    c11 = input$el1_c1
    b01 = input$el1_b0
    b11 = input$el1_b1
    
    c02 = input$el2_c0
    c12 = input$el2_c1
    b02 = input$el2_b0
    b12 = input$el2_b1
    
    lb = input$el_lb
    ub = input$el_ub
    
    fp <- seq(lb, ub, by = (ub - lb)/100)
    el_df <- data.frame(dose = rep(fp, 2), 
                        response = c(sapply(fp, function(x) exp_log(x, c01, c11, b01, b11)), sapply(fp, function(x) exp_log(x, c02, c12, b02, b12))), 
                        model = c(rep("model1", 101), rep("model2", 101)))
    values$el1$drplot <- ggplot(data = el_df, aes(x = dose, y = response, group = model)) + 
      geom_line(aes(color = model))
  })
  
  output$el_plot <- renderPlot({
    drplot = values$el1$drplot
    drplot
  })
  
  observeEvent(input$el1_pso, {
    waiter5 <- waiter::Waiter$new(
      id = "el1_pso_out", 
      html = spin_circle(), 
      color = "lightgrey"
    )
    waiter5$show()
    on.exit(waiter5$hide())
    
    c0 = input$el1_c0
    c1 = input$el1_c1
    b0 = input$el1_b0
    b1 = input$el1_b1
    approx_nswarm = input$el_approx_swarm
    approx_iter = input$el_approx_iter
    criterion = input$el_criterion
    lb = input$el_lb
    ub = input$el_ub
    ndp = input$el_dim
    
    el_par <- c(c0, c1, b0, b1)
    psoinfo_approx <- psoinfo_setting(nSwarms = approx_nswarm, Iters = approx_iter)
    approx_res <- pso_approx(model = "ExpLog", criterion = criterion, parms = el_par, 
                             upper = ub, lower = lb, psoinfo_approx = psoinfo_approx)
    
    values$el1$approx = approx_res
    
    approx <- approx_res$design |> round_approx()
    approx_d <- approx$Support
    approx_w <- approx$Weight
    el_n <- length(approx_d)
    el_grid <- seq(lb, ub, by = 0.0001)
    el_mat_list <- lapply(1:el_n, function(x) approx_w[x] * el_mat(approx_d[x], c0, c1, b0, b1))
    el_M <- Reduce("+", el_mat_list)
    diag(el_M) <- diag(el_M) + 1e-10
    el_M_inv <- solve(el_M)
    
    if (criterion == "D"){
      eld_get_val <- sapply(el_grid, function(x) eld_get(x, el_par, el_M_inv))
      get_df1 <- data.frame(dose = el_grid, val = eld_get_val)
      get_df2 <- data.frame(dose = approx_d, val = eld_get_val[(approx_d / 0.0001)+1])
      values$el1$getplot <- ggplot(data = get_df1, aes(x = dose, y = val)) +
        geom_line(stat = "identity") +
        geom_point(data = get_df2, aes(x = dose, y = val), shape = 1, size = 10)
    } else if (criterion == "tau"){
      text = paste("This function does not support tau-optimal designs.")
      values$el1$getplot <- ggplot() + 
        annotate("text", x = 4, y = 25, size = 8, label = text) + 
        theme_void()
    } else if (criterion == "h"){
      text = paste("This function does not support h-optimal designs.")
      values$el1$getplot <- ggplot() + 
        annotate("text", x = 4, y = 25, size = 8, label = text) + 
        theme_void()
    }
    
    compute$el1_approx = TRUE
    
    eff_round_design <- eff_round(approx = approx, model = "ExpLog", criterion = criterion, 
                                  parms = el_par, nPoints = ndp)
    values$el1$eff_round = eff_round_design
    
    compute$el1_effround = TRUE
  })
  
  observeEvent(input$el2_pso, {
    waiter7 <- waiter::Waiter$new(
      id = "el2_pso_out", 
      html = spin_circle(), 
      color = "lightgrey"
    )
    waiter7$show()
    on.exit(waiter7$hide())
    
    c0 = input$el2_c0
    c1 = input$el2_c1
    b0 = input$el2_b0
    b1 = input$el2_b1
    approx_nswarm = input$el_approx_swarm
    approx_iter = input$el_approx_iter
    criterion = input$el_criterion
    lb = input$el_lb
    ub = input$el_ub
    ndp = input$el_dim
    
    el_par <- c(c0, c1, b0, b1)
    psoinfo_approx <- psoinfo_setting(nSwarms = approx_nswarm, Iters = approx_iter)
    approx_res <- pso_approx(model = "ExpLog", criterion = criterion, parms = el_par, 
                             upper = ub, lower = lb, psoinfo_approx = psoinfo_approx)
    
    
    c01 = input$el1_c0
    c11 = input$el1_c1
    b01 = input$el1_b0
    b11 = input$el1_b1
    el_par1 <- c(c01, c11, b01, b11)
    if (criterion == "tau"){
      tau1 <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, c0 = c01, c1 = c11, b0 = b01, b1 = b11)$root
      el_par1 <- c(el_par1, tau1)
    }
    approx_d <- approx_res$design$Support
    approx_w <- approx_res$design$Weight
    approx_w1 <- approx_w[-length(approx_w)]
    approx_val <- values$el1$approx$val
    
    if (criterion == "D"){
      comp_val <- el_doptimal_approx(c(approx_d, approx_w), loc = el_par1)
      approx_res$comp_eff <- (comp_val / approx_val)^(1/4) |> as.numeric()
    }
    if (criterion == "tau"){
      comp_val <- el_tauoptimal_approx(c(approx_d, approx_w), loc = el_par1)
      approx_res$comp_eff <- (approx_val / comp_val) |> as.numeric()
    }
    if (criterion == "h"){
      comp_val <- el_hoptimal_approx(c(approx_d, approx_w), loc = el_par1)
      approx_res$comp_eff <- (approx_val / comp_val) |> as.numeric()
    }
    
    values$el2$approx = approx_res
    
    approx <- approx_res$design |> round_approx()
    approx_d <- approx$Support
    approx_w <- approx$Weight
    el_n <- length(approx_d)
    el_grid <- seq(lb, ub, by = 0.0001)
    el_mat_list <- lapply(1:el_n, function(x) approx_w[x] * el_mat(approx_d[x], c0, c1, b0, b1))
    el_M <- Reduce("+", el_mat_list)
    diag(el_M) <- diag(el_M) + 1e-10
    el_M_inv <- solve(el_M)
    
    if (criterion == "D"){
      eld_get_val <- sapply(el_grid, function(x) eld_get(x, el_par, el_M_inv))
      get_df1 <- data.frame(dose = el_grid, val = eld_get_val)
      get_df2 <- data.frame(dose = approx_d, val = eld_get_val[(approx_d / 0.0001)+1])
      values$el2$getplot <- ggplot(data = get_df1, aes(x = dose, y = val)) +
        geom_line(stat = "identity") +
        geom_point(data = get_df2, aes(x = dose, y = val), shape = 1, size = 10)
    } else if (criterion == "tau"){
      text = paste("This function does not support \n tau-optimal designs.")
      values$hb1$getplot <- ggplot() + 
        annotate("text", x = 4, y = 25, size = 8, label = text) + 
        theme_void()
    } else if (criterion == "h"){
      text = paste("This function does not support \n h-optimal designs.")
      values$el2$getplot <- ggplot() + 
        annotate("text", x = 4, y = 25, size = 8, label = text) + 
        theme_void()
    }
    
    compute$el2_approx = TRUE
    
    eff_round_design <- eff_round(approx = approx, model = "ExpLog", criterion = criterion, 
                                  parms = el_par, nPoints = ndp)
    effround_d <- lapply(1:nrow(eff_round_design$design), 
                         function(x) rep(eff_round_design$design$Support[x], eff_round_design$design$N[x])) |> unlist()
    approx_val <- values$el1$approx$val
    
    if (criterion == "D"){
      comp_val <- el_doptimal(effround_d, loc = el_par1)
      eff_round_design$comp_eff <- (comp_val / approx_val)^(1/4) |> as.numeric()
    }
    if (criterion == "tau"){
      comp_val <- el_tauoptimal(effround_d, loc = el_par1)
      eff_round_design$comp_eff <- (approx_val / comp_val) |> as.numeric()
    }
    if (criterion == "h"){
      comp_val <- el_hoptimal(effround_d, loc = el_par1)
      eff_round_design$comp_eff <- (approx_val / comp_val) |> as.numeric()
    }
    
    eff_round_design$comp_eff_erm <- eff_round_design$comp_eff / values$el1$eff_round$efficiency
    values$el2$eff_round = eff_round_design 
    
    compute$el2_effround = TRUE
  })
  
  observeEvent(input$el_comp, {
    if(compute$el1_approx && compute$el2_approx){
      el1_c0 = input$el1_c0
      el1_c1 = input$el1_c1
      el1_b0 = input$el1_b0
      el1_b1 = input$el1_b1
      el2_c0 = input$el2_c0
      el2_c1 = input$el2_c1
      el2_b0 = input$el2_b0
      el2_b1 = input$el2_b1
      approx_nswarm = input$el_approx_swarm
      approx_iter = input$el_approx_iter
      criterion = input$el_criterion
      lb = input$el_lb
      ub = input$el_ub
      ndp = input$el_dim
      
      el_par1 <- c(el1_c0, el1_c1, el1_b0, el1_b1)
      
      approx_d <- values$el2$approx$design$Support
      approx_w <- values$el2$approx$design$Weight[-length(values$el2$approx$design$Weight)]
      el1_approx_val <- values$el1$approx$val
      el1_erm_eff <- values$el1$eff_round$efficiency
      el2_erm_eff <- values$el2$eff_round$efficiency
      if (criterion == "D"){
        el2_approx_val <- el_doptimal_approx(c(approx_d, approx_w), loc = el_par1)
        approx_comp_eff <- (el2_approx_val / el1_approx_val)^(1/4) |> as.numeric()
        erm_approx_comp_eff <- el2_erm_eff * approx_comp_eff |> as.numeric()
        erm_comp_eff <- erm_approx_comp_eff / el1_erm_eff |> as.numeric()
      }
      if (criterion == "tau"){
        tau <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, c0 = el1_c0, c1 = el1_c1, b0 = el1_b0, b1 = el1_b1)$root
        el_par1 <- c(el_par1, tau)
        el2_approx_val <- el_tauoptimal_approx(c(approx_d, approx_w), loc = el_par1)
        approx_comp_eff <- (el1_approx_val / el2_approx_val) |> as.numeric()
        erm_approx_comp_eff <- el2_erm_eff * approx_comp_eff |> as.numeric()
        erm_comp_eff <- erm_approx_comp_eff / el1_erm_eff |> as.numeric()
      }
      if (criterion == "h"){
        el2_approx_val <- el_hoptimal_approx(c(approx_d, approx_w), loc = el_par1)
        approx_comp_eff <- (el1_approx_val / el2_approx_val) |> as.numeric()
        erm_approx_comp_eff <- el2_erm_eff * approx_comp_eff |> as.numeric()
        erm_comp_eff <- erm_approx_comp_eff / el1_erm_eff |> as.numeric()
      }
      
      values$el1$comp$approx = approx_comp_eff
      values$el1$comp$erm_approx = erm_approx_comp_eff
      values$el1$comp$erm = erm_comp_eff
      
      compute$el_comp = TRUE
    }
  })
  
  observeEvent(input$el_spec, {
    if(compute$el1_approx && compute$el2_approx){
      el1_c0 = input$el1_c0
      el1_c1 = input$el1_c1
      el1_b0 = input$el1_b0
      el1_b1 = input$el1_b1
      el2_c0 = input$el2_c0
      el2_c1 = input$el2_c1
      el2_b0 = input$el2_b0
      el2_b1 = input$el2_b1
      approx_nswarm = input$el_approx_swarm
      approx_iter = input$el_approx_iter
      criterion = input$el_criterion
      lb = input$el_lb
      ub = input$el_ub
      ndp = input$el_dim
      
      el_par1 <- c(el1_c0, el1_c1, el1_b0, el1_b1)
      el_par2 <- c(el2_c0, el2_c1, el2_b0, el2_b1)
      
      sup_error <- FALSE
      N_error <- FALSE
      tryCatch({support <- input$el_spec_sup |> strsplit(",") |> unlist() |> as.numeric()}, 
               warning = function(w) {sup_error <<- TRUE}, 
               error = function(e) {sup_error <<- TRUE})
      if (sup_error == FALSE) spec_error$el_support = FALSE
      else if (sup_error == TRUE) spec_error$el_support = TRUE
      
      tryCatch({N <- input$el_spec_N |> strsplit(",") |> unlist() |> as.numeric()}, 
               warning = function(w) {N_error <<- TRUE}, 
               error = function(e) {N_error <<- TRUE})
      if (N_error == FALSE) spec_error$el_N = FALSE
      else if (N_error == TRUE) spec_error$el_N = TRUE
      
      if(spec_error$el_support == FALSE & spec_error$el_N == FALSE){
        if (length(support) != length(N)){
          spec_error$el_align <- TRUE
        } else {
          spec_error$el_align <- FALSE
        }
        
        if (criterion == "D" & length(support) < 4){
          spec_error$el_dop <- TRUE
        } else {
          spec_error$el_dop <- FALSE
        }
        
        if(sum(N) != ndp){
          spec_error$el_Nrep <- TRUE
        } else {
          spec_error$el_Nrep <- FALSE
        }
        
        if(any(support < lb) || any(support > ub)){
          spec_error$el_space = TRUE
        }
      } 
      if(!(TRUE %in% c(spec_error$el_support, spec_error$el_N, spec_error$el_align))){
        el_dose <- lapply(1:length(N), function(x) rep(support[x], N[x])) |> unlist()
        
        if(criterion == "D"){
          spec_val1 <- el_doptimal(el_dose, el_par1)
          spec_val2 <- el_doptimal(el_dose, el_par2)
          approx1 <- values$el1$approx$val
          approx2 <- values$el2$approx$val
          eff1 <- (spec_val1 / approx1)^(1/4)
          eff2 <- (spec_val2 / approx2)^(1/4)
        } else if (criterion == "tau"){
          tau1 <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, c0 = el1_c0, c1 = el1_c1, b0 = el1_b0, b1 = el1_b1)$root
          tau2 <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, c0 = el2_c0, c1 = el2_c1, b0 = el2_b0, b1 = el2_b1)$root
          el_par1 <- c(el_par1, tau1)
          el_par2 <- c(el_par2, tau2)
          spec_val1 <- el_tauoptimal(el_dose, el_par1)
          spec_val2 <- el_tauoptimal(el_dose, el_par2)
          approx1 <- values$el1$approx$val
          approx2 <- values$el2$approx$val
          eff1 <- approx1 / spec_val1
          eff2 <- approx2 / spec_val2
        } else if (criterion == "h"){
          spec_val1 <- el_hoptimal(el_dose, el_par1)
          spec_val2 <- el_hoptimal(el_dose, el_par2)
          approx1 <- values$el1$approx$val
          approx2 <- values$el2$approx$val
          eff1 <- approx1 / spec_val1
          eff2 <- approx2 / spec_val2
        }
        
        values$el2$spec$design <- data.frame(Support = support, N = N)
        values$el2$spec$efficiency$approx1 <- eff1
        values$el2$spec$efficiency$erm1 <- eff1 /  values$el2$eff_round$efficiency
        values$el2$spec$efficiency$approx2 <- eff2
        values$el2$spec$efficiency$erm2 <- eff2 /  values$el2$eff_round$efficiency
        values$el2$spec$design <- data.frame(Support = support, N = N)
        
        compute$el_spec = TRUE
      }
    }
  })
  
  output$el1_pso_out <- renderPrint({
    if(compute$el1_approx == FALSE){
      cat("Press the button above to compute!")
    } else {
      cat("Approximate Design: ", "\n")
      approx_df <- values$el1$approx$design |> round_approx()
      print(approx_df, row.names = FALSE)
      cat("\nEfficient Rounding Method (ERM)-generated Design: ", "\n")
      erm_df <- values$el1$eff_round$design |> round(4)
      print(erm_df, row.names = F)
      cat("\nEfficiency of the ERM-generated design relative \nto the approximate design for Model 1: ", values$el1$eff_round$efficiency |> round(4))
      cat("\n\nThe value of the design efficiency represents its worth; 
a design with 0.5 efficiency means that it has to be 
replicated twice to do as well as the optimum design. 
Similar interpretation applies when the design is compared 
to a non-optimal design.")
    }
    
  })
  
  output$el2_pso_out <- renderPrint({
    if(compute$el2_approx == FALSE) {
      cat("Press the button above to compute!")
    } else {
      cat("Approximate Design: ", "\n")
      approx_df <- values$el2$approx$design |> round_approx()
      print(approx_df, row.names = FALSE)
      cat("\n\nEfficient Rounding Method (ERM)-generated Design: ", "\n")
      erm_df <- values$el2$eff_round$design |> round(4)
      print(erm_df, row.names = FALSE)
      cat("\nEfficiency of the ERM-generated design relative \nto the approximate design for Model 2: ", 
          values$el1$eff_round$efficiency |> round(4))
    }
    
  })
  
  output$el_spec_out <- renderPrint({
    if(compute$el1_approx == FALSE){
      cat("Please compute the optimal design for model 1 first!")
    } else if (compute$el2_approx == FALSE){
      cat("Please compute the optimal design for model 2 first!")
    } else if (compute$el_spec == FALSE){
      cat("Press the button above to compute!")
      cat("\nEach number must be seperated by a single\ncomma. e.g. 0,0.0375,0.075,0.1125,0.15")
    } else if (spec_error$el_support){
      cat("Error in the format of Support.\nEach number must be seperated by a single\ncomma. e.g. 2,2,2,2,2")
    } else if (spec_error$el_N){
      cat("Error in the format of N.\nEach number must be seperated by a single comma.\ne.g. 2,2,3")
    } else if (spec_error$el_align){
      cat("The number of support points must be\nthe same with the number of replicates!")
    } else {
      if (spec_error$el_dop) cat("Warning:\nFor D-optimal, at least 4 support is needed\nto have efficiency > 0.\n\n")
      if (spec_error$el_Nrep) cat("Warning:\nThe total number of replicates you entered\ndoes not equal to the optimal approximate\ndesign.\n\n")
      if (spec_error$el_space) cat("Warning:\nSome dose levels are outside the design space, please check!\n\n")
      cat("User Specified Design:\n")
      print(values$el2$spec$design, row.names = FALSE)
      cat("\nEfficiency of the user specified design relative \nto the approximate design for Model 1: ", 
          values$el2$spec$efficiency$approx1 |> round(4))
      cat("\n\nEfficiency of the user specified design relative \nto the ERM-generated design for Model 1: ", 
          values$el2$spec$efficiency$erm1 |> round(4))
      cat("\n\nEfficiency of the user specified design relative \nto the approximate design for Model 2: ", 
          values$el2$spec$efficiency$approx2 |> round(4))
      cat("\n\nEfficiency of the user specified design relative \nto the ERM-generated design for Model 2: ", 
          values$el2$spec$efficiency$erm2 |> round(4))
    }
    
  })
  
  output$el_comp_out <- renderPrint({
    if(compute$el1_approx == FALSE){
      cat("Please compute the optimal design for model 1 first!")
    } else if (compute$el2_approx == FALSE){
      cat("Please compute the optimal design for model 2 first!")
    } else if (compute$el_comp == FALSE){
      cat("Press the button above to compute!")
    } else {
      cat("Efficiency of the approximate design for Model 2 relative \nto the approximate design for Model 1: ", 
          values$el1$comp$approx |> round(4))
      cat("\n\nEfficiency of the ERM-generated design for Model 2 relative \nto the approximate design for Model 1: ", 
          values$el1$comp$erm_approx |> round(4))
      cat("\n\nEfficiency of the ERM-generated design for Model 2 relative \nto the ERM-generated design for Model 1: ", 
          values$el1$comp$erm |> round(4))
    }
    
  })
  
  output$el_get_plot <- renderPlot({
    md = input$el_get_md
    if (md == "md1") getplot = values$el1$getplot
    else if (md == "md2") getplot = values$el2$getplot
    getplot
  })
}



### Hunt-Bowman model

## Hunt-Bowman parameters
hb_parms <- function(c1, tau, b0, b1){
  c(c1, tau, b0, b1)
}

# Hunt-Bowman function
hunt_bowman <- function(d, c1, tau, b0, b1){
  if (d <= tau){
    res = c1 * d^2 - c1 * tau * d + 1 / (1 + exp(b0))
  } else {
    res = 1 / (1 + exp(b0 - b1 * (d - tau)))
  }
}

# Dose-response plot for Hunt-Bowman models
hunt_bowman_plot <- function(parms, upper_bound){
  c1 = parms[1]
  tau = parms[2]
  b0 = parms[3]
  b1 = parms[4]
  
  fp <- seq(0, upper_bound, by = (upper_bound - 0)/100)
  hb_df <- data.frame(dose = fp, response = sapply(fp, function(x) hunt_bowman(x, c1, tau, b0, b1)))
  ggplot(data = hb_df, aes(x = dose, y = response)) + 
    geom_line()
}

# Hunt-Bowman information matrix
hb_f <- function(d, c1, tau, b0, b1){
  if (d <= tau){
    f1 <- d^2 - tau * d
    f2 <- -c1 * d
    f3 <- -exp(b0) / (1 + exp(b0))^2
    f4 <- 0
    f <- matrix(c(f1, f2, f3, f4))
  } 
  else if (d > tau){
    f1 <- 0
    f2 <- -b1 * exp(b0 - b1 * (d - tau)) / (1 + exp(b0 - b1 * (d - tau)))^2
    f3 <- -exp(b0 - b1 * (d - tau)) / (1 + exp(b0 - b1 * (d - tau)))^2
    f4 <- (d - tau) * exp(b0 - b1 * (d - tau)) / (1 + exp(b0 - b1 * (d - tau)))^2
    f <- matrix(c(f1, f2, f3, f4))
  }
  
  f
}

hb_grad <- function(d, theta){
  c1 = theta[1]
  tau = theta[2]
  b0 = theta[3]
  b1 = theta[4]
  
  if (d <= tau){
    f1 <- d^2 - tau * d
    f2 <- -c1 * d
    f3 <- -exp(b0) / (1 + exp(b0))^2
    f4 <- 0
    return(c(f1, f2, f3, f4))
  } 
  else if (d > tau){
    f1 <- 0
    f2 <- -b1 * exp(b0 - b1 * (d - tau)) / (1 + exp(b0 - b1 * (d - tau)))^2
    f3 <- -exp(b0 - b1 * (d - tau)) / (1 + exp(b0 - b1 * (d - tau)))^2
    f4 <- (d - tau) * exp(b0 - b1 * (d - tau)) / (1 + exp(b0 - b1 * (d - tau)))^2
    return(c(f1, f2, f3, f4))
  }
}

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

# Hunt-Bowman Exact Design Criteria
hb_doptimal <- function(d, loc){
  
  # Hunt-Bowman parameters
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  
  # Number of experiment points
  n <- length(d)
  
  # Evaluate d-optimality criterion value
  mat_list <- lapply(d, function(x) 1/n * hb_mat(x, c1, tau, b0, b1))
  inf_mat <- Reduce("+", mat_list)
  -det(inf_mat)
}

hb_doptimal_approx <- function(input, loc){
  
  # Hunt-Bowman parameters
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  pen = 9e+10
  
  n <- (length(input)+1)/2
  d <- input[1:n]
  w <- input[(n+1):(2*n-1)]
  w[n] <- 1 - sum(w)
  
  # Evaluate d-optimality criterion value
  if (w[n] < 0) res <- pen
  else {
    mat_list <- lapply(1:n, function(x) w[x] * hb_mat(d[x], c1, tau, b0, b1))
    inf_mat <- Reduce("+", mat_list)
    res <- -det(inf_mat)
  }
  
  res
}

hbd_get <- function(x, parms, M_inv){
  c1 <- parms[1]
  tau <- parms[2]
  b0 <- parms[3]
  b1 <- parms[4]
  
  #print(M)
  fx <- hb_f(x, c1, tau, b0, b1)
  
  t(fx) %*% M_inv %*% fx
}

hb_tauoptimal <- function(d, loc){
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  n <- length(d)
  pen = 9e+10
  
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

hb_tauoptimal_approx <- function(input, loc){
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  n <- (length(input)+1)/2
  pen = 9e+10
  
  d <- input[1:n]
  w <- input[(n+1):(2*n-1)]
  w[n] <- 1 - sum(w)
  
  mat_list <- lapply(1:n, function(x) w[x] * hb_mat(d[x], c1, tau, b0, b1))
  M <- Reduce("+", mat_list)
  diag(M) <- diag(M) + 1e-10
  b <- matrix(c(0, 1, 0, 0))
  
  if (rcond(M) < 2.220446e-16 || w[n] < 0){
    res = pen
  } 
  else{
    M_inv <- solve(M)
    res <- t(b) %*% M_inv %*% b
  }
  
  res
}

hbtau_get <- function(x, parms, M_inv){
  c1 <- parms[1]
  tau <- parms[2]
  b0 <- parms[3]
  b1 <- parms[4]
  
  fx <- hb_f(x, c1, tau, b0, b1)
  b <- matrix(c(0, 1, 0, 0))
  
  (t(fx) %*% M_inv %*% b)^2 - t(b) %*% M_inv %*% b
}

hb_hoptimal <- function(d, loc){
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  n <- length(d)
  pen = 9e+10
  
  mat_list <- lapply(d, function(x) 1/n * hb_mat(x, c1, tau, b0, b1))
  M <- Reduce("+", mat_list)
  diag(M) <- diag(M) + 1e-10
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
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  n <- (length(input)+1)/2
  pen = 9e+10
  
  d <- input[1:n]
  w <- input[(n+1):(2*n-1)]
  w[n] <- 1 - sum(w)
  
  mat_list <- lapply(1:n, function(x) w[x] * hb_mat(d[x], c1, tau, b0, b1))
  M <- Reduce("+", mat_list)
  diag(M) <- diag(M) + 1e-10
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

hbh_get <- function(x, parms, M_inv){
  c1 <- parms[1]
  tau <- parms[2]
  b0 <- parms[3]
  b1 <- parms[4]
  
  fx <- hb_f(x, c1, tau, b0, b1)
  h <- matrix(c(-tau, -c1, 0, 0))
  
  (t(fx) %*% M_inv %*% h)^2 - t(h) %*% M_inv %*% h
}

### exp-log model

# exp-log model parameters
el_parms <- function(c0, c1, b0, b1){
  c(c0, c1, b0, b1)
}

# exp-log function
exp_log <- function(d, c0, c1, b0, b1){
  c0 * exp(-c1 * d) + 1 / (1 + exp(b0 - b1 * d))
}

# Dose-response plot for exp-log models
exp_log_plot <- function(parms, upper_bound){
  c0 = parms[1]
  c1 = parms[2]
  b0 = parms[3]
  b1 = parms[4]
  
  fp <- seq(0, upper_bound, by = (upper_bound - 0)/100)
  el_df <- data.frame(dose = fp, response = sapply(fp, function(x) exp_log(x, c0, c1, b0, b1)))
  ggplot(data = el_df, aes(x = dose, y = response)) + 
    geom_line()
}

# exp-log model information matrix
exp_log_f <- function(d, c0, c1, b0, b1){
  f1 <- exp(-c1 * d)
  f2 <- -c0 * d * exp(-c1 * d)
  f3 <- -exp(b0 - b1 * d) / (1 + exp(b0 - b1 * d))^2
  f4 <- d * exp(b0 - b1 * d) / (1 + exp(b0 - b1 * d))^2
  f <- matrix(c(f1, f2, f3, f4))
  f
}

exp_log_grad <- function(d, theta){
  c0 = theta[1] 
  c1 = theta[2] 
  b0 = theta[3] 
  b1 = theta[4]
  
  f1 <- exp(-c1 * d)
  f2 <- -c0 * d * exp(-c1 * d)
  f3 <- -exp(b0 - b1 * d) / (1 + exp(b0 - b1 * d))^2
  f4 <- d * exp(b0 - b1 * d) / (1 + exp(b0 - b1 * d))^2
  f <- c(f1, f2, f3, f4)
  f
}

el_mat <- function(d, c0, c1, b0, b1){
  n <- length(d)
  mat_list <- lapply(d, function(x) 1/n * exp_log_f(x, c0, c1, b0, b1) %*% t(exp_log_f(x, c0, c1, b0, b1)))
  
  M <- Reduce("+", mat_list)
  M
}

## exp-log D-optimal
el_doptimal <- function(d, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  
  M <- el_mat(d, c0, c1, b0, b1)
  
  -det(M)
}

el_doptimal_approx <- function(input, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  n <- (length(input)+1)/2
  pen = 9e+10
  
  d <- input[1:n]
  w <- input[(n+1):(2*n-1)]
  w[n] <- 1 - sum(w)
  
  mat_list <- lapply(1:n, function(x) w[x] * exp_log_f(d[x], c0, c1, b0, b1) %*% t(exp_log_f(d[x], c0, c1, b0, b1)))
  M <- Reduce("+", mat_list)
  
  if(w[n] < 0) res = pen
  else res = -det(M)
  
  res
}

eld_get <- function(x, parms, M_inv){
  c0 = parms[1]
  c1 = parms[2]
  b0 = parms[3]
  b1 = parms[4]
  
  fx <- exp_log_f(x, c0, c1, b0, b1)
  t(fx) %*% M_inv %*% fx
}

## exp-log model h-optimal
exp_log_h <- function(d, c0, c1, b0, b1){
  h1 <- -c1
  h2 <- -c0
  h3 <- b1 * exp(b0) * (1 - exp(b0)) / (exp(b0) + 1)^3
  h4 <- exp(b0) / (exp(b0) + 1)^2 
  h <- matrix(c(h1, h2, h3, h4))
  
  h
}

el_hoptimal <- function(d, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  pen = 9e+10
  
  n <- length(d)
  h1 <- -c1
  h2 <- -c0
  h3 <- b1 * exp(b0) * (1 - exp(b0)) / (exp(b0) + 1)^3
  h4 <- exp(b0) / (exp(b0) + 1)^2 
  h <- matrix(c(h1, h2, h3, h4))
  
  M <- el_mat(d, c0, c1, b0, b1)
  
  if (rcond(M) < 2.220446e-16){
    res = pen
  }
  else {
    M_inv <- solve(M)
    res = t(h) %*% M_inv %*% h
  }
  
  res
}

el_hoptimal_approx <- function(input, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  n = (length(input)+1)/2
  pen = 9e+10
  
  d <- input[1:n]
  w <- input[(n+1):(2*n-1)]
  w[n] <- 1 - sum(w)
  
  h1 <- -c1
  h2 <- -c0
  h3 <- b1 * exp(b0) * (1 - exp(b0)) / (exp(b0) + 1)^3
  h4 <- exp(b0) / (exp(b0) + 1)^2 
  h <- matrix(c(h1, h2, h3, h4))
  
  mat_list <- lapply(1:n, function(x) w[x] * exp_log_f(d[x], c0, c1, b0, b1) %*% t(exp_log_f(d[x], c0, c1, b0, b1)))
  M <- Reduce("+", mat_list)
  
  if (rcond(M) < 2.220446e-16 || w[n] < 0){
    res = pen
  }
  else {
    M_inv <- solve(M)
    res = t(h) %*% M_inv %*% h
  }
  
  res
}

elh_get <- function(x, parms, M_inv){
  c0 <- parms[1]
  c1 <- parms[2]
  b0 <- parms[3]
  b1 <- parms[4]
  
  fx <- exp_log_f(x, c0, c1, b0, b1)
  h <- exp_log_h(x, c0, c1, b0, b1)
  
  (t(fx) %*% M_inv %*% h)^2 - t(h) %*% M_inv %*% h
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

el_tauoptimal <- function(d, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  tau <- loc[5]
  pen = 9e+10
  
  n <- length(d)
  w <- rep(1/n, n)
  
  b <- exp_log_b(tau, c0, c1, b0, b1)
  
  M <- el_mat(d, c0, c1, b0, b1)
  diag(M) <- diag(M) + 1e-10
  if (rcond(M) < 2.220446e-16){
    res = pen
  }
  else {
    M_inv <- solve(M)
    res = t(b) %*% M_inv %*% b
  }
  
  res
}

el_tauoptimal_approx <- function(input, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  tau = loc[5]
  pen = 9e+10
  
  n <- (length(input)+1)/2
  
  d <- input[1:n]
  w <- input[(n+1):(2*n-1)]
  w[n] <- 1 - sum(w)
  
  b <- exp_log_b(tau, c0, c1, b0, b1)
  
  mat_list <- lapply(1:n, function(x) w[x] * exp_log_f(d[x], c0, c1, b0, b1) %*% t(exp_log_f(d[x], c0, c1, b0, b1)))
  M <- Reduce("+", mat_list)
  diag(M) <- diag(M) + 1e-10
  
  if (rcond(M) < 2.220446e-16 || w[n] < 0){
    res = pen
  }
  else {
    M_inv <- solve(M)
    res = t(b) %*% M_inv %*% b
  }
  
  res
}

eltau_get <- function(x, parms, M_inv){
  c0 <- parms[1]
  c1 <- parms[2]
  b0 <- parms[3]
  b1 <- parms[4]
  tau <- uniroot(tau_func, c(1e-10, 0.15), tol = 1e-10, 
                 c0 = c0, c1 = c1, b0 = b0, b1 = b1)$root
  
  fx <- exp_log_f(x, c0, c1, b0, b1)
  b <- exp_log_b(tau, c0, c1, b0, b1)
  
  (t(fx) %*% M_inv %*% b)^2 - t(b) %*% M_inv %*% b
}

### PSO function

## Initalize PSO settings

psoinfo_setting <- function(nSwarms = 64, Iters = 1000){
  getPSOInfo(nSwarm = nSwarms, maxIter = Iters, w0 = 0.9, w1 = 0.4, w_var = 1)
}

# PSO for exact design

hormesis_pso <- function(model, criterion, parms, upper, lower, nPoints, nRep = 1, psoinfo_exact, psoinfo_approx){
  start_time <- Sys.time()
  if (model == "HuntBowman"){
    nDim = 4
    if (criterion == "D"){
      obj_exact <- hb_doptimal
      obj_approx <- hb_doptimal_approx
    } 
    else if(criterion == "tau"){
      obj_exact <- hb_tauoptimal
      obj_approx <- hb_tauoptimal_approx
      nDim = 2
    } 
    else if (criterion == "h"){
      obj_exact <- hb_hoptimal
      obj_approx <- hb_hoptimal_approx
    } 
  } 
  else if (model == "ExpLog"){
    nDim = 4
    if (criterion == "D"){
      obj_approx <- el_doptimal_approx
      obj_exact <- el_doptimal
    } 
    else if (criterion == "tau"){
      obj_approx <- el_tauoptimal_approx
      obj_exact <- el_tauoptimal
      tau <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, 
                     c0 = parms[1], c1 = parms[2], b0 = parms[3], b1 = parms[4])$root
      parms <- c(parms, tau)
      nDim = 2
    } 
    else if (criterion == "h"){
      obj_approx <- el_hoptimal_approx
      obj_exact <- el_hoptimal
    } 
  }
  
  #psoinfo_approx <- psoinfo_setting(256, 2000)
  approx_lb <- c(rep(lower, nDim), rep(0, nDim-1))
  approx_ub <- c(rep(upper, nDim), rep(1, nDim-1))
  approx_result <- globpso(objFunc = obj_approx, lower = approx_lb, upper = approx_ub, PSO_INFO = psoinfo_approx, 
                           loc = parms, verbose = F)

  approx_d <- approx_result$par[1:nDim]
  approx_w <- approx_result$par[(nDim+1):(2*nDim-1)]
  approx_w <- c(approx_w, 1 - sum(approx_w))
  approx_design <- data.frame(support = round(approx_d, 4), weight = round(approx_w, 3))
  approx_design <- approx_design[order(approx_design$support),] 
  approx_val <- approx_result$val
  
  idx0 <- which(approx_design$weight == 0)
  if(length(idx0) != 0) approx_design <- approx_design[-idx0, ]
  
  if (nrow(approx_design) != length(unique(approx_design$support))){
    cnt <- approx_design |> count(support)
    ext <- cnt$support[which(cnt$n != 1)]
    apprep <- approx_design$support == ext
    w <- sum(approx_design$weight[apprep])
    approx_design <- approx_design[!apprep,]
    approx_design <- rbind(approx_design, c(ext, w))
  }
  
  effr <- efficient.rounding(approx_w, nPoints)
  mu_pso <- lapply(1:length(effr), function(x) rep(approx_d[x], effr[x])) |> unlist()
  mvnorm_sd = (upper - lower)/4
  nswarm <- psoinfo_exact$nSwarm/2
  mvnorm_pso <- mvrnorm(n = (nswarm-1), mu = mu_pso, Sigma = diag(nPoints) * mvnorm_sd)
  mvnorm_pso[mvnorm_pso < 0] <- 0
  mvnorm_pso[mvnorm_pso > 0.15] <- 0.15
  mvnorm_pso <- rbind(mvnorm_pso, mu_pso)
  
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- globpso(objFunc = obj_exact, lower = rep(lower, nPoints), 
                                                                        upper = rep(upper, nPoints), init = mvnorm_pso, 
                                                                        PSO_INFO = psoinfo_exact, loc = parms, verbose = F))
  val_list <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  best_idx <- which.min(val_list)
  exact_design <- pso_results[[best_idx]]$par |> round(4) |> table() |> data.frame()
  colnames(exact_design) <- c("Support", "N")
  exact_val <- pso_results[[best_idx]]$val
  
  # pso_mvnorm <- globpso(objFunc = obj_exact, lower = rep(lower, nPoints), upper = rep(upper, nPoints), 
  #                       init = mvnorm_pso, PSO_INFO = psoinfo, loc = parms, verbose = F)
  # exact_design <- pso_mvnorm$par |> round(4) |> table() |> data.frame()
  # colnames(exact_design) <- c("Support", "N")
  # exact_val <- pso_mvnorm$val
  
  eff_round_res <- mu_pso |> round(4) |> table() |> data.frame()
  colnames(eff_round_res) <- c("Support", "N")
  if (model == "HuntBowman"){
    if (criterion == "D") eff_round_eff <- (hb_doptimal(d = mu_pso, loc = parms) / approx_val)^(1/4) |> round(4)
    else if(criterion == "tau") eff_round_eff <- (approx_val / hb_tauoptimal(d = mu_pso, loc = parms)) |> round(4)
    else if (criterion == "h") eff_round_eff <- (approx_val / hb_hoptimal(d = mu_pso, loc = parms)) |> round(4)
  }
  else if (model == "ExpLog"){
    if (criterion == "D") eff_round_eff <- (el_doptimal(d = mu_pso, loc = parms) / approx_val) ^ (1/4) |> round(4)
    else if (criterion == "tau") eff_round_eff <- (approx_val / el_tauoptimal(d = mu_pso, loc = parms)) |> round(4)
    else if (criterion == "h") eff_round_eff <- (approx_val / el_hoptimal(d = mu_pso, loc = parms)) |> round(4)
  }
  
  nEff <- length(parms)
  if (criterion == "D") eff <- (exact_val/approx_val) ^ (1 / nEff)
  else eff <- (approx_val / exact_val)
  eff = round(eff, 4)
  end_time <- Sys.time()
  
  result <- list(approx_design = approx_design, exact_design = exact_design, 
                 efficiency = eff, eff_round_design = eff_round_res, eff_round_eff = eff_round_eff,
                 runtime = end_time - start_time)
  result
}

pso_approx <- function(model, criterion, parms, upper, lower, psoinfo_approx){
  if (model == "HuntBowman"){
    nDim = 4
    if (criterion == "D"){
      obj_exact <- hb_doptimal
      obj_approx <- hb_doptimal_approx
    } 
    else if(criterion == "tau"){
      obj_exact <- hb_tauoptimal
      obj_approx <- hb_tauoptimal_approx
      nDim = 2
    } 
    else if (criterion == "h"){
      obj_exact <- hb_hoptimal
      obj_approx <- hb_hoptimal_approx
    } 
  } 
  else if (model == "ExpLog"){
    nDim = 4
    if (criterion == "D"){
      obj_approx <- el_doptimal_approx
      obj_exact <- el_doptimal
    } 
    else if (criterion == "tau"){
      obj_approx <- el_tauoptimal_approx
      obj_exact <- el_tauoptimal
      tau <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, 
                     c0 = parms[1], c1 = parms[2], b0 = parms[3], b1 = parms[4])$root
      parms <- c(parms, tau)
      nDim = 2
    } 
    else if (criterion == "h"){
      obj_approx <- el_hoptimal_approx
      obj_exact <- el_hoptimal
    } 
  }
  
  #psoinfo_approx <- psoinfo_setting(256, 2000)
  approx_lb <- c(rep(lower, nDim), rep(0, nDim-1))
  approx_ub <- c(rep(upper, nDim), rep(1, nDim-1))
  approx_result <- globpso(objFunc = obj_approx, lower = approx_lb, upper = approx_ub, PSO_INFO = psoinfo_approx, 
                           loc = parms, verbose = T)
  
  approx_d <- approx_result$par[1:nDim]
  approx_w <- approx_result$par[(nDim+1):(2*nDim-1)]
  approx_w <- c(approx_w, 1 - sum(approx_w))
  approx_design <- data.frame(Support = approx_d, Weight = approx_w)
  approx_design <- approx_design[order(approx_design$Support),] 
  approx_val <- approx_result$val
  
  # idx0 <- which(round(approx_design$weight, 3) == 0)
  # if(length(idx0) != 0) approx_design <- approx_design[-idx0, ]
  # 
  # if (nrow(approx_design) != length(unique(approx_design$Support))){
  #   cnt <- approx_design |> count(support)
  #   ext <- cnt$Support[which(cnt$n != 1)]
  #   apprep <- approx_design$Support == ext
  #   w <- sum(approx_design$Weight[apprep])
  #   approx_design <- approx_design[!apprep,]
  #   approx_design <- rbind(approx_design, c(ext, w))
  # }
  
  approx_res <- list(design = approx_design, val = approx_val)
  approx_res
}

eff_round <- function(approx, model, criterion, parms, nPoints){
  approx_d <- approx$Support
  approx_w <- approx$Weight
  approx_w1 <- approx_w[-length(approx_w)]
  
  if(sum(approx_w) != 1) print(approx_w)
  round_n <- efficient.rounding(approx_w, nPoints)
  eff_design <- data.frame(Support = approx_d, N = round_n)
  eff_d <- lapply(1:length(approx_d), function(x) rep(eff_design$Support[x], eff_design$N[x])) |> unlist()
  
  if (model == "HuntBowman"){
    if (criterion == "D"){
      approx_val <- hb_doptimal_approx(input = c(approx_d, approx_w1), loc = parms)
      eff_val <- hb_doptimal(d = eff_d, loc = parms)
      eff <- (eff_val / approx_val)^(1/4)
    } 
    else if(criterion == "tau"){
      approx_val <- hb_tauoptimal_approx(input = c(approx_d, approx_w1), loc = parms)
      eff_val <- hb_tauoptimal(d = eff_d, loc = parms)
      eff <- (approx_val / eff_val)
    } 
    else if (criterion == "h"){
      approx_val <- hb_hoptimal_approx(input = c(approx_d, approx_w1), loc = parms)
      eff_val <- hb_hoptimal(d = eff_d, loc = parms)
      eff <- (approx_val / eff_val)
    } 
  } 
  else if (model == "ExpLog"){
    nDim = 4
    if (criterion == "D"){
      approx_val <- el_doptimal_approx(input = c(approx_d, approx_w1), loc = parms)
      eff_val <- el_doptimal(d = eff_d, loc = parms)
      eff <- (eff_val / approx_val)^(1/4)
    } 
    else if (criterion == "tau"){
      tau <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, 
                     c0 = parms[1], c1 = parms[2], b0 = parms[3], b1 = parms[4])$root
      parms <- c(parms, tau)
      approx_val <- el_tauoptimal_approx(input = c(approx_d, approx_w1), loc = parms)
      eff_val <- el_tauoptimal(d = eff_d, loc = parms)
      eff <- (approx_val / eff_val)
    } 
    else if (criterion == "h"){
      approx_val <- el_hoptimal_approx(input = c(approx_d, approx_w1), loc = parms)
      eff_val <- el_hoptimal(d = eff_d, loc = parms)
      eff <- (approx_val / eff_val)
    } 
  }
  
  eff_res <- list(design = eff_design, efficiency = eff)
}

pso_exact <- function(approx, model, criterion, parms, upper, lower, psoinfo_exact, nPoints){
  if (model == "HuntBowman"){
    nDim = 4
    if (criterion == "D"){
      obj_exact <- hb_doptimal
      obj_approx <- hb_doptimal_approx
    } 
    else if(criterion == "tau"){
      obj_exact <- hb_tauoptimal
      obj_approx <- hb_tauoptimal_approx
      nDim = 2
    } 
    else if (criterion == "h"){
      obj_exact <- hb_hoptimal
      obj_approx <- hb_hoptimal_approx
    } 
  } 
  else if (model == "ExpLog"){
    nDim = 4
    if (criterion == "D"){
      obj_approx <- el_doptimal_approx
      obj_exact <- el_doptimal
    } 
    else if (criterion == "tau"){
      obj_approx <- el_tauoptimal_approx
      obj_exact <- el_tauoptimal
      tau <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, 
                     c0 = parms[1], c1 = parms[2], b0 = parms[3], b1 = parms[4])$root
      parms <- c(parms, tau)
      nDim = 2
    } 
    else if (criterion == "h"){
      obj_approx <- el_hoptimal_approx
      obj_exact <- el_hoptimal
    } 
  }
  
  approx_d <- approx$design$Support
  approx_w <- approx$design$Weight
  approx_w1 <- approx_w[-length(approx_w)]
  approx_val <- approx$val
  
  effr <- efficient.rounding(approx_w, nPoints)
  mu_pso <- lapply(1:length(effr), function(x) rep(approx_d[x], effr[x])) |> unlist()
  mvnorm_sd = (upper - lower)/4
  nswarm <- psoinfo_exact$nSwarm/2
  mvnorm_pso <- mvrnorm(n = (nswarm-1), mu = mu_pso, Sigma = diag(nPoints) * mvnorm_sd)
  mvnorm_pso[mvnorm_pso < lower] <- lower
  mvnorm_pso[mvnorm_pso > upper] <- upper
  mvnorm_pso <- rbind(mvnorm_pso, mu_pso)
  
  pso_results <- globpso(objFunc = obj_exact, lower = rep(lower, nPoints), 
                         upper = rep(upper, nPoints), init = mvnorm_pso, 
                         PSO_INFO = psoinfo_exact, loc = parms, verbose = F)

  exact_design <- pso_results$par |> round(4) |> table() |> data.frame()
  colnames(exact_design) <- c("Support", "N")
  exact_design$Support <- as.numeric(levels(exact_design$Support))[exact_design$Support]
  
  exact_val <- pso_results$val
  if (criterion == "D") eff <- (exact_val / approx_val)^(1/4)
  else eff <- (approx_val / exact_val)
  
  exact_res <- list(design = exact_design, efficiency = eff)
  exact_res
}


round_approx <- function(approx){
  approx <- approx |> round(4)
  
  while(sum(approx$Weight) != 1){
    idx <- sample(1:length(approx$Weight), 1)
    if (sum(approx$Weight) > 1) approx$Weight[idx] <- approx$Weight[idx] - 0.0001
    else approx$Weight[idx] <- approx$Weight[idx] + 0.0001
  }
  
  idx0 <- which(approx$Weight == 0)
  if(length(idx0) != 0) approx <- approx[-idx0, ]
  
  if (nrow(approx) != length(unique(approx$Support))){
    cnt <- approx |> count(Support)
    ext <- cnt$Support[which(cnt$n != 1)]
    apprep <- approx$Support == ext
    w <- sum(approx$Weight[apprep])
    approx <- approx[!apprep,]
    approx <- rbind(appro, c(ext, w))
  }
  
  approx
}


# Run the application 
shinyApp(ui = ui, server = server)
