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

hb_doptimal_approx <- function(dw, loc){
  
  # Hunt-Bowman parameters
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  pen = 9e+10
  
  n <- length(dw)/2
  d <- dw[1:n]
  w <- dw[(n+1):(2*n)]
  
  # Evaluate d-optimality criterion value
  if (sum(w) > 1) res <- pen
  else {
    mat_list <- lapply(1:n, function(x) w[x] * hb_mat(d[x], c1, tau, b0, b1))
    inf_mat <- Reduce("+", mat_list)
    res <- -det(inf_mat)
  }
  
  res
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
  n <- length(input)/2
  pen = 9e+10
  
  d <- input[1:n]
  w <- input[(n+1):(2*n)]
  
  mat_list <- lapply(1:n, function(x) w[x] * hb_mat(d[x], c1, tau, b0, b1))
  M <- Reduce("+", mat_list)
  diag(M) <- diag(M) + 1e-10
  b <- matrix(c(0, 1, 0, 0))
  
  if (rcond(M) < 2.220446e-16 || sum(w) > 1){
    res = pen
  } 
  else{
    M_inv <- solve(M)
    res <- t(b) %*% M_inv %*% b
  }
  
  res
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
  n <- length(input)/2
  pen = 9e+10
  
  d <- input[1:n]
  w <- input[(n+1):(2*n)]
  
  mat_list <- lapply(1:n, function(x) w[x] * hb_mat(d[x], c1, tau, b0, b1))
  M <- Reduce("+", mat_list)
  diag(M) <- diag(M) + 1e-10
  h <- matrix(c(-tau, -c1, 0, 0))
  
  if (rcond(M) < 2.220446e-16 || sum(w) > 1){
    res = pen
  } 
  else{
    M_inv <- solve(M)
    res <- t(h) %*% M_inv %*% h
  }
  
  res
}

### exp-log model

# exp-log model parameters
el_params <- function(c0, c1, b0, b1){
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

exp_log_mat <- function(d, c0, c1, b0, b1){
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
  
  M <- exp_log_mat(d, c0, c1, b0, b1)
  
  -det(M)
}

el_doptimal_approx <- function(input, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  n <- length(input)/2
  pen = 9e+10
  
  d <- input[1:n]
  w <- input[(n+1):(2*n)]
  
  mat_list <- lapply(1:n, function(x) w[x] * exp_log_f(d[x], c0, c1, b0, b1) %*% t(exp_log_f(d[x], c0, c1, b0, b1)))
  M <- Reduce("+", mat_list)
  
  if(sum(w) > 1) res = pen
  else res = -det(M)
  
  res
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

el_hoptimal_approx <- function(input, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  n = length(input)/2
  pen = 9e+10
  
  d <- input[1:n]
  w <- input[(n+1):(2*n)]
  
  h1 <- -c1
  h2 <- -c0
  h3 <- b1 * exp(b0) * (1 - exp(b0)) / (exp(b0) + 1)^3
  h4 <- exp(b0) / (exp(b0) + 1)^2 
  h <- matrix(c(h1, h2, h3, h4))
  
  mat_list <- lapply(1:n, function(x) w[x] * exp_log_f(d[x], c0, c1, b0, b1) %*% t(exp_log_f(d[x], c0, c1, b0, b1)))
  M <- Reduce("+", mat_list)
  
  if (rcond(M) < 2.220446e-16 || sum(w) > 1){
    res = pen
  }
  else {
    M_inv <- solve(M)
    res = t(h) %*% M_inv %*% h
  }
  
  res
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
  
  M <- exp_log_mat(d, c0, c1, b0, b1)
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
  
  n <- length(input)/2
  
  d <- input[1:n]
  w <- input[(n+1):(2*n)]
  
  b <- exp_log_b(tau, c0, c1, b0, b1)
  
  mat_list <- lapply(1:n, function(x) w[x] * exp_log_f(d[x], c0, c1, b0, b1) %*% t(exp_log_f(d[x], c0, c1, b0, b1)))
  M <- Reduce("+", mat_list)
  diag(M) <- diag(M) + 1e-10
  
  if (rcond(M) < 2.220446e-16 || sum(w) > 1){
    res = pen
  }
  else {
    M_inv <- solve(M)
    res = t(b) %*% M_inv %*% b
  }
  
  res
}

## Simple logistic model D-optimal
logistic_params <- function(alpha, beta){
  c(alpha, beta)
}

logistic <- function(d, alpha, beta){
  eta <- alpha + beta * d
  exp(eta) / (1 + exp(eta))
}

# Dose-response plot for logistic models
logistic_plot <- function(parms, bound){
  alpha = parms[1]
  beta = parms[2]
  lb = bound * -1
  ub = bound
  
  fp <- seq(lb, ub, by = (ub - lb)/100)
  log_df <- data.frame(dose = fp, response = sapply(fp, function(x) logistic(x, alpha, beta)))
  ggplot(data = log_df, aes(x = dose, y = response)) + 
    geom_line()
}

logit_doptimal <- function(d, loc){
  a <- loc[1]
  b <- loc[2]
  n <- length(d)
  
  theta <- a + b * d
  omega <- exp(theta) / (1 + exp(theta))^2
  
  inf_mat <- matrix(c(sum(omega), sum(d * omega), sum(d * omega), sum(d^2 * omega)), nrow = 2)
  -det(inf_mat) / n^2
}

logit_doptimal_approx <- function(input, loc){
  a <- loc[1]
  b <- loc[2]
  pen = 9e+10
  n <- length(input)/2
  
  d <- input[1:n]
  w <- input[(n+1):(2*n)]
  
  theta <- a + b * d
  omega <- exp(theta) / (1 + exp(theta))^2
  
  if(sum(w) > 1){
    res = pen
  }
  else{
    M <- matrix(c(sum(w * omega), sum(w * omega * d), 
                  sum(w * omega * d), sum(w * omega * d^2)), 
                nrow = 2)
    res = -det(M)
  }
  
  res
}

## Quadratic logistic D-optimal
qlogistic_params <- function(alpha, beta1, beta2){
  c(alpha, beta1, beta2)
}

qlogistic <- function(d, alpha, beta1, beta2){
  eta <- alpha + beta1 * d + beta2 * d^2
  exp(eta) / (1 + exp(eta))
}

# Dose-response plot for quadratic logistic models
qlogistic_plot <- function(parms, bound){
  alpha = parms[1]
  beta1 = parms[2]
  beta2 = parms[3]
  lb = bound * -1
  ub = bound
  
  fp <- seq(lb, ub, by = (ub - lb)/100)
  qlog_df <- data.frame(dose = fp, response = sapply(fp, function(x) qlogistic(x, alpha, beta1, beta2)))
  ggplot(data = qlog_df, aes(x = dose, y = response)) + 
    geom_line()
}

qlogit_doptimal <- function(d, loc){
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
  
  -det(inf_mat) / n^3
}

qlogit_doptimal_approx <- function(input, loc){
  a <- loc[1]
  b1 <- loc[2]
  b2 <- loc[3]
  n <- length(input)/2
  pen <- 9e+10
  
  d <- input[1:n]
  w <- input[(n+1):(2*n)]
  
  theta <- a + b1 * d + b2 * d^2
  omega <- exp(theta) / (1 + exp(theta))^2
  
  if (sum(w) > 1) res = pen
  else {
    M <- matrix(c(sum(w * omega), sum(w * d * omega), sum(w * d^2 * omega), 
                  sum(w * d * omega), sum(w * d^2 * omega), sum(w * d^3 * omega), 
                  sum(w * d^2 * omega), sum(w * d^3 * omega), sum(w * d^4 * omega)), 
                nrow = 3)
    res = -det(M)
  }
  
  res
}

## Cubic logistic models
clogistic_params <- function(alpha, beta1, beta2, beta3){
  c(alpha, beta1, beta2, beta3)
}

clogistic <- function(d, alpha, beta1, beta2, beta3){
  eta <- alpha + beta1 * d + beta2 * d^2 + beta3 * d^3
  exp(eta) / (1 + exp(eta))
}

# Dose-response plot for cubic logistic models
clogistic_plot <- function(parms, bound){
  alpha = parms[1]
  beta1 = parms[2]
  beta2 = parms[3]
  beta3 = parms[4]
  lb = bound * -1
  ub = bound
  
  fp <- seq(lb, ub, by = (ub - lb)/100)
  clog_df <- data.frame(dose = fp, response = sapply(fp, function(x) clogistic(x, alpha, beta1, beta2, beta3)))
  ggplot(data = clog_df, aes(x = dose, y = response)) + 
    geom_line()
}

clogit_doptimal <- function(d, loc){
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

clogit_doptimal_approx <- function(input, loc){
  a <- loc[1]
  b1 <- loc[2]
  b2 <- loc[3]
  b3 <- loc[4]
  pen = 9e+10
  n <- length(input)/2
  
  d <- input[1:n]
  w <- input[(1+n):(2*n)]
  
  theta <- a + b1 * d + b2 * d^2 + b3 * d^3
  omega <- exp(theta) / (1 + exp(theta))^2
  
  if(sum(w) > 1) res = pen
  else{
    M <- matrix(c(sum(omega), sum(d * omega), sum(d^2 * omega), sum(d^3 * omega), 
                  sum(d * omega), sum(d^2 * omega), sum(d^3 * omega),  sum(d^4 * omega), 
                  sum(d^2 * omega), sum(d^3 * omega), sum(d^4 * omega), sum(d^5 * omega), 
                  sum(d^3 * omega), sum(d^4 * omega), sum(d^5 * omega), sum(d^6 * omega)), 
                nrow = 4)
    res = -det(M)
  }
  
  res
}

### PSO function

## Initalize PSO settings

psoinfo_setting <- function(nSwarms = 64, Iters = 1000){
  getPSOInfo(nSwarm = nSwarms, maxIter = Iters, w0 = 0.9, w1 = 0.4, w_var = 1)
}

# PSO for exact design
hormesis_exact <- function(model, criterion, nPoints, parms, psoinfo, upper, lower){
  
  # set the lower and upper bounds for PSO
  lb <- rep(lower, nPoints)
  ub <- rep(upper, nPoints)
  
  if (model == "HuntBowman"){
    if (criterion == "D") obj <- hb_doptimal
    else if (criterion == "tau") obj <- hb_tauoptimal
    else if (criterion == "h") obj <- hb_hoptimal
  } else if (model == "ExpLog"){
    if (criterion == "D") obj <- el_doptimal
    else if (criterion == "tau"){
      obj <- el_tauoptimal
      tau <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, 
                     c0 = parms[1], c1 = parms[2], b0 = parms[3], b1 = parms[4])$root
      parms <- c(parms, tau)
    } 
    else if (criterion == "h") obj <- el_hoptimal
  } 
  else if (model == "logistic") obj <- logit_doptimal
  else if (model == "qlogistic") obj <- qlogit_doptimal
  else if (model == "clogistic") obj <- clogit_doptimal
  
  
  pso_res <- globpso(objFunc = obj, 
                     lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, 
                     loc = parms)
  
  pso_res$exact <- pso_res$par |> round(4) |> table() |> data.frame()
  colnames(pso_res$exact) <- c("Support", "N")
  pso_res
}

hormesis_approx <- function(model, criterion, parms, psoinfo, upper, lower){
  
  if (model == "HuntBowman"){
    nDim = 4
    if (criterion == "D") obj <- hb_doptimal_approx
    else if (criterion == "tau") obj <- hb_tauoptimal_approx
    else if (criterion == "h") obj <- hb_hoptimal_approx
  } 
  else if (model == "ExpLog"){
    nDim = 4
    if (criterion == "D") obj <- el_doptimal_approx
    else if (criterion == "tau"){
      obj <- el_tauoptimal_approx
      tau <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, 
                     c0 = parms[1], c1 = parms[2], b0 = parms[3], b1 = parms[4])$root
      parms <- c(parms, tau)
    } 
    else if (criterion == "h") obj <- el_hoptimal_approx
  } 
  else if (model == "logistic"){
    obj <- logit_doptimal_approx
    nDim = 2
  } 
  else if (model == "qlogistic"){
    obj <- qlogit_doptimal_approx
    nDim = 4
  } 
  else if (model == "clogistic"){
    obj <- clogit_doptimal_approx
    nDim = 5  
  } 
  
  lb <- c(rep(lower, nDim), rep(0, nDim))
  ub <- c(rep(upper, nDim), rep(1, nDim))
  
  pso_res <- globpso(objFunc = obj, 
                     lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, 
                     loc = parms)
  
  pso_res$approx <- data.frame(support = pso_res$par[1:nDim], weight = pso_res$par[(nDim+1):(2*nDim)]) |> 
    arrange(support)
  pso_res$approx$support <- pso_res$approx$support |> round(4)
  pso_res$approx$weight <- pso_res$approx$weight |> round(3)
  
  idx0 <- which(pso_res$approx$weight == 0)
  if(length(idx0) != 0) pso_res$approx <- pso_res$approx[-idx0, ]
  
  
  if (nrow(pso_res$approx) != length(unique(pso_res$approx$support))){
    cnt <- pso_res$approx |> count(support)
    ext <- cnt$support[which(cnt$n != 1)]
    apprep <- pso_res$approx$support == ext
    w <- sum(pso_res$approx$weight[apprep])
    pso_res$approx <- pso_res$approx[!apprep,]
    pso_res$approx <- rbind(pso_res$approx, c(ext, w))
  }
  
  pso_res
}

hormesis_pso <- function(model, criterion, parms, psoinfo, upper, lower, nPoints, nRep){
  exact_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- hormesis_exact(model = model, criterion = criterion, 
                                                                               nPoints = nPoints, parms = parms, 
                                                                               psoinfo = psoinfo, upper = upper, lower = lower))
  val_list <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  best_idx <- which.min(val_list)
  exact_design <- pso_results[[best_idx]]$exact
  exact_val <- pso_results[[best_idx]]$val
  
  psoinfo_approx <- psoinfo_setting(256, 2000)
  approx_result <- hormesis_approx(model = model, criterion = criterion, parms = parms, 
                                   psoinfo = psoinfo_approx, upper = upper, lower = lower)
  approx_design <- approx_result$approx
  approx_val <- approx_result$val
  
  nDim <- length(parms)
  if (criterion == "D") eff <- (exact_val/approx_val) ^ (1 / nDim)
  else eff <- (approx_val / exact_val)
  eff = round(eff, 4)
  
  list(exact_design = exact_design, approximate_design = approx_design, efficiency = eff)
}
