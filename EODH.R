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
  
  n <- (length(dw)+1)/2
  d <- dw[1:n]
  w <- dw[(n+1):(2*n-1)]
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

## Simple logistic model D-optimal
logistic_parms <- function(alpha, beta){
  c(alpha, beta)
}

logistic <- function(d, alpha, beta){
  eta <- alpha + beta * d
  exp(eta) / (1 + exp(eta))
}

# Dose-response plot for logistic models
logistic_plot <- function(parms, lb, ub){
  alpha = parms[1]
  beta = parms[2]
  
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
  n <- (length(input)+1)/2
  
  d <- input[1:n]
  w <- input[(n+1):(2*n-1)]
  w[n] <- 1 - sum(w)
  
  theta <- a + b * d
  omega <- exp(theta) / (1 + exp(theta))^2
  
  if(w[n] < 0){
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
qlogistic_parms <- function(alpha, beta1, beta2){
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
  n <- (length(input)+1)/2
  pen <- 9e+10
  
  d <- input[1:n]
  w <- input[(n+1):(2*n-1)]
  w[n] <- 1 - sum(w)
  
  theta <- a + b1 * d + b2 * d^2
  omega <- exp(theta) / (1 + exp(theta))^2
  
  if (w[n] < 0) res = pen
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
clogistic_parms <- function(alpha, beta1, beta2, beta3){
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
  n <- (length(input)+1)/2
  
  d <- input[1:n]
  w <- input[(n+1):(2*n-1)]
  w[n] <- 1 - sum(w)
  
  theta <- a + b1 * d + b2 * d^2 + b3 * d^3
  omega <- exp(theta) / (1 + exp(theta))^2
  
  if(w[n] < 0) res = pen
  else{
    M <- matrix(c(sum(w * omega), sum(w * d * omega), sum(w * d^2 * omega), sum(w * d^3 * omega), 
                  sum(w * d * omega), sum(w * d^2 * omega), sum(w * d^3 * omega),  sum(w * d^4 * omega), 
                  sum(w * d^2 * omega), sum(w * d^3 * omega), sum(w * d^4 * omega), sum(w * d^5 * omega), 
                  sum(w * d^3 * omega), sum(w * d^4 * omega), sum(w * d^5 * omega), sum(w * d^6 * omega)), 
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
  else if (model == "logistic"){
    obj_approx <- logit_doptimal_approx
    obj_exact <- logit_doptimal
    nDim = 2
  } 
  else if (model == "qlogistic"){
    obj_approx <- qlogit_doptimal_approx
    obj_exact <- qlogit_doptimal
    nDim = 4
  } 
  else if (model == "clogistic"){
    obj_approx <- clogit_doptimal_approx
    obj_exact <- clogit_doptimal
    nDim = 5  
  } 
   
  #psoinfo_approx <- psoinfo_setting(256, 2000)
  approx_lb <- c(rep(lower, nDim), rep(0, nDim-1))
  approx_ub <- c(rep(upper, nDim), rep(1, nDim-1))
  approx_result <- globpso(objFunc = obj_approx, lower = approx_lb, upper = approx_ub, PSO_INFO = psoinfo_approx, 
                           loc = parms, verbose = T)
  #print(approx_result$par)
  approx_d <- approx_result$par[1:nDim] |> round(4)
  approx_w <- approx_result$par[(nDim+1):(2*nDim-1)] |> round(3)
  approx_w <- c(approx_w, 1 - sum(approx_w))
  approx_design <- data.frame(support = approx_d, weight = approx_w)
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
  
  effr <- efficient.rounding(approx_design$weight, nPoints)
  mu_pso <- lapply(1:length(effr), function(x) rep(approx_design$support[x], effr[x])) |> unlist()
  mvnorm_sd = (upper - lower)/4
  nswarm <- psoinfo_exact$nSwarm/2
  mvnorm_pso <- mvrnorm(n = (nswarm-1), mu = mu_pso, Sigma = diag(nPoints) * mvnorm_sd)
  mvnorm_pso[mvnorm_pso < 0] <- 0
  mvnorm_pso[mvnorm_pso > 0.15] <- 0.15
  mvnorm_pso <- rbind(mvnorm_pso, mu_pso)
  
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- globpso(objFunc = obj_exact, lower = rep(lower, nPoints), 
                                                                        upper = rep(upper, nPoints), init = mvnorm_pso, 
                                                                        PSO_INFO = psoinfo_exact, loc = parms, verbose = T))
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
  
  nEff <- length(parms)
  if (criterion == "D") eff <- (exact_val/approx_val) ^ (1 / nEff)
  else eff <- (approx_val / exact_val)
  eff = round(eff, 4)
  end_time <- Sys.time()
  
  result <- list(approx_design = approx_design, exact_design = exact_design, 
                 efficiency = eff, runtime = end_time - start_time)
  result
}

hormesis_de <- function(model, criterion, parms, upper, lower, nPoints, nRep = 1, deinfo_exact, deinfo_approx){
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
  else if (model == "logistic"){
    obj_approx <- logit_doptimal_approx
    obj_exact <- logit_doptimal
    nDim = 2
  } 
  else if (model == "qlogistic"){
    obj_approx <- qlogit_doptimal_approx
    obj_exact <- qlogit_doptimal
    nDim = 4
  } 
  else if (model == "clogistic"){
    obj_approx <- clogit_doptimal_approx
    obj_exact <- clogit_doptimal
    nDim = 5  
  } 
  
  #psoinfo_approx <- psoinfo_setting(256, 2000)
  approx_lb <- c(rep(lower, nDim), rep(0, nDim-1))
  approx_ub <- c(rep(upper, nDim), rep(1, nDim-1))
  approx_result <- diffevo(objFunc = obj_approx, lower = approx_lb, upper = approx_ub, DE_INFO = deinfo_approx, 
                           loc = parms, verbose = T)
  #print(approx_result$par)
  approx_d <- approx_result$par[1:nDim] |> round(4)
  approx_w <- approx_result$par[(nDim+1):(2*nDim-1)] |> round(3)
  approx_w <- c(approx_w, 1 - sum(approx_w))
  approx_design <- data.frame(support = approx_d, weight = approx_w)
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
  
  effr <- efficient.rounding(approx_design$weight, nPoints)
  mu_de <- lapply(1:length(effr), function(x) rep(approx_design$support[x], effr[x])) |> unlist()
  mvnorm_sd = (upper - lower)/4
  nswarm <- deinfo_exact$nPop/2
  mvnorm_de <- mvrnorm(n = (nswarm-1), mu = mu_de, Sigma = diag(nPoints) * mvnorm_sd)
  mvnorm_de[mvnorm_de < lower] <- lower
  mvnorm_de[mvnorm_de > upper] <- upper
  mvnorm_de <- rbind(mvnorm_de, mu_de)
  
  de_results <- list()
  de_results <- lapply(1:nRep, function(x) de_results[[x]] <- diffevo(objFunc = obj_exact, lower = rep(lower, nPoints), 
                                                                        upper = rep(upper, nPoints), init = mvnorm_de, 
                                                                        DE_INFO = deinfo_exact, loc = parms, verbose = T))
  val_list <- sapply(1:nRep, function(x) de_results[[x]]$val)
  best_idx <- which.min(val_list)
  exact_design <- de_results[[best_idx]]$par |> round(4) |> table() |> data.frame()
  colnames(exact_design) <- c("Support", "N")
  exact_val <- de_results[[best_idx]]$val
  
  # pso_mvnorm <- globpso(objFunc = obj_exact, lower = rep(lower, nPoints), upper = rep(upper, nPoints), 
  #                       init = mvnorm_pso, PSO_INFO = psoinfo, loc = parms, verbose = F)
  # exact_design <- pso_mvnorm$par |> round(4) |> table() |> data.frame()
  # colnames(exact_design) <- c("Support", "N")
  # exact_val <- pso_mvnorm$val
  
  nEff <- length(parms)
  if (criterion == "D") eff <- (exact_val/approx_val) ^ (1 / nEff)
  else eff <- (approx_val / exact_val)
  eff = round(eff, 4)
  end_time <- Sys.time()
  
  result <- list(approx_design = approx_design, exact_design = exact_design, 
                 efficiency = eff, runtime = end_time - start_time)
  result
}
