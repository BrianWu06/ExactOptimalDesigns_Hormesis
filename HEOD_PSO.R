library(globpso)

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
  
  # Number of experiment points
  d <- c(0, d)
  w <- length(d)
  
  # Evaluate d-optimality criterion value
  mat_list <- lapply(d, function(x) 1/w * hb_mat(x, c1, tau, b0, b1))
  inf_mat <- Reduce("+", mat_list)
  -det(inf_mat)
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
  
  pso_res$val <- -1 * pso_res$val
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
  
  # Criterion Value of each replication.
  hb_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val) 
  # Design points found in each replication.
  hb_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par) 
  # Best result of each iteration in each replication. (For finding the appropriate nRep.)
  hb_list$history <- sapply(1:nRep, function(x) pso_results[[x]]$history) 
  
  best_idx <- which.max(hb_list$val)
  # Best efficiency among all replications.
  hb_list$result$best_val <- hb_list$val[best_idx]
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

# Find the exact design points that minimizes the h-optimality efficiency
# of exp-log model by PSO
exp_log_pso <- function(nPoints, parms, psoinfo){
  
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

# Replicate m PSO results of exp-log model h-optimality
exp_log_h_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo){
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- exp_log_pso(nPoints, parms, psoinfo))
  
  exp_log_list <- list(efficiency = c(), design_points = c(), history = c(), 
                       result = list(best_eff = 0, best_points = c()))
  
  # Criterion value of each replication.
  exp_log_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  # Design points found in each replication.
  exp_log_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par)
  # Best result of each iteration in each replication. (For finding the appropriate nRep.)
  exp_log_list$history <- sapply(1:nRep, function(x) pso_results[[x]]$history)
  
  best_idx <- which.min(exp_log_list$val)
  # Best efficiency among all replications.
  exp_log_list$result$best_val <- exp_log_list$val[best_idx]
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
  res
}

exp_log_tau_pso <- function(nPoints, parms, psoinfo){
  
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
  pso_res$history <- pso_res$history
  pso_res
}

exp_log_tau_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo){
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- exp_log_tau_pso(nPoints, parms, psoinfo))
  
  exp_log_list <- list(design_points = c(), weight = c(), val = c())
  
  # Efficiency of each replication.
  exp_log_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  # Design points found in each replication.
  exp_log_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par)
  # Best result of each iteration in each replication. (For finding the appropriate nRep.)
  exp_log_list$history <- sapply(1:nRep, function(x) pso_results[[x]]$history)
  
  best_idx <- which.min(exp_log_list$val)
  # Best efficiency among all replications.
  exp_log_list$result$best_val <- exp_log_list$val[best_idx]
  # Design points correspond to the best efficiency.
  exp_log_list$result$best_points <- exp_log_list$design_points[,best_idx]
  
  exp_log_list
}

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

logistic_d_pso <- function(npoints, parms, psoinfo){
  lb <- rep(-10, npoints)
  ub <- rep(10, npoints)
  
  pso_res <- globpso(objFunc = logistic_doptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$val <- -1 * pso_res$val
  pso_res
}

logistic_d_pso_rep <- function(nRep, nPoints = 2, parms, psoinfo){
  
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- logistic_d_pso(nPoints, parms, psoinfo))
  
  logistic_list <- list(design_points = c(), val = c())
  logistic_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par)
  logistic_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  logistic_list$history <- sapply(1:nRep, function(x) pso_results[[x]]$history)
  
  best_idx <- which.max(logistic_list$val)
  logistic_list$result$best_val <- logistic_list$val[best_idx]
  logistic_list$result$best_points <- logistic_list$design_points[, best_idx] |> sort() |> round(4)
  
  logistic_list
}

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

qlogistic_pso <- function(npoints, parms, psoinfo){
  lb <- rep(-10, npoints)
  ub <- rep(10, npoints)
  
  pso_res <- globpso(objFunc = qlogistic_doptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$val <- -1 * pso_res$val
  pso_res
}

qlogistic_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo){
  
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- qlogistic_pso(nPoints, parms, psoinfo))
  
  qlogistic_list <- list(design_points = c(), val = c())
  qlogistic_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par)
  qlogistic_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  qlogistic_list$history <- sapply(1:nRep, function(x) pso_results[[x]]$history)
  
  best_idx <- which.max(qlogistic_list$val)
  qlogistic_list$result$best_val <- qlogistic_list$val[best_idx]
  qlogistic_list$result$best_points <- qlogistic_list$design_points[, best_idx] |> sort() |> round(4)
  
  qlogistic_list
}