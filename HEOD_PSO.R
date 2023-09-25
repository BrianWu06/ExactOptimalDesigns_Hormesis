library(globpso)
library(dplyr)

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
logistic_doptimal_pso_rep <- function(nRep, nPoints = 2, parms, psoinfo){
  
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
