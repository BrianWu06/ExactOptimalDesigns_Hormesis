## Initalize PSO settings

psoinfo_setting <- function(nSwarms = 64, Iters = 1000){
  getPSOInfo(nSwarm = nSwarms, maxIter = Iters, w0 = 0.9, w1 = 0.4, w_var = 1)
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

## Hunt-Bowman D-optimal

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

# Find the D-optimal exact design for Hunt-Bowman model
hb_doptimal_pso <- function(nPoints, parms, psoinfo, upper){
  
  # set the lower and upper bounds for PSO
  lb <- rep(0, nPoints)
  ub <- rep(upper, nPoints)
  
  #run PSO
  pso_res <- globpso(objFunc = hb_doptimal, 
                     lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, 
                     loc = parms)
  
  pso_res$val <- -1 * pso_res$val
  pso_res$par <- pso_res$par |> sort() |> round(4)
  pso_res$history <- pso_res$history * -1
  pso_res
}

# Find D-optimal approximate design for the Hunt-Bowman model
hb_doptimal_approx <- function(parms, ub){
  psoinfo <- psoinfo_setting(Iters = 400)
  approx_design <- hb_doptimal_pso(nPoints = 4, parms, psoinfo, ub)
  approx_design
}

# Replicate m PSO results of the D-optimal exact design for Hunt-Bowman model
hb_doptimal_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo, upper){
  
  # The D-optimal approximate design for the Hunt-Bowman model under certain parameter set.
  approx_design <- hb_doptimal_approx(parms, upper)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- hb_doptimal_pso(nPoints, parms, psoinfo, upper))
  
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

hb_tauoptimal_approx <- function(input, loc){
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  n <- length(input)/2
  pen = 999999999
  
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

# Find tau-optimal exact design for the Hunt-Bowman model
hb_tauoptimal_pso <- function(nPoints, parms, psoinfo, upper){
  lb <- rep(0, nPoints)
  ub <- rep(upper, nPoints)
  
  pso_res <- globpso(objFunc = hb_tauoptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$par <- pso_res$par |> round(4) |> sort()
  pso_res
}

# Find tau-optimal approximate design for the Hunt-Bowman model
hb_tauoptimal_approx <- function(parms, upper){
  psoinfo <- psoinfo_setting()
  approx_design <- hb_tauoptimal_pso(nPoints = 2, parms, psoinfo, upper)
  approx_design
}

# Replicate m PSO results of the tau-optimal exact design for Hunt-Bowman model
hb_tauoptimal_pso_rep <- function(nRep, nPoints = 2, parms, psoinfo, upper){
  
  # The tau-optimal approximate design for the Hunt-Bowman model under certain parameter set.
  approx_design <- hb_tauoptimal_approx(parms, upper)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- hb_tauoptimal_pso(nPoints, parms, psoinfo, upper))
  
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
  l <- length(input)
  n <- l/2
  d <- input[1:n]
  w <- input[(n+1):(2*n)]
  
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  pen = 999999999
  
  mat_list <- lapply(1:n, function(x) w[x] * hb_mat(d[x], c1, tau, b0, b1))
  M <- Reduce("+", mat_list)
  diag(M) = diag(M) + 1e-10
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

# Find h-optimal exact design for the Hunt-Bowman model
hb_hoptimal_pso <- function(nPoints, parms, psoinfo, upper){
  lb <- rep(0, nPoints)
  ub <- rep(upper, nPoints)
  
  pso_res <- globpso(objFunc = hb_hoptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$par <- pso_res$par |> round(4) |> sort()
  pso_res
}

# Find h-optimal approximate design for the Hunt-Bowman model
hb_hoptimal_approx_pso <- function(nPoints, parms, upper){
  psoinfo <- psoinfo_setting(nSwarms =512, Iters = 5000)
  lb <- rep(0, 2 * nPoints)
  ub <- c(rep(upper, nPoints), rep(1, nPoints))
  
  pso_res <- globpso(objFunc = hb_hoptimal_approx, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$par <- pso_res$par |> round(4)
  pso_res
}

# Replicate m PSO results of the h-optimal exact design for Hunt-Bowman model
hb_hoptimal_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo, upper){
  
  # The h-optimal approximate design for the Hunt-Bowman model under certain parameter set.
  approx_design <- hb_hoptimal_approx_pso(nPoints = 4, parms = parms, upper)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- hb_hoptimal_pso(nPoints, parms, psoinfo, upper))
  
  hb_list <- list(val = c(), design_points = c(), history = c(), 
                  result = list(best_val = 0),
                  approximate_design = data.frame(support_points = approx_design$par[1:4],
                                                    weight = approx_design$par[5:8]))
  hb_list$approximate_design <- hb_list$approximate_design %>% arrange(support_points)
  idx0 <- print(match(0, hb_list$approximate_design$weight, -1))
  if (idx0 != -1) hb_list$approximate_design <- hb_list$approximate_design[-idx0, ]

  if (nrow(hb_list$approximate_design) != length(unique(hb_list$approximate_design$support_points))){
    cnt <- hb_list$approximate_design |> count(support_points)
    ext <- cnt$support_points[which(cnt$n != 1)]
    apprep <- hb_list$approximate_design$support_points == ext
    w <- sum(hb_list$approximate_design$weight[apprep])
    hb_list$approximate_design <- hb_list$approximate_design[!apprep,]
    hb_list$approximate_design <- rbind(hb_list$approximate_design, c(ext, w))
   }
  
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
exp_log_doptimal <- function(d, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  
  M <- exp_log_mat(d, c0, c1, b0, b1)
  
  -det(M)
}

# Find D-optimal exact design for the exp-log model
exp_log_doptimal_pso <- function(nPoints, parms, psoinfo, upper){
  lb <- rep(0, nPoints)
  ub <- rep(upper, nPoints)
  
  pso_res <- globpso(objFunc = exp_log_doptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$val <- pso_res$val * -1
  pso_res$par <- pso_res$par |> round(4) |> sort()
  pso_res
}

# Find D-optimal approximate design for the exp-log model
exp_log_doptimal_approx <- function(parms, upper){
  psoinfo <- psoinfo_setting()
  approx_design <- exp_log_doptimal_pso(nPoints = 4, parms, psoinfo, upper)
  approx_design
}

# Replicate m PSO results of the D-optimal exact design for exp-log model
exp_log_doptimal_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo, upper){
  
  # The D-optimal approximate design for the exp-log model under certain parameter set.
  approx_design <- exp_log_doptimal_approx(parms, upper)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- exp_log_doptimal_pso(nPoints, parms, psoinfo, upper))
  
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
exp_log_h <- function(d, c0, c1, b0, b1){
  h1 <- -c1
  h2 <- -c0
  h3 <- b1 * exp(b0) * (1 - exp(b0)) / (exp(b0) + 1)^3
  h4 <- exp(b0) / (exp(b0) + 1)^2 
  h <- matrix(c(h1, h2, h3, h4))
  
  h
}

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
  
  mat_list <- lapply(1:n, function(x) w[x] * exp_log_mat(d[x], c0, c1, b0, b1))
  M <- Reduce("+", mat_list)
  #print(M)
  
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
exp_log_hoptimal_pso <- function(nPoints, parms, psoinfo, upper){
  
  lb <- rep(0, nPoints)
  ub <- rep(upper, nPoints)
  
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
exp_log_hoptimal_approx_pso <- function(nPoints, parms, upper){
  psoinfo <- psoinfo_setting(nSwarms = 256, Iters = 2000)
  lb <- rep(0, 2 * nPoints - 1)
  ub <- c(rep(upper, nPoints), rep(1, nPoints - 1))
  
  pso_res <- globpso(objFunc = exp_log_hoptimal_approx, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$par <- c(pso_res$par, 1 - sum(pso_res$par[(nPoints+1):(2*nPoints-1)]))
  pso_res$par <- pso_res$par |> round(4)
  pso_res
}

# Replicate m PSO results of the h-optimal exact design for exp-log model
exp_log_hoptimal_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo, upper){
  
  # The h-optimal approximate design for the exp-log model under certain parameter set.
  approx_design <- exp_log_hoptimal_approx_pso(4, parms, upper)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- exp_log_hoptimal_pso(nPoints, parms, psoinfo, upper))
  
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

# Find tau-optimal exact design for the exp-log model
exp_log_tauoptimal_pso <- function(nPoints, parms, psoinfo, upper){
  
  lb <- c(rep(0, nPoints))
  ub <- c(rep(upper, nPoints))
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
exp_log_tauoptimal_approx <- function(parms, upper){
  psoinfo <- psoinfo_setting()
  approx_design <- exp_log_tauoptimal_pso(nPoints = 2, parms, psoinfo, upper)
  approx_design
}

# Replicate m PSO results of the tau-optimal exact design for exp-log model
exp_log_tauoptimal_pso_rep <- function(nRep, nPoints = 2, parms, psoinfo, upper){
  
  # The tau-optimal approximate design for the exp-log model under certain parameter set.
  approx_design <- exp_log_tauoptimal_approx(parms, upper)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- exp_log_tauoptimal_pso(nPoints, parms, psoinfo, upper))
  
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
logistic_doptimal_pso <- function(nPoints, parms, psoinfo, bd){
  lb <- rep(-1 * bd, nPoints)
  ub <- rep(bd, nPoints)
  
  pso_res <- globpso(objFunc = logistic_doptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$val <- -1 * pso_res$val
  pso_res$par <- pso_res$par |> round(4) |> sort()
  pso_res
}

# Find D-optimal approximate design for the simple logistic model
logistic_doptimal_approx <- function(parms, bound){
  psoinfo <- psoinfo_setting()
  approx_design <- logistic_doptimal_pso(nPoints = 2, parms, psoinfo, bound)
  approx_design
}

# Replicate m PSO results of the D-optimal exact design for simple logistic model
logistic_pso_rep <- function(nRep, nPoints = 2, parms, psoinfo, bound){
  
  # The D-optimal approximate design for the simple logistic model under certain parameter set.
  approx_design <- logistic_doptimal_approx(parms, bound)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- logistic_doptimal_pso(nPoints, parms, psoinfo, bound))
  
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
qlogistic_pso <- function(nPoints, parms, psoinfo, bd){
  lb <- rep(-1 * bd, nPoints)
  ub <- rep(bd, nPoints)
  
  pso_res <- globpso(objFunc = qlogistic_doptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$val <- -1 * pso_res$val
  pso_res$par <- pso_res$par |> round(4) |> sort()
  pso_res
}

# Find D-optimal approximate design for the quadratic logistic model
qlogistic_approx_pso <- function(parms, bd){
  lb <- c(rep(-1 * bd, 4), rep(0, 3))
  ub <- c(rep(bd, 4), rep(1, 3))
  psoinfo <- psoinfo_setting()
  
  pso_res <- globpso(objFunc = qlogistic_doptimal_approx, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$val <- -1 * pso_res$val
  pso_res$par <- c(pso_res$par, 1 - sum(pso_res$par[5:7])) |> round(4)
  pso_res
}

# Replicate m PSO results of the D-optimal exact design for quadratic logistic model
qlogistic_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo, bound){
  
  # The D-optimal approximate design for the quadratic logistic model under certain parameter set.
  approx_design <- qlogistic_approx_pso(parms, bound)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- qlogistic_pso(nPoints, parms, psoinfo, bound))
  
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
  
  if (nrow(qlogistic_list$approximate_design) != length(unique(qlogistic_list$approximate_design$support_points))){
    cnt <- qlogistic_list$approximate_design |> count(support_points)
    ext <- cnt$support_points[which(cnt$n != 1)]
    apprep <- qlogistic_list$approximate_design$support_points == ext
    w <- sum(qlogistic_list$approximate_design$weight[apprep])
    qlogistic_list$approximate_design <- qlogistic_list$approximate_design[!apprep,]
    qlogistic_list$approximate_design <- rbind(qlogistic_list$approximate_design, c(ext, w))
  }
  
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
clogistic_pso <- function(npoints, parms, psoinfo, bd){
  lb <- rep(-1 * bd, npoints)
  ub <- rep(bd, npoints)
  
  pso_res <- globpso(objFunc = clogistic_doptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$val <- pso_res$val * -1
  pso_res$par <- pso_res$par |> round(4) |> sort()
  pso_res
}

# Find D-optimal approximate design for the cubic logistic model
clogistic_approx_pso <- function(parms, bd){
  lb <- c(rep(-1 * bd, 5), rep(0, 4))
  ub <- c(rep(bd, 5), rep(1, 4))
  psoinfo <- psoinfo_setting()
  
  pso_res <- globpso(objFunc = clogistic_doptimal_approx, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$val <- pso_res$val * -1
  pso_res$par <- c(pso_res$par, 1 - sum(pso_res$par[6:9])) |> round(4)
  pso_res
}

# Replicate m PSO results of the D-optimal exact design for cubic logistic model
clogistic_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo, bound){
  
  # The D-optimal approximate design for the cubic logistic model under certain parameter set.
  approx_design <- clogistic_approx_pso(parms, bound)
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- clogistic_pso(nPoints, parms, psoinfo, bound))
  
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
  
  if (nrow(clogistic_list$approximate_design) != length(unique(clogistic_list$approximate_design$support_points))){
    cnt <- clogistic_list$approximate_design |> count(support_points)
    ext <- cnt$support_points[which(cnt$n != 1)]
    apprep <- clogistic_list$approximate_design$support_points == ext
    w <- sum(clogistic_list$approximate_design$weight[apprep])
    clogistic_list$approximate_design <- clogistic_list$approximate_design[!apprep,]
    clogistic_list$approximate_design <- rbind(clogistic_list$approximate_design, c(ext, w))
  }
  
  best_idx <- which.max(clogistic_list$val)
  # Best criterion value among all replications.
  clogistic_list$result$best_val <- clogistic_list$val[best_idx]
  # Design points correspond to the best criterion value.
  clogistic_list$result$design_points <- clogistic_list$design_points[, best_idx] |> sort() |> round(4)
  # Efficiency correspond to the best criterion value.
  clogistic_list$result$efficiency <- (clogistic_list$val[best_idx] / approx_design$val)^(1/nPoints)
  
  clogistic_list
}