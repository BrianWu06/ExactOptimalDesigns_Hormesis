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

hb_pso <- function(criterion, nPoints, parms, psoinfo, upper){
  
  # set the lower and upper bounds for PSO
  lb <- rep(0, nPoints)
  ub <- rep(upper, nPoints)
  
  if (criterion == "D") obj <- hb_doptimal
  else if (criterion == "tau") obj <- hb_tauoptimal
  else if (criterion == "h") obj <- hb_hoptimal
  
  pso_res <- globpso(objFunc = obj, 
                     lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, 
                     loc = parms)
  
  
  pso_res$val <- -1 * pso_res$val
  pso_res$par <- pso_res$par |> sort() |> round(4)
  pso_res$history <- pso_res$history * -1
  pso_res
}

hb_approx_pso

hb_doptimal_pso_rep <- function(criterion, nRep, nPoints = 4, parms, psoinfo, upper){
  
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