library(globpso)
source("HEOD_PSO.R")

#Initialize PSO settings
psoinfo <- psoinfo_setting()


# Hunt-Bowman d-optimal
hb1_parms <- hb_parms(c1 = 170, tau = 0.04, b0 = 1.46, b1 = 40) # Set Hunt-Bowman model parameters.
pso_test <- hb_pso(nPoints = 4, parms = hb1_parms, psoinfo = psoinfo) # Set PSO settings
pso_test$par 
pso_test$val 

hb_doptimal(d = c(0, 0, 0.02, 0.02, 0.02, 0.04, 0.04, 0.04, 0.0991, 0.0991, 0.0991, 
                  0, 0.02, 0.04, 0.0991, 0, 0.02, 0.04, 0.0991), 
            loc = hb1_parms)

hb_res <- hb_pso_rep(mRep = 10, nPoints = 5, parms = hb1_parms, psoinfo = psoinfo)

hb_res$efficiency
hb_res$best_point

hb_results <- list()
psoinfo <- psoinfo_setting(nSwarms = 256, Iters = 10000)

for(i in 16:20){
  hb_results[[i]] <- hb_pso_rep(nRep = 10, nPoints = i, 
                                          parms = hb1_parms, psoinfo = psoinfo)
  print(i)
  print(hb_results[[i]]$result$best_eff)
  print(hb_results[[i]]$result$best_points)
}

saveRDS(hb_results, "doptimal.rds")
hb_results <- readRDS("doptimal.rds")
dop <- sapply(4:20, function(x) hb_results[[x]]$result$best_eff^(1/x))
plot(x = c(4:20), y = dop, type = "l", xlab = "N (Number of Observations)", ylab = "Efficiency")
points(x = c(4:20), y = dop, pch = 19)


# exp-log h-optimal
el_par <- exp_log_params(0.15, 89, 3.2, 41)
1 / exp_log_hoptimal(d = c(0.0124, 0.0620, 0.1243), el_par)


h_eff_base = 243738.4
elh_10_1 <- exp_log_hoptimal(d = c(rep(0, 4), rep(0.0108, 5), rep(0.0526, 1), rep(0.1187, 0)), 
                             el_par)
elh_10_2 <- exp_log_hoptimal(d = c(rep(0, 4), rep(0.0108, 3), rep(0.0526, 2), rep(0.1187, 1)), 
                             el_par)
elh_10_3 <- exp_log_hoptimal(d = c(rep(0, 2), rep(0.0108, 5), rep(0.0526, 2), rep(0.1187, 1)), 
                             el_par)
elh_10_4 <- exp_log_hoptimal(d = c(rep(0, 3), rep(0.0108, 5), rep(0.0526, 1), rep(0.1187, 1)), 
                             el_par)
elh_10_5 <- exp_log_hoptimal(d = c(rep(0, 4), rep(0.0108, 4), rep(0.0526, 1), rep(0.1187, 1)), 
                             el_par)
elh_10 <- exp_log_hoptimal(d = c(rep(0, 3), rep(0.0116, 5), rep(0.0561, 1), rep(0.1242, 1)), 
                             el_par)
elh_4_1 <- exp_log_hoptimal(d = c(rep(0, 1), rep(0.0108, 2), rep(0.0526, 1), rep(0.1187, 0)), 
                             el_par)
elh_5 <- exp_log_hoptimal(d = c(rep(0, 1), rep(0.0108, 1), rep(0.0526, 1), rep(0.1187, 2)), 
                          el_par)

5 * c(0.371, 0.501, 0.087, 0.041) 

(h_eff_base/elh_5)^(1/5)

exp_log_1$val
exp_log_1$par
plot(x = c(1:1001), y = exp_log_1$history, type = "l", 
     xlab = "iterations", ylab = "criterion value", 
     main = "Optimal value from each iteration.")

exp_log_6 <- exp_log_pso_rep(nRep = 3, nPoints = 6, parms = el_par, psoinfo = psoinfo)
exp_log_6$efficiency
exp_log_6$design_points
exp_log_6$result$best_eff
exp_log_6$result$best_points


exp_log_results <- list()
psoinfo <- psoinfo_setting(nSwarms = 512, Iters = 30000)

for(i in 15:17){
  exp_log_results[[i]] <- exp_log_pso_rep(nRep = 10, nPoints = i, 
                                          parms = el_par, psoinfo = psoinfo)
  print(i)
}


psoinfo <- psoinfo_setting(nSwarms = 128, Iters = 2000)

for(i in 4:10){
  exp_log_results[[i]] <- exp_log_pso_rep(nRep = 1, nPoints = i, 
                                          parms = el_par, psoinfo = psoinfo)
  
  print(i)
  print(exp_log_results[[i]]$result$best_eff)
  print(exp_log_results[[i]]$result$best_points)
}

psoinfo <- psoinfo_setting(nSwarms = 128, Iters = 30000)

for(i in 15:17){
  exp_log_results[[i]] <- exp_log_pso_rep(nRep = 10, nPoints = i, 
                                          parms = el_par, psoinfo = psoinfo)
  print(i)
}

exp_log_results[[15]]$efficiency
exp_log_results[[15]]$design_points
exp_log_results[[17]]$result$best_eff
exp_log_results[[17]]$result$best_points

saveRDS(exp_log_results, "hoptimal.rds")


ndesigns <- matrix(ncol = 4, nrow = 17)
rownames(ndesigns) <- c(4:20)
colnames(ndesigns) <- c(0.000, 0.0108, 0.0526, 0.1187)
for(j in 4:20){
  ndesigns[j-3, 1] = .371 * j
  ndesigns[j-3, 2] = .501 * j
  ndesigns[j-3, 3] = .087 * j
  ndesigns[j-3, 4] = .041 * j
}


h_round <- list(points = matrix(ncol = 4, nrow = 17), eff = c(0,0,0))
rownames(h_round$points) <- c(4:20)
colnames(h_round$points) <- c(0.000, 0.0108, 0.0526, 0.1187)
h_round$eff[4] <- 1 / exp_log_hoptimal(d = c(rep(0, 1), rep(0.0108, 2), rep(0.0526, 1), rep(0.1187, 0)), 
                 loc = el_par)
h_round$points[4-3,] <- c(1, 2, 1, 0)
saveRDS(h_round, "hround.rds")

ndesigns[5-3,]
1 / exp_log_hoptimal(d = c(rep(0, 1), rep(0.0108, 2), rep(0.0526, 1), rep(0.1187, 0)), 
                     loc = el_par)


h_best <- list(points = matrix(ncol = 4, nrow = 17), eff = c(0,0,0))
rownames(h_best$points) <- c(4:20)
colnames(h_best$points) <- c(0.000, 0.0108, 0.0526, 0.1187)
for(i in 4:20){
  res = 0
  points = c(0,0,0,0)
  for(a in 0:i){
    for(b in 0:(i-a)){
      for(c in 0:(i-a-b)){
        eff <- 1 / exp_log_hoptimal(d = c(rep(0, a), rep(0.0108, b), rep(0.0526, c), rep(0.1187, i-a-b-c)), 
                                    loc = el_par)
        if(eff > res){
          res = eff
          points <- c(a, b, c, i-a-b-c)
        }
      }
    }
  }
  h_best$points[i-3,] <- points
  h_best$eff[i] <- res
}
h_best$points
h_best$eff
saveRDS(h_best, "hbest.rds")

hbest2 <- data.frame("0" = c(0), "0.013" = c(0), "0.0591" = c(0), "0.1242" = c(0), val = c(0))
for (i in 0:5){
  for(j in 0:(5-i)){
    for (k in 0:(5-i-j)){
      hbest2 <- rbind(hbest2, c(i, j, k, 5-i-j-k, 
                                (elhapprox$val / exp_log_hoptimal(d = c(rep(0, i), rep(0.013, j), rep(0.0591, k), rep(0.1242, 5-i-j-k)), loc = el_par))^(1/5)))
    }
  }
}
hbest2
head(hbest2[order(hbest2$val, decreasing = T), ])

hbest3 <- data.frame("0" = c(0), "0.0116" = c(0), "0.0561" = c(0), "0.1242" = c(0), val = c(0))
for (i in 0:10){
  for(j in 0:(10-i)){
    for (k in 0:(10-i-j)){
      hbest3 <- rbind(hbest3, c(i, j, k, 10-i-j-k, 
                                (elhapprox$val / exp_log_hoptimal(d = c(rep(0, i), rep(0.0116, j), rep(0.0561, k), rep(0.1242, 10-i-j-k)), loc = el_par))^(1/10)))
    }
  }
}
head(hbest3[order(hbest3$val, decreasing = T), ])

  hbest4 <- data.frame("0" = c(0), "0.0106" = c(0), "0.0553" = c(0), "0.1242" = c(0), val = c(0))
for (i in 0:15){
  for(j in 0:(15-i)){
    for (k in 0:(15-i-j)){
      hbest4 <- rbind(hbest4, c(i, j, k, 15-i-j-k, 
                                (elhapprox$val / exp_log_hoptimal(d = c(rep(0, i), rep(0.0106, j), rep(0.0553, k), rep(0.1242, 15-i-j-k)), loc = el_par))^(1/15)))
    }
  }
}
head(hbest4[order(hbest4$val, decreasing = T), ])

hbest5 <- data.frame("0" = c(0), "0.0113" = c(0), "0.0532" = c(0), "0.1191" = c(0), val = c(0))
for (i in 0:20){
  for(j in 0:(20-i)){
    for (k in 0:(20-i-j)){
      hbest5 <- rbind(hbest5, c(i, j, k, 20-i-j-k, 
                                (elhapprox$val / exp_log_hoptimal(d = c(rep(0, i), rep(0.0113, j), rep(0.0532, k), rep(0.1191, 20-i-j-k)), loc = el_par))^(1/20)))
    }
  }
}
head(hbest5[order(hbest5$val, decreasing = T), ])



exp_log_loc_mat <- function(d, c0, c1, b0, b1){
  d <- c(0, d)
  n <- length(d)
  w <- c(0.371, 0.501, 0.087, 0.041)
  mat_list <- lapply(1:n, function(x) w[x] * exp_log_f(d[x], c0, c1, b0, b1))
  
  #print(mat_list)
  M <- Reduce("+", mat_list)
  M
}

exp_log_loc_hoptimal <- function(d, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  pen = 999999999
  
  n <- length(d) + 1
  h1 <- -c1
  h2 <- -c0
  h3 <- b1 * exp(b0) * (1 - exp(b0)) / (exp(b0) + 1)^3
  h4 <- exp(b0) / (exp(b0) + 1)^2 
  h <- matrix(c(h1, h2, h3, h4))
  
  M <- exp_log_loc_mat(d, c0, c1, b0, b1)
  if (rcond(M) < 2.220446e-16){
    res = pen
  } 
  else {
    M_inv <- solve(M)
    res = t(h) %*% M_inv %*% h
  }
  
  res
}

efb <- exp_log_loc_hoptimal(d = c(0.0108, 0.0526, 0.1187), el_par)
efb

psoinfo <- psoinfo_setting(nSwarms = 256, Iters = 2000)
exp_log_tau <- exp_log_tau_pso(4, el_par, psoinfo)
exp_log_tau$points
exp_log_tau$weight

exp_log_tau_res <- list()

for(j in 4:10){
  exp_log_tau_res[[j]] <- exp_log_loc_tau_pso_rep(3, j, el_par, psoinfo)
  print(j)
  print(exp_log_tau_res[[j]]$design_points)
  print(exp_log_tau_res[[j]]$weight)
}

psoinfo <- psoinfo_setting(nSwarms = 512, Iters = 10000)
exp_log_tau_res1 <- list()
for(j in 7:10){
  exp_log_tau_res1[[j]] <- exp_log_tau_pso_rep(10, j, el_par, psoinfo)
  print(j)
  print(exp_log_tau_res1[[j]]$design_points)
  print(exp_log_tau_res1[[j]]$val)
}

tauop <- c()
tauop[19] <- exp_log_tauoptimal(d = c(rep(0,10), rep(0.04197718,10)), loc = c(el_par, 0.04197718))
tauoptimal <- exp_log_tauoptimal(d = c(rep(0,2), rep(0.04197718,2)), loc = c(el_par, 0.04197718)) / tauop
saveRDS(tauoptimal, "tauop.rds")

tauoptimall <- readRDS("tauop.rds")
tauoptimall
tau <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, 
               c0 = 0.15, c1 = 89, b0 = 3.2, b1 = 41)$root
tau

tauoptimall <- sapply(1:19, function(x) tauoptimal[x]^(1/(x+1)))
plot(x = c(2:20), y = tauoptimall, xlab = "N (Number of Observations)", ylab = "Efficiency", type = "l")
points(x = c(2:20), y = tauoptimall, pch = 19, cex = .75)

h_op <- list()
for (j in 16:20){
  if(j <= 10){
    psoinfo <- psoinfo_setting(nSwarms = 128, Iters = 2000)
    h_op[[j]] <- exp_log_pso(nPoints = j, parms = el_par, psoinfo = psoinfo)
  } 
  else if (j <= 15){
    psoinfo <- psoinfo_setting(nSwarms = 256, Iters = 10000)
    h_op[[j]] <- exp_log_pso(nPoints = j, parms = el_par, psoinfo = psoinfo)
  } 
  else {
    psoinfo <- psoinfo_setting(nSwarms = 1024, Iters = 30000)
    h_op[[j]] <- exp_log_pso(nPoints = j, parms = el_par, psoinfo = psoinfo)
  }
  print(psoinfo$nSwarm)
  print(j)
}

par(mfrow = c(2,2))
plot(x = c(0:2000), y = h_op[[5]]$history^(1/5), type = "l", 
     xlab = "Number of Iterations", ylab = "Efficiency", main = "5-Point Exact Design")
plot(x = c(0:2000), y = h_op[[10]]$history^(1/10), type = "l", 
     xlab = "Number of Iterations", ylab = "Efficiency", main = "10-Point Exact Design")
plot(x = c(0:10000), y = h_op[[15]]$history^(1/15), type = "l", 
     xlab = "Number of Iterations", ylab = "Efficiency", main = "15-Point Exact Design")
plot(x = c(0:30000), y = h_op[[20]]$history^(1/20), type = "l", 
     xlab = "Number of Iterations", ylab = "Efficiency", main = "20-Point Exact Design")




h_cpu <- lapply(4:20, function(x) h_op[[x]]$cputime)
plot(x = c(4:20), y = h_cpu, type = "l",
                xlab = "Number of Design Points", ylab = "CPU Runtime(sec.)")

d_op <- list()
for (j in 4:20){
  if(j <= 10){
    psoinfo <- psoinfo_setting(nSwarms = 64, Iters = 1000)
    d_op[[j]] <- hb_pso(nPoints = j, parms = el_par, psoinfo = psoinfo)
  } 
  else if (j <= 15){
    psoinfo <- psoinfo_setting(nSwarms =128, Iters = 5000)
    d_op[[j]] <- hb_pso(nPoints = j, parms = el_par, psoinfo = psoinfo)
  } 
  else {
    psoinfo <- psoinfo_setting(nSwarms = 256, Iters = 10000)
    d_op[[j]] <- hb_pso(nPoints = j, parms = el_par, psoinfo = psoinfo)
  }
  print(psoinfo$nSwarm)
  print(j)
}

par(mfrow = c(2,2))
plot(x = c(0:2000), y = h_op[[5]]$history^(1/5), type = "l", 
     xlab = "Number of Iterations", ylab = "Efficiency", main = "5-Point Exact Design")
plot(x = c(0:2000), y = h_op[[10]]$history^(1/10), type = "l", 
     xlab = "Number of Iterations", ylab = "Efficiency", main = "10-Point Exact Design")
plot(x = c(0:10000), y = h_op[[15]]$history^(1/15), type = "l", 
     xlab = "Number of Iterations", ylab = "Efficiency", main = "15-Point Exact Design")
plot(x = c(0:30000), y = h_op[[20]]$history^(1/20), type = "l", 
     xlab = "Number of Iterations", ylab = "Efficiency", main = "20-Point Exact Design")




d_cpu <- lapply(4:20, function(x) d_op[[x]]$cputime)
plot(x = c(4:20), y = d_cpu, type = "l",
     xlab = "Number of Design Points", ylab = "CPU Runtime(sec.)")


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
  
  pso_res
}

psoinfo <- psoinfo_setting(nSwarms = 128, Iters = 2000)
logistic_base <- logistic_pso(2, parms = c(2, 1), psoinfo = psoinfo)
sort(2.145 + 0.658 * logistic_best$par)
logistic_base$val
logistic_best_list <- list(points = list(), val = c())

for(k in 2:20){
  logistic_best <- logistic_pso(k, parms = c(2, 1), psoinfo = psoinfo)
  logistic_best_list$points[[k]] <- sort(2 + 1 * logistic_best$par)
  logistic_best_list$val <- c(logistic_best_list$val, (logistic_best$val/logistic_base$val)^(1/k))
}

logistic_best_list$val
plot(x = c(2:20), y = logistic_best_list$val, type = "l", 
     xlab = "N (Number of Observations)", ylab = "Efficiency")
points(x = c(2:20), y = logistic_best_list$val, pch = 19, cex = .75)

qlogsitic_mat <- function(){
  
}

qlogistic_doptimal1 <- function(d, loc){
  a <- loc[1]
  b <- loc[2]
  mu <- loc[3]
  n <- length(d)
  
  eta <- a + b * (x - mu)^2
  pi <- exp(eta) / (1 + exp(eta))
  v <- pi * (1 - pi)
  
  
  
}

qlogistic_doptimal2 <- function(d, loc){
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

qlogistic_doptimal3 <- function(d, loc){
  a <- loc[1]
  b1 <- loc[2]
  b2 <- loc[3]
  d <- d
  n <- length(d)
  
  theta <- a + b1 * d + b2 * d^2
  omega <- exp(theta) / (1 + exp(theta))^2
  w <- c(0.297, 0.203, 0.203, 0.297)
  
  inf_mat <- matrix(c(sum(w * omega), sum(w * d * omega), sum(w * d^2 * omega), 
                      sum(w * d * omega), sum(w * d^2 * omega), sum(w * d^3 * omega), 
                      sum(w * d^2 * omega), sum(w * d^3 * omega), sum(w * d^4 * omega)), 
                    nrow = 3)
  #print(inf_mat)
  
  -det(inf_mat)
}

qlogistic_doptimal4 <- function(d, loc){
  a <- loc[1]
  b1 <- loc[2]
  b2 <- loc[3]
  n <- length(d) * 3
  
  d1 <- rep(d,3)
  theta <- rep(a + b1 * d + b2 * d^2, 3)
  omega <- exp(theta) / (1 + exp(theta))^2
  
  inf_mat <- matrix(c(sum(omega), sum(d1 * omega), sum(d1^2 * omega), 
                      sum(d1 * omega), sum(d1^2 * omega), sum(d1^3 * omega), 
                      sum(d1^2 * omega), sum(d1^3 * omega), sum(d1^4 * omega)), 
                    nrow = 3)
  #print(inf_mat)
  
  -det(inf_mat) / n^3
}

qlogistic_pso1 <- function(npoints, parms, psoinfo){
  lb <- rep(-3, npoints-1)
  ub <- rep(3, npoints-1)
  
  pso_res <- globpso(objFunc = qlogistic_doptimal1, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res
}

qlogistic_pso2 <- function(npoints, parms, psoinfo){
  lb <- rep(-3, npoints)
  ub <- rep(3, npoints)
  
  pso_res <- globpso(objFunc = qlogistic_doptimal2, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res
}

qlogistic_pso4 <- function(npoints, parms, psoinfo){
  lb <- rep(-3, npoints)
  ub <- rep(3, npoints)
  
  pso_res <- globpso(objFunc = qlogistic_doptimal4, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res
}

psoinfo <- psoinfo_setting(nSwarms = 1024, Iters = 10000)

qlogistic_doptimal2(rep(c(-2.0876, -1.3938, 1.3938, 2.0876), 3), c(3, 0, -1))

qlogistic_best <- qlogistic_pso4(4, parms = c(3, 0, -1), psoinfo = psoinfo)
sort(qlogistic_best$par) |> round(4)
qlogistic_best$val

qlogistic_best1 <- list(points = list(), val = c())
for(j in 2:10){
  qlogistic_best <- qlogistic_pso1(j, parms = c(3, 0, -1), psoinfo = psoinfo)
  qlogistic_best1$points[[j]] <- sort(qlogistic_best$par) |> round(4)
  qlogistic_best1$val <- c(qlogistic_best1$val, qlogistic_best$val)
}  

qlogistic_best1$points
qlogistic_best1$val

qlogistic_best2 <- list()
psoinfo <- psoinfo_setting(nSwarms = 256, Iters = 3000)
for(j in 4:20){
  qlogistic_best2[[j]] <- list()
  qlogistic_best <- qlogistic_pso2(j, parms = c(3, 0 ,-1), psoinfo = psoinfo)
  qlogistic_best2[[j]]$points <- sort(qlogistic_best$par) |> round(4)
  qlogistic_best2[[j]]$val <- c(qlogistic_best2$val, qlogistic_best$val)
}  

qlogistic_best2

qlogistic_best <- qlogistic_pso2(12, parms = c(3, 0 ,-1), psoinfo = psoinfo)
qlogistic_best2[[12]]$points <- sort(qlogistic_best$par) |> round(4)
qlogistic_best2[[12]]$val <- c(qlogistic_best2$val, qlogistic_best$val)
qlogistic_best2[[12]]

qlogistic_best2_points <- lapply(4:20, function(x) qlogistic_best2[[x]]$points)
qlogistic_best2_val <- sapply(4:20, function(x) qlogistic_best2[[x]]$val)
qlogistic_best2_points
qlogistic_best2_val
qlogistic_base <- qlogistic_doptimal3(d = c(-2.061, -1.324, 1.324, 2.061), loc = c(3, 0, -1))
eff2 <- qlogistic_best2_val / qlogistic_base
eff2

d_best <- list(points = matrix(ncol = 4, nrow = 17), eff = c(0,0,0))
rownames(d_best$points) <- c(4:20)
colnames(d_best$points) <- c(-2.061, -1.324, 1.324, 2.061)
for(i in 4:20){
  res = 0
  points = c(0,0,0,0)
  for(a in 0:i){
    for(b in 0:(i-a)){
      for(c in 0:(i-a-b)){
        eff <- qlogistic_doptimal2(d = c(rep(-2.061, a), rep(-1.324, b), rep(1.324, c), rep(2.061, i-a-b-c)), 
                                    loc = c(3, 0, -1))
        if(eff < res){
          res = eff
          points <- c(a, b, c, i-a-b-c)
        }
      }
    }
  }
  d_best$points[i-3,] <- points
  d_best$eff[i] <- res
}
d_best$points
d_best$eff

d_best_eff <- sapply(4:20, function(x) (d_best$eff[x]/qlogistic_base)^(1/x))

plot(x = c(4:20), y = sapply(4:20, function(x) eff2[x-3]^(1/x)), type = "l", 
     xlab = "Number of Design Points", ylab = "Efficiency")
lines(x = c(4:20), y = d_best_eff, col = "red")
points(x = c(3:20), y = sapply(3:20, function(x) eff2[x-3]^(1/x)), pch = 19, cex = .75)


qlogistic_best01 <- list(points = list(), val = c())
for(j in 2:10){
  qlogistic_best <- qlogistic_pso1(j, parms = c(0, 0 ,-1), psoinfo = psoinfo)
  qlogistic_best01$points[[j]] <- sort(c(qlogistic_best$par,0)) |> round(4)
  qlogistic_best01$val <- c(qlogistic_best01$val, qlogistic_best$val)
}  

qlogistic_best01$points
qlogistic_best01$val



qlogistic_best02 <- list(points = list(), val = c())
for(j in 2:10){
  qlogistic_best <- qlogistic_pso2(j, parms = c(0, 0 ,-1), psoinfo = psoinfo)
  qlogistic_best02$points[[j]] <- sort(qlogistic_best$par) |> round(4)
  qlogistic_best02$val <- c(qlogistic_best02$val, qlogistic_best$val)
}  

qlogistic_best02$points
qlogistic_best02$val



qlogistic_best11 <- list(points = list(), val = c())
for(j in 3:20){
  qlogistic_best <- qlogistic_pso1(j, parms = c(-3, 0, -1), psoinfo = psoinfo)
  qlogistic_best11$points[[j]] <- sort(c(qlogistic_best$par,0)) |> round(4)
  qlogistic_best11$val <- c(qlogistic_best11$val, qlogistic_best$val)
}  
qlogistic_best11$points
qlogistic_best11$val



qlogistic_best12 <- list(points = list(), val = c())
for(j in 3:20){
  qlogistic_best <- qlogistic_pso2(j, parms = c(-3, 0, -1), psoinfo = psoinfo)
  qlogistic_best12$points[[j]] <- sort(qlogistic_best$par) |> round(4)
  qlogistic_best12$val <- c(qlogistic_best12$val, qlogistic_best$val)
}  
qlogistic_best <- qlogistic_pso2(9, parms = c(3, 0, -1), psoinfo = psoinfo)
qlogistic_best$par |> sort() |> round(4)
qlogistic_best$val
qlogistic_best12$points
qlogistic_best12$val
qlogistic_base1 <- qlogistic_best12$val[1]
eff1 <- qlogistic_best12$val / qlogistic_base1
sapply(3:20, function(x) eff1[x-2]^(1/x))
plot(x = c(3:20), y = sapply(3:20, function(x) eff1[x-2]^(1/x)), type = "l", 
     xlab = "Number of Design Points", ylab = "Efficiency")
points(x = c(3:20), y = sapply(3:20, function(x) eff1[x-2]^(1/x)), pch = 19, cex = .75)

(qlogistic_doptimal2(d = c(qlogistic_best$par, 1.237986), loc = c(-3, 0, -1)) / qlogistic_base1)^(1/4)

qlog1 <- c()
qlog1 <- c(qlog1, (qlogistic_doptimal(d = c(rep(0,7),rep(-1.238,7),rep(1.238,6)),loc = c(-3,0,-1)) / 
                     qlogistic_doptimal(c(-1.238, 0, 1.238),loc = c(-3,0,-1)))^(1/20))
plot(x = c(3:20), y = qlog1, type = "l", 
     xlab = "N (Number of Observations)", ylab = "Efficiency")
points(x = c(3:20), y = qlog1, pch = 19, cex = .75)

qlogistic_pso(19, c(3,0,-1), psoinfo, 5)$par

qlog2 <- c()
qlog2 <- c(qlog2, (qlogistic_doptimal(d = c(rep(-2.0592,6),rep(-1.3187,4),rep(1.3187,4),rep(2.0592,6)),loc = c(3,0,-1)) / 
                     qlogistic_doptimal_approx(input = c(-2.061,-1.324,1.324,2.061,0.297,0.203,0.203),loc = c(3,0,-1)))
           ^(1/20))
qlog2
plot(x = c(4:20), y = qlog2, type = "l", 
     xlab = "N (Number of Observations)", ylab = "Efficiency")
points(x = c(4:20), y = qlog2, pch = 19, cex = .75)


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
  #print(inf_mat)
  
  -det(inf_mat) / n^4
}

clogistic_doptimal_loc <- function(input, loc){
  a <- loc[1]
  b1 <- loc[2]
  b2 <- loc[3]
  b3 <- loc[4]
  #print(c(a, b1, b2, b3))
  n <- 4
  d <- input[1:4]
  weight <- input[5:7]
  weight[4] <- 1- sum(weight)
  #print(weight)
  
  theta <- a + b1 * d + b2 * d^2 + b3 * d^3
  omega <- exp(theta) / (1 + exp(theta))^2
  #print(omega)
  
  inf_mat <- matrix(c(sum(weight * omega), sum(weight * d * omega), 
                      sum(weight * d^2 * omega), sum(weight * d^3 * omega), 
                      sum(weight * d * omega), sum(weight * d^2 * omega), 
                      sum(weight * d^3 * omega),  sum(weight * d^4 * omega), 
                      sum(weight * d^2 * omega), sum(weight * d^3 * omega), 
                      sum(weight * d^4 * omega), sum(weight * d^5 * omega), 
                      sum(weight * d^3 * omega), sum(weight * d^4 * omega), 
                      sum(weight * d^5 * omega), sum(weight * d^6 * omega)), 
                    nrow = 4)
  #print(inf_mat)
  
  if (weight[4] < 0){
    res = 9999999
  } else {
    res = -det(inf_mat)
  }
  
  res
}

clogistic_pso <- function(npoints, parms, psoinfo){
  lb <- rep(-5, npoints)
  ub <- rep(5, npoints)
  
  pso_res <- globpso(objFunc = clogistic_doptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res
}

clogistic_loc_pso <- function(parms, psoinfo){
  lb <- c(rep(-5, 4), rep(0, 3))
  ub <- c(rep(5, 4), rep(1, 3))
  
  pso_res <- globpso(objFunc = clogistic_doptimal_loc, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res
}

psoinfo <- psoinfo_setting(nSwarms = 128, Iters = 1000)
clogistic_list1 <- list(points = list(), val = c())
for(j in 4:12){
  clogistic_best <- clogistic_pso(j, parms = c(-3, 0, 0, -1), psoinfo = psoinfo)

  clogistic_list1$points[[j]] <- sort(clogistic_best$par)
  clogistic_list1$val <- c(clogistic_list1$val, clogistic_best$val)
}  
clogistic_list1$points
clogistic_list1$val

clogistic1_base <- clogistic_list1$val[1]
clogistic1_eff <- sapply(4:12, function(x) (clogistic_list1$val[x-3]/clogistic1_base)^(1/x))
plot(x = c(4:12), y = clogistic1_eff, type = "l", 
     xlab = "N (Number of Observations)", ylab = "Efficiency")
points(x = c(4:12), y = clogistic1_eff, pch = 19, cex = .75)

-3 - clogistic_list1$points[[4]]^3

clogistic_list2 <- list(points = list(), val = c())
for(j in 4:12){
  clogistic_best <- clogistic_pso(j, parms = c(-3, 0, 0, 1), psoinfo = psoinfo)
  
  clogistic_list2$points[[j]] <- sort(clogistic_best$par) |> round(4)
  clogistic_list2$val <- c(clogistic_list2$val, clogistic_best$val)
}  
clogistic_list2$points
clogistic_list2$val

-3 + clogistic_list2$points[[4]]^3

clogistic_list4 <- list(points = list(), val = c())
for(j in 4:12){
  clogistic_best <- clogistic_pso(j, parms = c(0, 0, 0, 1), psoinfo = psoinfo)
  
  clogistic_list4$points[[j]] <- sort(clogistic_best$par) |> round(4)
  clogistic_list4$val <- c(clogistic_list4$val, clogistic_best$val)
}  
clogistic_list4$points
clogistic_list4$val

clogistic_list4$points[[4]]^3

psoinfo <- psoinfo_setting(nSwarms = 512, Iters = 20000)
clogistic_list <- list(points = list(), val = c())
for(j in 14:20){
  clogistic_best <- clogistic_pso_rep(5, j, c(1, 3, 2, -1), psoinfo, 5)
  
  clogistic_list$points[[j]] <- sort(clogistic_best$result$design_points) |> round(4)
  clogistic_list$val <- c(clogistic_list$val, clogistic_best$val)
}
clogistic_list$points
clogistic_list$val

clogistic_doptimal_loc(input = c(-1.1495, -0.4978, 0.1262, 2.9709, 3.1796, 0.2, 0.2, 0.2, 0.2), 
                       loc = c(1, 3, 2, -1))

psoinfo <- psoinfo_setting(nSwarms = 256, Iters = 10000)
clogistic_loc_res <- clogistic_loc_pso(parms = c(1, 3, 2, -1), psoinfo = psoinfo)
clogistic_loc_res$par |> round(4)
clogistic_loc_res$val

clogistic3_base <- clogistic_loc_res$val
clogistic3_eff <- sapply(5:15, function(x) (clogistic_list3$val[x-4]/clogistic3_base)^(1/x))
plot(x = c(5:15), y = clogistic3_eff, type = "l", 
     xlab = "Number of Design Points", ylab = "Efficiency")
points(x = c(5:15), y = clogistic3_eff, pch = 19, cex = .75)

clogistic_loc_res1 <- clogistic_loc_pso(parms = c(-3, 0, 0, -1), psoinfo = psoinfo)
clogistic_loc_res1$par
clogistic_loc_res1$val

clog2_base <- clogistic_doptimal_approx(input = c(-1.1073,-0.3497,0.0005,2.9576,3.1874,0.2454,0.1082,0.1669,0.2377), 
                                        loc = c(1,3,2,-1))

clog2 <- c()
clog2 <- c(clog2, (clogistic_doptimal(d = c(rep(-1.1328,3),rep(-0.4231,2),rep(0.1190,2),rep(2.9604,2),rep(3.1779,3)),loc = c(1,3,2,-1)) / 
                     clog2_base)
           ^(1/12))
clog2
plot(x = c(5:12), y = clog2, type = "l", 
     xlab = "N (Number of Observations)", ylab = "Efficiency")
points(x = c(5:12), y = clog2, pch = 19, cex = .75)

clogistic_get <- function(x, parms, d, w){
  a <- parms[1]
  b1 <- parms[2]
  b2 <- parms[3]
  b3 <- parms[4]
  
  fx <- matrix(c(1, x, x^2, x^3))
  wx <- exp(a + b1 * x + b2 * x^2 + b3 * x^3) / (1 + exp(a + b1 * x + b2 * x^2 + b3 * x^3))^2
  
  theta <- a + b1 * d + b2 * d^2 + b3 * d^3
  omega <- exp(theta) / (1 + exp(theta))^2
  
  inf_mat <- matrix(c(sum(w * omega), sum(w * d * omega), sum(w * d^2 * omega), sum(w * d^3 * omega), 
                      sum(w * d * omega), sum(w * d^2 * omega), sum(w * d^3 * omega),  sum(w * d^4 * omega), 
                      sum(w * d^2 * omega), sum(w * d^3 * omega), sum(w * d^4 * omega), sum(w * d^5 * omega), 
                      sum(w * d^3 * omega), sum(w * d^4 * omega), sum(w * d^5 * omega), sum(w * d^6 * omega)), 
                    nrow = 4)
  mat_inv <- solve(inf_mat)
  
  wx * t(fx) %*% mat_inv %*% fx
}


clogistic1 <- clogistic_pso(4, parms = c(-3, 0, 0, -1), psoinfo = psoinfo)
clogistic1$par |> round(4)


eval_reg1 <- seq(from = -2, to = 2, by = 1e-4)
clogistic_get1 <- sapply(1:40001, function(i) clogistic_get(x = eval_reg1[i], parms = c(-3, 0, 0, -1), d = clogistic1$par, 
                                                             w = rep(1/4, 4)))
plot(x = eval_reg1, y = clogistic_get1, type = "l", xlab = "x")
abline(h = 4, lty = "dashed")
points(x = clogistic1$par |> round(4), y = rep(4, 4), cex = 1.2)


eval_reg2 <- seq(from = -1.5, to = 3.5, by = 1e-4)
clogistic_get2 <- sapply(1:50001, function(i) clogistic_get(x = eval_reg2[i], parms = c(1, 3, 2, -1), 
                                                             d = c(-1.1073, -0.3497, 0.0005, 2.9576, 3.1874), 
                                                             w = c(0.2454, 0.1082, 0.1669, 0.2377, 0.2418)))
plot(x = eval_reg2, y = clogistic_get2, type = "l", xlab = "x", ylab = TeX('$\\bar{d}(w, \\xi)$'))
abline(h = 4, lty = "dashed")
points(x = c(-1.1073, -0.3497, 0.0005, 2.9576, 3.1874), y = rep(4, 5), cex = 1.2)



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

hb_tauoptimal_pso <- function(nPoints, parms, psoinfo){
  lb <- rep(0, nPoints)
  ub <- rep(0.15, nPoints)
  
  pso_res <- globpso(objFunc = hb_tauoptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$par <- pso_res$par |> round(4) |> sort()
  pso_res
}

hb_tau_res <- list()
psoinfo <- psoinfo_setting(nSwarms = 256, Iters = 4000)
for(j in 2:20){
  hb_tau_pso <- hb_tauoptimal_pso(j, parms = c(170, 0.04, 1.46, 40), psoinfo = psoinfo)
  hb_tau_res[[j]] <- hb_tau_pso
  print(hb_tau_pso$val)
  print(hb_tau_pso$par)
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

hb_hoptimal_loc <- function(d, w, loc){
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  n <- length(d)
  pen = 999999999
  
  mat_list <- lapply(1:n, function(x) w[x] * hb_mat(d[x], c1, tau, b0, b1))
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

hb_hoptimal_pso <- function(nPoints, parms, psoinfo){
  lb <- rep(0, nPoints)
  ub <- rep(0.15, nPoints)
  
  pso_res <- globpso(objFunc = hb_hoptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$par <- pso_res$par |> round(4) |> sort()
  pso_res
}

hb_h_res <- list()
psoinfo <- psoinfo_setting(nSwarms = 1024, Iters = 30000)
for(j in 18:20){
  hb_h_pso <- hb_hoptimal_pso(j, parms = c(170, 0.04, 1.46, 40), psoinfo = psoinfo)
  hb_h_res[[j]] <- hb_h_pso
  print(hb_h_pso$val)
  print(hb_h_pso$par)
}


hb_h_res[[18]]$par
hb_h_base <- hb_hoptimal_loc(d = c(0, 0.020, 0.040), w = c(0.359, 0.5, 0.141), 
                             loc = c(170, 0.04, 1.46, 40))

hb_h_res1 <- list()
psoinfo <- psoinfo_setting(nSwarms = 512, Iters = 10000)
for(j in 16){
  hb_h_pso <- hb_hoptimal_pso(j, parms = c(170, 0.04, 1.86, 40), psoinfo = psoinfo)
  hb_h_res1[[j]] <- hb_h_pso
  print(hb_h_pso$val)
  print(hb_h_pso$par)
}

for(j in 4:20){
  print(hb_h_res[[j]]$par)
}

hb_h_point <- lapply(4:20, function(x) hb_h_res1[[x]]$par)
hb_h_val <- sapply(4:20, function(x) hb_h_res1[[x]]$val)
hb_h_base <- hb_hoptimal_loc(d = c(0, 0.020, 0.0479, 0.125), w = c(0.3174, 0.5029, 0.1644, 0.0153), 
                             loc = c(170, 0.04, 1.86, 40))
hb_h_res1[[16]] <- hb_hoptimal(d = c(rep(0, 5), rep(0.0207, 7), rep(0.0534, 2), rep(0.1314, 1)), loc = c(170, 0.04, 1.86, 40))
hb_h_val <- sapply(c(4:20), function(x) if(x == 15 || x == 16){
  res = hb_h_res1[[x]]
} else {
  res = hb_h_res1[[x]]$val
}
)

hb_h_eff <- sapply(4:20, function(x) (hb_h_base / hb_h_val[x-3]) ^ (1/x))
hb_h_eff
plot(x = c(4:20), y = hb_h_eff, type = "l", xlab = "Number of Design Points", ylab = "Efficiency")
points(x = c(4:20), y = hb_h_eff, pch = 19, cex = .75)

hb_best <- list(points = matrix(ncol = 4, nrow = 17), val = c(0,0,0))
rownames(hb_best$points) <- c(4:20)
colnames(hb_best$points) <- c(0, 0.020, 0.0479, 0.125)
for(i in 4:20){
  res = 9999999999
  points = c(0,0,0,0)
  for(a in 0:i){
    for(b in 0:(i-a)){
      for(c in 0:(i-a-b)){
        val <- hb_hoptimal(d = c(rep(0, a), rep(0.020, b), rep(0.0479, c), rep(0.125, i-a-b-c)), 
                                   loc = c(170, 0.04, 1.86, 40))
        if(val < res){
          res = val
          points <- c(a, b, c, i-a-b-c)
        }
      }
    }
  }
  hb_best$points[i-3,] <- points
  hb_best$val[i] <- res
}
hb_best$points
hb_best$val

hb_best_eff <- sapply(4:20, function(x) (hb_h_base / hb_best$val[x]) ^ (1/x))

plot(x = c(4:20), y = hb_best_eff, type = "l", xlab = "Number of Design Points", ylab = "Efficiency", col = "red")
lines(x = c(4:20), y = hb_h_eff, col = "blue")
points(x = c(4:20), y = hb_h_eff, pch = 16, cex = .75, col = "blue")
points(x = c(4:20), y = hb_best_eff, pch = 17, cex = .75, col = "red")
legend("bottomright", c("Best Assignment", "PSO"), lty = c(1, 3), pch = c(16, 17), 
       col = c("red", "blue"))

hbtau_base <- hb_tauoptimal(d = c(0, 0.04), loc = c(170, 0.04, 1.46, 40))
hbtau_eff <- c()
hbtau_eff <- c(hbtau_eff, (hbtau_base / 
                             hb_tauoptimal(d = c(rep(0,10), rep(0.04,10)), loc = c(170, 0.04, 1.46, 40)))^(1/20))
hbtau_eff 
plot(x = c(2:20), y = hbtau_eff, type = "l", 
     xlab = "N (Number of Observations)", ylab = "Efficiency")
points(x = c(2:20), y = hbtau_eff, pch = 19, cex = .75)


hbh_base <- hb_hoptimal_approx(input = c(0, 0.02, 0.0479, 0.125, 0.3174, 0.5029, 0.1644), loc = c(170, 0.04, 1.86, 40))
hbh_eff <- c()
hbh_eff <- c(hbh_eff, (hbh_base / 
                             hb_hoptimal(d = c(rep(0,6), rep(0.0203,10), rep(0.0534,3), rep(0.1314,1)), loc = c(170, 0.04, 1.86, 40)))
               ^(1/20))
hbh_eff
plot(x = c(4:20), y = hbh_eff, type = "l", 
     xlab = "N (Number of Observations)", ylab = "Efficiency")
points(x = c(4:20), y = hbh_eff, pch = 19, cex = .75)

exp_log_doptimal <- function(d, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  
  M <- exp_log_mat(d, c0, c1, b0, b1)
  
  -det(M)
}

exp_log_doptimal_pso <- function(nPoints, parms, psoinfo){
  lb <- rep(0, nPoints)
  ub <- rep(0.15, nPoints)
  
  pso_res <- globpso(objFunc = exp_log_doptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res$val <- pso_res$val * -1
  pso_res$par <- pso_res$par |> round(4) |> sort()
  pso_res
}


el_d_res <- list()
psoinfo_hb <- psoinfo_setting(nSwarms = 128, Iters = 1000)
for(j in 4:20){
  el_d_pso <- exp_log_doptimal_pso(j, parms = c(0.15, 89, 3.2, 41), psoinfo = psoinfo, upper = 0.15)
  el_d_res[[j]] <- el_d_pso
  print(el_d_pso$val)
  print(el_d_pso$par)
}

for(j in 4:20){
  print(el_d_res[[j]]$par)
}

eld_base <- exp_log_doptimal(d = c(rep(0,1), rep(0.0109,1), rep(0.0558,1), rep(0.1051,1)), loc = c(0.15, 89, 3.2, 41))
eld_eff <- c()
eld_eff <- c(eld_eff, (exp_log_doptimal(d = c(rep(0,5), rep(0.0109,5), rep(0.0558,5), rep(0.1051,5)), loc = c(0.15, 89, 3.2, 41))
                       / eld_base) ^ (1/20))
eld_eff

plot(x = c(4:20), y = eld_eff, type = "l", 
     xlab = "N (Number of Observations)", ylab = "Efficiency")
points(x = c(4:20), y = eld_eff, pch = 19, cex = .75)


library(CVXR)

hb_test <- hb_doptimal_pso_rep(10, 7, hb_par, psoinfo_hb)
hb_test$result
unique(hb_test$result$design_points)
hb_test_ed <- hb_test$result$design_points |> table() |> data.frame() 
hb_test_ed
colnames(hb_test_ed) <- c("Support Points", "Replications")

qlog_par = c(-3,0,-1)
qlog_test <- qlogistic_pso_rep(1, 4, qlog_par, psoinfo_hb)
qlog_test$approximate_design
if (nrow(qlog_test$approximate_design) != length(unique(qlog_test$approximate_design$support_points))){
  cnt <- qlog_test$approximate_design |> count(support_points)
  ext <- cnt$support_points[which(cnt$n != 1)]
  apprep <- qlog_test$approximate_design$support_points == ext
  w <- sum(qlog_test$approximate_design$weight[apprep])
  qlog_test$approximate_design <- qlog_test$approximate_design[!apprep,]
  qlog_test$approximate_design <- rbind(qlog_test$approximate_design, c(ext, w))
}
qlog_test$approximate_design

qlog_test$approximate_design$weight[-apprep,]

print(hb_test$approximate_design, row.names = F)
cat("Support Point:", hb_test$approximate_design$support_points, "\n", 
    "Weight:", hb_test$approximate_design$weight)

hbd <- hb_doptimal_approx(hb_par)
plot(x = c(1:1001), y = hbd$history, type = "l")

library(ggplot2)
hunt_bowman <- function(d, c1, tau, b0, b1){
  if (d <= tau){
    res = c1 * d^2 - c1 * tau * d + 1 / (1 + exp(b0))
  } else {
    res = 1 / (1 + exp(b0 - b1 * (d - tau)))
  }
}
fp <- seq(0, 0.15, by = 0.001)
hb_df <- data.frame(dose = fp, response = sapply(fp, function(x) hunt_bowman(x, 170, 0.04, 1.16, 40)))
ggplot(data = hb_df, aes(x = dose, y = response)) + 
  geom_line()

exp_log <- function(d, c0, c1, b0, b1){
  c0 * exp(-c1 * d) + 1 / (1 + exp(b0 - b1 * d))
}
el_df <- data.frame(dose = fp, response = sapply(fp, function(x) exp_log(x, 0.15, 89, 3.2, 41)))
plot(x = el_df$x, y = el_df$y, type = "l")
ggplot(data = el_df, aes(x = dose, y = response)) + 
  geom_line()

logistic <- function(d, alpha, beta){
  eta <- alpha + beta * d
  exp(eta) / (1 + exp(eta))
}
fp <- seq(-3, 3, by = 6/100)
log_df <- data.frame(dose = fp, response = sapply(fp, function(x) logistic(x, -1, 2)))
ggplot(data = log_df, aes(x = dose, y = response)) + 
  geom_line()

qlogistic <- function(d, alpha, beta1, beta2){
  eta <- alpha + beta1 * d + beta2 * d^2
  exp(eta) / (1 + exp(eta))
}


clogistic <- function(d, alpha, beta1, beta2, beta3){
  eta <- alpha + beta1 * d + beta2 * d^2 + beta3 * d^3
  exp(eta) / (1 + exp(eta))
}


psoinfo_qlog <- psoinfo_setting(nSwarms = 128, Iters = 1000)
###
qlog_par1 <- qlogistic_params(alpha = 3, beta1 = 0, beta2 = -1)
qlog_par2 <- qlogistic_params(alpha = 2, beta1 = 0, beta2 = -1)
qlog_par3 <- qlogistic_params(alpha = 1.32, beta1 = 0, beta2 = -1)
qlog_par4 <- qlogistic_params(alpha = 1.33, beta1 = 0, beta2 = -1)

qlog_app1 <- qlogistic_approx_pso(qlog_par1, bd = 5)
qlog_app2 <- qlogistic_approx_pso(qlog_par2, bd = 5)
qlog_app3 <- qlogistic_approx_pso(qlog_par3, bd = 5)
qlog_app4 <- qlogistic_approx_pso(qlog_par4, bd = 5)
qlog_app2$par
qlog_app3$par
qlog_app4$par

qlog_appcomp <- (-1 * qlogistic_doptimal_approx(input = qlog_app2$par[1:7], loc = qlog_par1) / qlog_app1$val) ^ (1/4)
qlog_appcomp

qlog_appcomp2 <- (-1 * qlogistic_doptimal_approx(input = qlog_app3$par[1:7], loc = qlog_par4) / qlog_app4$val) ^ (1/4)
qlog_appcomp2

#qlog_comp <- list()

qlog_true <- c()
qlog_mis <- c()

for(j in 3:15){
  qlog_ex1 <- qlogistic_pso(j, qlog_par1, psoinfo_qlog, 5)
  qlog_ex2 <- qlogistic_pso(j, qlog_par2, psoinfo_qlog, 5)
  qlog_true <- c(qlog_true, (-1 * qlogistic_doptimal(d = qlog_ex1$par, loc = qlog_par1) / qlog_app1$val) ^ (1/j))
  qlog_mis <- c(qlog_mis, (-1 * qlogistic_doptimal(d = qlog_ex2$par, loc = qlog_par1) / qlog_app1$val) ^ (1/j))
}

plot(x = c(4:15), y = sapply(4:15, function(x) qlog_comp[[x]]), type = "l", col = "blue", ylim = c(0.8999, 0.97), 
     xlab = "Number of Observations", ylab = "Relative Efficiency", main = "D-optimal Designs for Quadratic Logistic Models")
abline(h = qlog_appcomp, col = "red")
legend("topleft", c("Exact designs", "Approximate design"), col = c("blue", "red"), lty = c(1, 1))

qlog_comp2 <- sapply(4:15, function(x) qlog_comp[[x]] * (qlog_ex1$val / qlog_app1$val) ^ (1/x))
qlog_comp2

plot(x = c(4:15), y = qlog_mis[-1], type = "l", col = "red", ylim = c(0.9, 1),
     xlab = "Number of Observations", ylab = "Relative Efficiency", main = "D-optimal Designs for Quadratic Logistic Models")
lines(x = c(4:15), y = qlog_true[-1], col = "blue")
legend("bottomright", c("True design", "Misspecified design"), col = c("blue", "red"), lty = c(1, 1))

###
log_par1 <- logistic_params(alpha = 2, beta = 1)
log_par2 <- logistic_params(alpha = 1, beta = 1)

log_app1 <- logistic_doptimal_approx(log_par1, bound = 5)
log_app2 <- logistic_doptimal_approx(log_par2, bound = 5)
log_appcomp <- (-1 * logistic_doptimal(d = log_app2$par, loc = log_par1) / log_app1$val) ^ (1/2)
log_appcomp

#log_comp <- list()

for(j in 2:15){
  log_ex1 <- logistic_doptimal_pso(j, log_par1, psoinfo_log, 5)
  log_ex2 <- logistic_doptimal_pso(j, log_par2, psoinfo_log, 5)
  log_comp[[j]] <- (-1 * logistic_doptimal(d = log_ex2$par, loc = log_par1) / log_ex1$val) ^ (1/j)
  print(log_comp[[j]])
}
log_comp

plot(x = c(2:15), y = sapply(2:15, function(x) log_comp[[x]]), type = "l", col = "blue",  
     xlab = "Number of Observations", ylab = "Relative Efficiency", main = "D-optimal Designs for Logistic Models")
abline(h = log_appcomp, col = "red")
legend("topleft", c("Exact designs", "Approximate design"), col = c("blue", "red"), lty = c(1, 1))

###
clog_par1 <- clogistic_params(-3, 0, 0, -1)
clog_par2 <- clogistic_params(-2, 0, 0, -1)

clog_app1 <- clogistic_approx_pso(clog_par1, bd = 5)
clog_app2 <- clogistic_approx_pso(clog_par2, bd = 5)
clog_appcomp <- (-1 * clogistic_doptimal_approx(input = clog_app2$par[1:9], loc = clog_par1) / clog_app1$val) ^ (1/4)
clog_appcomp

#clog_comp <- list()

for(j in 4:15){
  clog_ex1 <- clogistic_pso(j, clog_par1, psoinfo_clog, 5)
  clog_ex2 <- clogistic_pso(j, clog_par2, psoinfo_clog, 5)
  clog_comp[[j]] <- (-1 * clogistic_doptimal(d = clog_ex2$par, loc = clog_par1) / clog_ex1$val) ^ (1/j)
  print(clog_comp[[j]])
}
clog_comp

plot(x = c(4:15), y = sapply(4:15, function(x) clog_comp[[x]]), type = "l", col = "blue",  
     xlab = "Number of Observations", ylab = "Relative Efficiency", main = "D-optimal Designs for Cubic Logistic Models")
abline(h = clog_appcomp, col = "red")
legend("topleft", c("Exact designs", "Approximate design"), col = c("blue", "red"), lty = c(1, 1))

clog_comp2 <- sapply(4:15, function(x) clog_comp[[x]] * (clog_ex1$val / clog_app1$val) ^ (1/x))
clog_comp2

plot(x = c(4:15), y = clog_comp2, type = "l", col = "red",  
     xlab = "Number of Observations", ylab = "Relative Efficiency", main = "D-optimal Designs for Cubic Logistic Models")
lines(x = c(4:15), y = sapply(4:15, function(x) clog_comp[[x]]), col = "blue")
legend("topleft", c("True design", "Misspecified design"), col = c("blue", "red"), lty = c(1, 1))

hb_par1 <- hb_parms(c1 = 170, tau = 0.04, b0 = 1.46, b1 = 40) 
hb_par2 <- hb_parms(c1 = 171, tau = 0.05, b0 = 1.47, b1 = 41) 

###
hbd_app1 <- hb_doptimal_approx(hb_par1, 0.15)
hbd_app2 <- hb_doptimal_approx(hb_par2, 0.15)
hbd_appcomp <- (-1 * hb_doptimal(d = hbd_app2$par, loc = hb_par1) / hbd_app1$val) ^ (1/4)
hbd_appcomp

#hbd_comp <- list()

for(j in 4:15){
  hbd_ex1 <- hb_doptimal_pso(j, hb_par1, psoinfo_hb, 0.15)
  hbd_ex2 <- hb_doptimal_pso(j, hb_par2, psoinfo_hb, 0.15)
  hbd_comp[[j]] <- (-1 * hb_doptimal(d = hbd_ex2$par, loc = hb_par1) / hbd_ex1$val) ^ (1/j)
}
print(hbd_comp)

plot(x = c(4:15), y = sapply(4:15, function(x) hbd_comp[[x]]), type = "l", col = "blue",  
     xlab = "Number of Observations", ylab = "Relative Efficiency", main = "D-optimal Designs for Hunt-Bowman Models")
abline(h = hbd_appcomp, col = "red")
legend("topleft", c("Exact designs", "Approximate design"), col = c("blue", "red"), lty = c(1, 1))

###
hbt_app1 <- hb_tauoptimal_approx(hb_par1, 0.15)
hbt_app2 <- hb_tauoptimal_approx(hb_par2, 0.15)
hbt_appcomp <- (hbt_app1$val / hb_tauoptimal(d = hbt_app2$par, loc = hb_par1)) ^ (1/2)
hbt_appcomp

#hbt_comp <- list()

for(j in 11:12){
  hbt_ex1 <- hb_tauoptimal_pso(j, hb_par1, psoinfo_hb, 0.15)
  hbt_ex2 <- hb_tauoptimal_pso(j, hb_par2, psoinfo_hb, 0.15)
  hbt_comp[[j]] <- (hbt_ex1$val / hb_tauoptimal(d = hbt_ex2$par, loc = hb_par1)) ^ (1/j)
  print(hbt_comp[[j]])
}
print(hbt_comp)

plot(x = c(2:15), y = sapply(2:15, function(x) hbt_comp[[x]]), type = "l", col = "blue",  
     xlab = "Number of Observations", ylab = "Relative Efficiency", main = "tau-optimal Designs for Hunt-Bowman Models")
abline(h = hbt_appcomp, col = "red")
legend("topleft", c("Exact designs", "Approximate design"), col = c("blue", "red"), lty = c(1, 1))



###
hbh_app1 <- hb_hoptimal_approx_pso(4, hb_par1, 0.15)
hbh_app2 <- hb_hoptimal_approx_pso(4, hb_par2, 0.15)
hbh_appcomp <- (hbh_app1$val / hb_hoptimal_approx(input = hbh_app2$par[1:7], loc = hb_par1)) ^ (1/3)
hbh_appcomp

#hbh_comp <- list()

for(j in 4:15){
  hbh_ex1 <- hb_hoptimal_pso(j, hb_par1, psoinfo_hb, 0.15)
  hbh_ex2 <- hb_hoptimal_pso(j, hb_par2, psoinfo_hb, 0.15)
  hbh_comp[[j]] <- (hbh_ex1$val / hb_hoptimal(d = hbh_ex2$par, loc = hb_par1)) ^ (1/j)
}
print(hbh_comp)

plot(x = c(4:15), y = sapply(4:15, function(x) qlog_comp[[x]]), type = "l", col = "blue", ylim = c(0.8999, 0.97), 
     xlab = "Number of Observations", ylab = "Relative Efficiency", main = "h-optimal Designs for Hunt-Bowman Models")
abline(h = qlog_appcomp, col = "red")
legend("topleft", c("Exact designs", "Approximate design"), col = c("blue", "red"), lty = c(1, 1))

hbh_true <- c()
hbh_mis <- c()

for(j in 4:15){
  hbh_ex1 <- hb_hoptimal_pso(j, hb_par1, psoinfo_hb, 0.15)
  hbh_ex2 <- hb_hoptimal_pso(j, hb_par2, psoinfo_hb, 0.15)
  hbh_true <- c(hbh_true, (hbh_app1$val / hbh_ex1$val) ^ (1/j))
  hbh_mis <- c(hbh_mis, (hbh_app1$val / hb_hoptimal(d = hbh_ex2$par, loc = hb_par1)) ^ (1/j))
}

plot(x = c(4:15), y = hbh_mis, type = "l", col = "red", ylim = c(0.85, 1),
     xlab = "Number of Observations", ylab = "Relative Efficiency", main = "h-optimal Designs for Hunt-Bowman Models")
lines(x = c(4:15), y = hbh_true, col = "blue")
legend("bottomright", c("True design", "Misspecified design"), col = c("blue", "red"), lty = c(1, 1))

###
el_par1 <- exp_log_params(0.15, 89, 3.2, 41) 
el_par2 <- exp_log_params(0.15, 110, 2.4, 30) 

#eld_comp <- list()

eld_app1 <- exp_log_doptimal_approx(el_par1, 0.15)
eld_app2 <- exp_log_doptimal_approx(el_par2, 0.15)
eld_appcomp <- (-1 * exp_log_doptimal(d = eld_app2$par, loc = el_par1) / eld_app1$val) ^ (1/4)

for(j in 7:15){
  eld_ex1 <- exp_log_doptimal_pso(j, el_par1, psoinfo_el, 0.15)
  eld_ex2 <- exp_log_doptimal_pso(j, el_par2, psoinfo_el, 0.15)
  eld_comp[[j]] <- (-1 * exp_log_doptimal(d = eld_ex2$par, loc = el_par1) / eld_ex1$val) ^ (1/j)
  #print(eld_comp[[j]])
}

print(eld_comp)

plot(x = c(4:15), y = sapply(4:15, function(x) eld_comp[[x]]), type = "l", col = "blue",  
     xlab = "Number of Observations", ylab = "Relative Efficiency", main = "D-optimal Designs for exp-log Models")
abline(h = eld_appcomp, col = "red")
legend("topleft", c("Exact designs", "Approximate design"), col = c("blue", "red"), lty = c(1, 1))

###
#elt_comp <- list()

elt_app1 <- exp_log_tauoptimal_approx(el_par1, 0.15)
elt_app2 <- exp_log_tauoptimal_approx(el_par2, 0.15)
tau_app <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, 
                   c0 = el_par1[1], c1 = el_par1[2], b0 = el_par1[3], b1 = el_par1[4])$root
elt_appcomp <- (elt_app1$val / exp_log_tauoptimal(d = elt_app2$par, loc = c(el_par1, tau))) ^ (1/2)
elt_appcomp
uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, 
        c0 = el_par2[1], c1 = el_par2[2], b0 = el_par2[3], b1 = el_par2[4])$root

for(j in 2:15){
  elt_ex1 <- exp_log_tauoptimal_pso(j, el_par1, psoinfo_el, 0.15)
  elt_ex2 <- exp_log_tauoptimal_pso(j, el_par2, psoinfo_el, 0.15)
  elt_comp[[j]] <- (elt_ex1$val / exp_log_tauoptimal(d = elt_ex2$par, loc = c(el_par1, tau_app))) ^ (1/j)
  #print(elt_comp[[j]])
}
tau2 <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, 
                c0 = el_par2[1], c1 = el_par2[2], b0 = el_par2[3], b1 = el_par2[4])$root
elt_ex1 <- exp_log_tauoptimal(d = c(rep(0, 8), rep(tau_app, 7)), c(el_par1, tau_app))
elt_ex2 <- exp_log_tauoptimal(d = c(rep(0, 8), rep(tau2, 7)), c(el_par1, tau_app))
elt_comp[[15]] <- (elt_ex1 / elt_ex2) ^ (1/15)

print(elt_comp)

plot(x = c(2:15), y = sapply(2:15, function(x) elt_comp[[x]]), type = "l", col = "blue",  
     xlab = "Number of Observations", ylab = "Relative Efficiency", main = "tau-optimal Designs for exp-log Models")
abline(h = elt_appcomp, col = "red")
legend("topleft", c("Exact designs", "Approximate design"), col = c("blue", "red"), lty = c(1, 1))


exp_log_tauoptimal(d = c(0, 0.0871), loc = el_par1)


elt_true <- c()
elt_mis <- c()
elt_true <- c(elt_true,  (exp_log_tauoptimal(d = c(rep(0, 8), rep(tau_app, 7)), c(el_par1, tau_app)) / elt_app1$val)^(1/15))
elt_mis <- c(elt_mis,  (exp_log_tauoptimal(d = c(rep(0, 8), rep(tau2, 7)), c(el_par1, tau_app)) / elt_app1$val)^(1/15))

plot(x = c(2:15), y = 1 / elt_mis, type = "l", col = "red", ylim = c(0, 1),
     xlab = "Number of Observations", ylab = "Relative Efficiency", main = "tau-optimal Designs for exp-log Models")
lines(x = c(2:15), y = 1 / elt_true, col = "blue")
legend("bottomright", c("True design", "Misspecified design"), col = c("blue", "red"), lty = c(1, 1))

###
#elh_comp <- list()

elh_app1 <- exp_log_hoptimal_approx_pso(4, el_par1, 0.15)
elh_app2 <- exp_log_hoptimal_approx_pso(4, el_par2, 0.15)
elh_appcomp <- (elh_app1$val / exp_log_hoptimal(d = elh_app2$par[1:7], loc = el_par1)) ^ (1/4)
elh_appcomp


for(j in 13:15){
  elh_ex1 <- exp_log_hoptimal_pso(j, el_par1, psoinfo_el, 0.15)
  elh_ex2 <- exp_log_hoptimal_pso(j, el_par2, psoinfo_el, 0.15)
  elh_comp[[j]] <- (elh_ex1$val / exp_log_hoptimal(d = elh_ex2$par, loc = el_par1)) ^ (1/j)
  #print(elh_comp[[j]])
}
print(elh_comp)

plot(x = c(4:15), y = sapply(4:15, function(x) elh_comp[[x]]), type = "l", col = "blue", ylim = c(0.74, 1),  
     xlab = "Number of Observations", ylab = "Relative Efficiency", main = "h-optimal Designs for exp-log Models")
abline(h = elh_appcomp, col = "red")
legend(x = 9.5, y = 0.85, c("Exact designs", "Approximate design"), col = c("blue", "red"), lty = c(1, 1))

poly_infmat <- function(x, k){
  n = length(x)
  xi <- lapply(1:(n/k), function(j) matrix(x[(k*j-(k-1)):(k*j)]))
  M <- lapply(1:(n/k), function(j) xi[[j]] %*% t(xi[[j]]) * k/n)
  Reduce("+", M)
}

poly_infmat(x = c(0,0,1,2,3,5), k = 3)

poly_dop <- function(x, loc){
  inf_mat <- poly_infmat(x, loc)
  -det(inf_mat)
}

poly_dop_pso <- function(Np, Nf, l, u, psoinfo){
  ub <- rep(u, Np*Nf)
  lb <- rep(l, Np*Nf)
  pso_res <- globpso(poly_dop, lower = lb, upper = ub, PSO_INFO = psoinfo, verbose = F, loc = Nf)
  pso_res$par <- sapply(1:Np, function(x) pso_res$par[(Nf*x-(Nf-1)):(Nf*x)])
  pso_res
}

poly_aop <- function(x, loc){
  inf_mat <- poly_infmat(x, loc)
  
  if (rcond(inf_mat) < 2.220446e-16){
    res = 9e+10
  } 
  else{
    M_inv <- solve(inf_mat)
    res <- sum(diag(M_inv))
  }
  
  res
}

poly_aop_pso <- function(Np, Nf, l, u, psoinfo){
  ub <- rep(u, Np*Nf)
  lb <- rep(l, Np*Nf)
  pso_res <- globpso(poly_aop, lower = lb, upper = ub, PSO_INFO = psoinfo, verbose = F, loc = Nf)
  pso_res$par <- sapply(1:Np, function(x) pso_res$par[(Nf*x-(Nf-1)):(Nf*x)])
  pso_res
}

psoinfo <- psoinfo_setting(nSwarms = 256, Iters = 5000)

poly_dop_pso(Np = 8, Nf = 4, l = -1, u = 1, psoinfo = psoinfo)$par

poly_aop_pso(Np = 4, Nf = 3, l = -1, u = 1, psoinfo = psoinfo)$par


hb_hbase <- hb_hoptimal_approx(input = c(0, 0.02, 0.04, 0.359, 0.5), loc = hb_par)
for(j in c(4:20)){
  if (j == 11) psoinfo_hb <- psoinfo_setting(nSwarms = 512, Iters = 10000)
  if (j == 16) psoinfo_hb <- psoinfo_setting(nSwarms = 1024, Iters = 30000)
  hb_hoptiaml_res <- hb_hoptimal_pso_rep(nRep = 5, nPoints = j, parms = hb_par, psoinfo = psoinfo_hb, upper = 0.15)
  print(j)
  print(hb_hoptiaml_res$result)
  print("==============================")
}

hbh_eff <- c(0.8913, 0.9435, 0.9635, 0.9787, 0.9817, 0.9881, 0.9906, 0.9914, 0.9922, 0.9937, 0.9945, 0.9957, 0.9960, 
             0.9967, 0.9969, 0.9972, 0.9973)
plot(x = c(4:20), y = hbh_eff, type = "l", xlab = "N (Number of Observations)", ylab = "Efficiency")
points(x = c(4:20), y = hbh_eff, pch = 19, cex =)

library(matlib)
library(MASS)

hbd_get <- function(x, parms, M_inv){
  c1 <- parms[1]
  tau <- parms[2]
  b0 <- parms[3]
  b1 <- parms[4]
  
  #print(M)
  fx <- hb_f(x, c1, tau, b0, b1)
  
  t(fx) %*% M_inv %*% fx
}

hbd_dose <- c(0, 0.02, 0.04, 0.0991)
hb_grid <- seq(0, 0.15, by = 0.0001)
hbd_mat_list <- lapply(c(0, 0.02, 0.04, 0.0991), function(x) 1/4 * hb_mat(x, 170, 0.04, 1.46, 40))
hbd_M <- Reduce("+", hbd_mat_list)
diag(hbd_M) <- diag(hbd_M) + 1e-10
hbd_M_inv <- solve(hbd_M)

hbd_get_val <- sapply(hb_grid, function(x) hbd_get(x, hb_par, hbd_M_inv))
get_df <- data.frame(dose = hb_grid, val = hbd_get_val)
get2 <- data.frame(dose = hbd_dose, val = hbd_get_val[(hbd_dose / 0.0001)+1])
hbd_get_val[(hbd_dose / 0.0001)+1]
ggplot(data = get_df, aes(x = dose, y = val)) + 
  geom_line(stat = "identity") + 
  geom_point(data = get2, aes(x = dose, y = val), shape = 1, size = 10)
plot(x = hb_grid, y = hbd_get_val, type = "l")
abline(v = c(0, 0.02, 0.04, 0.0991))
abline(h = 4)

hbtau_get <- function(x, parms, M_inv){
  c1 <- parms[1]
  tau <- parms[2]
  b0 <- parms[3]
  b1 <- parms[4]
  
  fx <- hb_f(x, c1, tau, b0, b1)
  b <- matrix(c(0, 1, 0, 0))
  
  #print(t(b) %*% M_inv %*% b)
  (t(fx) %*% M_inv %*% b)^2 - t(b) %*% M_inv %*% b
}

hbtau_get1 <- function(x, parms){
  c1 <- parms[1]
  tau <- parms[2]
  b0 <- parms[3]
  b1 <- parms[4]
  
  fx <- hb_f(x, c1, tau, b0, b1)
  q <- matrix(c(0, -2/(tau * c1), (1 + exp(b0))^2/exp(b0), 0.1))
  
  abs(t(q) %*% fx)
}

hbtaud <- c(0, 0.04)
hbtauw <- c(0.5, 0.5)

hbtaud <- c(0, 0.0399, 0.04, 0.0401)
hbtauw <- c(0.5, 0.1, 0.3, 0.1)

hbtaud <- c(-0.0001, 0, 0.0001, 0.04)
hbtauw <- c(0.1,0.3,0.1,0.5)
hbtau_mat_list <- lapply(1:4, function(x) hbtauw[x] * hb_mat(hbtaud[x], 170, 0.04, 1.46, 40))
hbtau_M <- Reduce("+", hbtau_mat_list)
diag(hbtau_M) <- diag(hbtau_M) + 1e-10
hbtau_M_inv <- solve(hbtau_M)

hbtau_get_val <- sapply(hb_grid, function(x) hbtau_get(x, hb_par, hbtau_M_inv))
plot(x = hb_grid, y = hbtau_get_val, type = "l")
abline(v = c(0, 0.04))
abline(h = 0)

hb_hoptimal_approx(input = c(0, 0.02, 0.04, 0.359, 0.5), loc = c(170, 0.04, 1.46, 40))

hb_hoptimal(d = c(0, 0.022, 0.022, 0.04), loc = c(170, 0.04, 1.46, 40))

(hb_hoptimal_approx(input = c(0, 0.02, 0.04, 0.359, 0.5), loc = c(170, 0.04, 1.46, 40)) / 
    hb_hoptimal(d = c(0, 0.022, 0.022, 0.04), loc = c(170, 0.04, 1.46, 40))) ^ (1/4)

hbh_get <- function(x, parms, M_inv){
  c1 <- parms[1]
  tau <- parms[2]
  b0 <- parms[3]
  b1 <- parms[4]
  
  fx <- hb_f(x, c1, tau, b0, b1)
  h <- matrix(c(-tau, -c1, 0, 0))
  
  (t(fx) %*% M_inv %*% h)^2 - t(h) %*% M_inv %*% h
}


hbhd <- c(0, 0.02, 0.04, 0.08, 0.125)
hbhw <- c(0.359, 0.5, 0.1399, 0.001, 0.001)
# hbhd <- c(0, 0.02, 0.0399, 0.0401)
# hbhw <- c(0.359, 0.5, 0.0705, 0.0705)
hbh_mat_list <- lapply(1:5, function(x) hbhw[x] * hb_mat(hbhd[x], 170, 0.04, 1.86, 40))
hbh_M <- Reduce("+", hbh_mat_list)
diag(hbh_M) <- diag(hbh_M) + 1e-10
hbh_M_inv <- solve(hbh_M)

hbh_get_val <- sapply(hb_grid, function(x) hbh_get(x, c(170, 0.04, 1.86, 40), hbh_M_inv))
plot(x = hb_grid, y = hbh_get_val, type = "l")
abline(v = hbhd)
abline(h = 0)

hb_hoptimal_approx(input = c(0, 0.019, 0.021, 0.04, 0.359, 0.25, 0.25), loc = c(170, 0.04, 1.86, 40))
hb_hoptimal_approx(input = c(0, 0.02, 0.04, 0.359, 0.5), loc = c(170, 0.04, 1.86, 40))

eld_get <- function(x, parms, M_inv){
  c0 = parms[1]
  c1 = parms[2]
  b0 = parms[3]
  b1 = parms[4]
  
  fx <- exp_log_f(x, c0, c1, b0, b1)
  t(fx) %*% M_inv %*% fx
}

el_grid <- seq(0, 0.15, by = 0.0001)
eld_mat_list <- lapply(c(0, 0.0109, 0.0558, 0.1051), function(x) 1/4 * exp_log_mat(x, 0.15, 89, 3.2, 41))
eld_M <- Reduce("+", eld_mat_list)
eld_M_inv <- solve(eld_M)

eld_get_val <- sapply(el_grid, function(x) eld_get(x, el_par, eld_M_inv))
plot(x = el_grid, y = eld_get_val, type = "l")
abline(v = c(0, 0.0109, 0.0558, 0.1051))
abline(h = 4)


# exp_log_b <- function(d, c0, c1, b0, b1){
#   bn <- c0 * c1 * exp(-c1*d) - ((b1 * exp(b0 - b1*d)) / (1 + exp(b0 - b1*d))^2)
#   
#   h1 <- (exp(-c1*d) - 1) / bn
#   h2 <- (-c0 * d * exp(-c1*d)) / bn
#   h3 <- ((-exp(b0 - b1*d) / (1 + exp(b0 - b1*d))^2) + -exp(b0) / (1 + exp(b0))^2) / bn
#   h4 <- ((d * exp(b0 - b1*d)) / (1 + exp(b0 - b1*d))^2) / bn
# 
#   matrix(c(h1, h2, h3, h4))
# }

eltau_get <- function(x, parms, M_inv){
  c0 <- parms[1]
  c1 <- parms[2]
  b0 <- parms[3]
  b1 <- parms[4]
  tau <- uniroot(tau_func, c(1e-10, 0.15), tol = 1e-10, 
                    c0 = 0.15, c1 = 89, b0 = 3.2, b1 = 41)$root
  
  fx <- exp_log_f(x, c0, c1, b0, b1)
  b <- exp_log_b(tau, c0, c1, b0, b1)
  
  (t(fx) %*% M_inv %*% b)^2 - t(b) %*% M_inv %*% b
}

#exp_log_tauoptimal_pso(2, el_par, psoinfo, 0.15)$par

tauval <- uniroot(tau_func, c(1e-10, 0.15), tol = 1e-10, 
        c0 = 0.15, c1 = 89, b0 = 3.2, b1 = 41)$root
eltaud <- c(0, tauval)
eltauw <- c(0.5, 0.5)
eltau_mat_list <- lapply(1:2, function(x) eltauw[x] * exp_log_mat(eltaud[x], 0.15, 89, 3.2, 41))
eltau_M <- Reduce("+", eltau_mat_list)
diag(eltau_M) <- diag(eltau_M) + 1e-5
eltau_M_inv <- solve(eltau_M)

eltau_get_val <- sapply(el_grid, function(x) eltau_get(x, el_par, eltau_M_inv))
plot(x = el_grid, y = eltau_get_val, type = "l")
abline(v = c(0, tauval))
abline(h = 0)

elh_get <- function(x, parms, M_inv){
  c0 <- parms[1]
  c1 <- parms[2]
  b0 <- parms[3]
  b1 <- parms[4]
  
  fx <- exp_log_f(x, c0, c1, b0, b1)
  h <- exp_log_h(x, c0, c1, b0, b1)
  
  (t(fx) %*% M_inv %*% h)^2 - t(h) %*% M_inv %*% h
}

elhd <- c(0, 0.0108, 0.0526, 0.1187)
elhw <- c(0.371, 0.501, 0.087, 0.041)
elh_mat_list <- lapply(1:4, function(x) elhw[x] * exp_log_mat(elhd[x], 0.15, 89, 3.2, 41))
elh_M <- Reduce("+", elh_mat_list)
elh_M_inv <- solve(elh_M)

elh_get_val <- sapply(el_grid, function(x) elh_get(x, c(0.15, 89, 3.2, 41), elh_M_inv))
plot(x = el_grid, y = elh_get_val, type = "l")
abline(v = elhd)
abline(h = 0)

library(OptimalDesign)
library(gurobi)

hb_grid <- seq(0, 0.15, by = 0.00025)
F_HB <- sapply(hb_grid, function(x) hb_f(x, 170, 0.04, 1.46, 40)) |> t()

HB_IP4 <- od_MISOCP(F_HB, b3 = 4, crit="D", bin = FALSE, t.max = 1000, gap = 0)
HB_IP4$t.act #194.60
HB_IP4_od <- data.frame(dose = hb_grid, replications = HB_IP4$w.best) |> subset(replications != 0)
HB_IP4_od

hbd_base <- hb_doptimal(d = c(0, 0.02, 0.04, 0.0991), hb_par)
(hb_doptimal(d = c(0, 0.051, 0.102, 0.150), hb_par)/hbd_base)^(1/4)
(hb_doptimal(d = c(0, 0.0198, 0.04, 0.0993), hb_par) / hbd_base)^(1/4)

HB_IP10 <- od_MISOCP(F_HB, b3 = 10, crit="D", bin = FALSE, t.max = 6000, gap = 0)
HB_IP10$t.act
HB_IP10_od <- data.frame(dose = hb_grid, replications = HB_IP10$w.best) |> subset(replications != 0)
HB_IP10_od
(hb_doptimal(d = c(rep(0,3), rep(0.051,3), rep(0.102,3), rep(0.15,1)), hb_par)/hb_doptimal(d = c(0, 0.02, 0.04, 0.0991), hb_par))^(1/10)

hbh_func <- c(-0.04, -170, 0, 0)
HBh_IP4 <- od_MISOCP(F_HB, b3 = 4, crit="c", h = hbh_func, bin = FALSE, t.max = Inf, gap = 0)
HBh_IP4$t.act
HBh_IP4_od <- data.frame(dose = hb_grid, replications = HBh_IP4$w.best) |> subset(replications != 0)
HBh_IP4_od
# (0, 0.022, 0.04); (1, 2, 1)
hbh_base <- hb_hoptimal_approx(input = c(0, 0.02, 0.04, 0.359, 0.5), loc = hb_par)
(hbh_base / hb_hoptimal(d = c(0, 0.021, 0.021, 0.039), loc = hb_par)) ^ (1/4)

HBh_IP10 <- od_MISOCP(F_HB, b3 = 10, crit="c", h = hbh_func, bin = FALSE, t.max = Inf, gap = 0)
HBh_IP10$t.act
HBh_IP10_od <- data.frame(dose = hb_grid, replications = HBh_IP10$w.best) |> subset(replications != 0)
HBh_IP10_od
(hb_hbase / hb_hoptimal(d = c(rep(0, 4), rep(0.018, 5), 0.039), loc = hb_par)) ^ (1/10)

eld_base <- exp_log_doptimal(d = c(0, 0.0109, 0.0558, 0.1051), loc = el_par)
EL_grid <- seq(0, 0.15, by = 0.003)
F_EL <- sapply(EL_grid, function(x) exp_log_f(x, 0.15, 89, 3.2, 41)) |> t()
EL_IP4 <- od_MISOCP(F_EL, b3 = 4, crit="D", bin = FALSE, t.max = 10000, gap = 0)
EL_IP4$t.act
EL_IP4_od <- data.frame(dose = EL_grid, replications = EL_IP4$w.best) |> subset(replications != 0)
EL_IP4_od
( exp_log_doptimal(d = c(0, 0.009, 0.057, 0.105), loc = el_par)/eld_base) ^ (1/4)

EL_IP10 <- od_MISOCP(F_EL, b3 = 10, crit="D", bin = FALSE, t.max = Inf, gap = 0)
EL_IP10$t.act
EL_IP10_od <- data.frame(dose = EL_grid, replications = EL_IP10$w.best) |> subset(replications != 0)
EL_IP10_od
(exp_log_doptimal(d = c(rep(0,4), rep(0.012,2), rep(0.057,2), rep(0.1065,2)), loc = el_par)/eld_base) ^ (1/10)

elh_func <- c(exp_log_h(tauval, 0.15, 89, 3.2, 41))
ELh_IP4 <- od_MISOCP(F_EL, b3 = 4, crit="c", h = elh_func, bin = FALSE, t.max = Inf, gap = 0)
ELh_IP4$t.act
ELh_IP4_od <- data.frame(dose = EL_grid, replications = ELh_IP4$w.best) |> subset(replications != 0)
ELh_IP4_od
elh_base <- exp_log_hoptimal_approx(input = c(0, 0.0108, 0.0526, 0.1187, 0.371, 0.501, 0.087), loc = el_par)
(elh_base / exp_log_hoptimal(d = c(0, 0.012, 0.06, 0.126), loc = el_par))^(1/4)

ELh_IP10 <- od_MISOCP(F_EL, b3 = 10, crit="c", h = elh_func, bin = FALSE, t.max = Inf, gap = 0)
ELh_IP10$t.act
ELh_IP10_od <- data.frame(dose = EL_grid, replications = ELh_IP10$w.best) |> subset(replications != 0)
ELh_IP10_od
(elh_base / exp_log_hoptimal(d = c(rep(0,3), rep(0.012,5), 0.057, 0.123), loc = el_par))^(1/10)

Flogis <- Fx_glm(~x1, c(2, 1), glm.model = "bin-logit", lower = -5, upper = 5, n.levels = 51)
logis_IP <- od_MISOCP(Flogis, b3 = 4, crit = "D", bin = FALSE, t.max = Inf, gap = 0)
logis_IP$t.act
Fx_lin <- Fx_cube(~x1, lower = -5, upper = 5, n.levels = 51)
logis_od <- data.frame(dose = Fx_lin[,2], replications = logis_IP$w.best) |> subset(replications != 0)
logis_od
log_base <- logistic_doptimal(d = c(-3.5434, -0.4566), loc = c(2,1))
(logistic_doptimal(d = c(-3.5, -3.5, -0.4, -0.4), loc = c(2,1)) / log_base) ^ (1/4)

logis2_IP <- od_MISOCP(Flogis, b3 = 10, crit = "D", bin = FALSE, t.max = Inf, gap = 0)
logis2_IP$t.act
logis2_od <- data.frame(dose = Fx_lin[,2], replications = logis2_IP$w.best) |> subset(replications != 0)
logis2_od
(logistic_doptimal(d = c(rep(-3.6,5), -0.6, rep(-0.4,4)), loc = c(2,1)) / log_base) ^ (1/2)

Fqlogis <- Fx_glm(~x1+I(x1^2), c(3, 0, -1), glm.model = "bin-logit", lower = -5, upper = 5, n.levels = 51)
qlogis_IP <- od_MISOCP(Fqlogis, b3 = 4, crit = "D", bin = FALSE, t.max = 3600, gap = 0)
qlogis_IP$t.act
qlogis_od <- data.frame(dose = Fx_lin[,2], replications = qlogis_IP$w.best) |> subset(replications != 0)
qlogis_od
qlog_base <- qlogistic_doptimal_approx(input = c(-2.061, -1.324, 1.324, 2.061, 0.297, 0.203, 0.203), loc = c(3, 0, -1))
(qlogistic_doptimal(d = c(-2.0, -1.2, 1.4, 2.0), loc = c(3, 0, -1)) / qlog_base) ^ (1/3)

qlogis_IP2 <- od_MISOCP(Fqlogis, b3 = 10, crit = "D", bin = FALSE, t.max = 3600, gap = 0)
qlogis_IP2$t.act
qlogis_od2 <- data.frame(dose = Fx_lin[,2], replications = qlogis_IP2$w.best) |> subset(replications != 0)
qlogis_od2
(qlogistic_doptimal(d = c(rep(-2.0,3), rep(-1.2,2), rep(1.2,2), rep(2.0,3)), loc = c(3, 0, -1)) / qlog_base) ^ (1/3)

Fqlogis2 <- Fx_glm(~x1+I(x1^2), c(-3, 0, -1), glm.model = "bin-logit", lower = -5, upper = 5, n.levels =51)
qlogis2_IP <- od_MISOCP(Fqlogis2, b3 = 4, crit = "D", bin = FALSE, t.max = 10800, gap = 0)
qlogis2_IP$t.act
qlogis2_od <- data.frame(dose = Fx_lin[,2], replications = qlogis2_IP$w.best) |> subset(replications != 0)
qlogis2_od

qlogis2_IP2 <- od_MISOCP(Fqlogis2, b3 = 10, crit = "D", bin = FALSE, t.max = 10800, gap = 0)
qlogis2_IP2$t.act
qlogis2_od2 <- data.frame(dose = Fx_lin[,2], replications = qlogis2_IP2$w.best) |> subset(replications != 0)
qlogis2_od2

qpso <- qlogistic_pso(10, c(3, 0, -1), psoinfo, 5)
qpso$par
qpso$cputime

Fclogis <- Fx_glm(~x1+I(x1^2)+I(x1^3), c(-1, 3, 2, -1), glm.model = "bin-logit", lower = -5, upper = 5, n.levels = 51)
clogis_IP <- od_MISOCP(Fclogis, b3 = 5, crit = "D", bin = FALSE, t.max = 10800, gap = 0, echo = FALSE)
clogis_IP$t.act
clogis_od <- data.frame(dose = Fx_lin[,2], replications = clogis_IP$w.best) |> subset(replications != 0)
clogis_od

clogis_IP2 <- od_MISOCP(Fclogis, b3 = 10, crit = "D", bin = FALSE, t.max = 10800, gap = 0)
clogis_IP2$t.act
clogis_od2 <- data.frame(dose = Fx_lin[,2], replications = clogis_IP2$w.best) |> subset(replications != 0)
clogis_od2


Fclogis2 <- Fx_glm(~x1+I(x1^2)+I(x1^3), c(-3, 0, 0, -1), glm.model = "bin-logit", lower = -5, upper = 5, n.levels = 51)
clogis2_IP <- od_MISOCP(Fclogis2, b3 = 4, crit = "D", bin = FALSE, t.max = 10800, gap = 0, echo = FALSE)
clogis2_IP$t.act
clogis2_od <- data.frame(dose = Fx_lin[,2], replications = clogis2_IP$w.best) |> subset(replications != 0)
clogis2_od
clog_base <- clogistic_doptimal(d = c(-1.6909, -1.2436, 0.2388, 1.1354), loc = c(-3, 0, 0, -1))
(clogistic_doptimal(d = c(-1.6, -1.2, 0.2, 1.2), loc = c(-3, 0, 0, -1)) / clog_base) ^ (1/4)

clogis2_IP2 <- od_MISOCP(Fclogis2, b3 = 10, crit = "D", bin = FALSE, t.max = 10800, gap = NULL)
clogis2_IP2$t.act
clogis2_od2 <- data.frame(dose = Fx_lin[,2], replications = clogis2_IP2$w.best) |> subset(replications != 0)
clogis2_od2
(clogistic_doptimal(d = c(-1.8,rep(-1.6,2),rep(-1.2,2),rep(0.2,2),1.0,rep(1.2,2)), loc = c(-3, 0, 0, -1)) / clog_base) ^ (1/4)

mm_f <- function(t, a, b){
  c(t/(b+t), -a/(b+t^2))
}

mm_corr <- function(t, sigma, lambda){
  
}


clog_pso <- clogistic_pso(10, c(-3, 0, 0, -1), psoinfo, 5)
clog_pso$par
clog_pso$cputime

clog_base <- clogistic_doptimal(d = c(-1.6909, -1.2436, 0.2388, 1.1354), loc = c(-3, 0, 0, -1))
(clogistic_doptimal(d = c(rep(-1.6909, 3), rep(-1.2436, 2), rep(0.2388, 3), rep(1.1354, 2)), loc = c(-3, 0, 0, -1))/clog_base)^(1/10)

eld_val <- c()
eld_val <- c(eld_val, 
             exp_log_doptimal(d = c(rep(0, 5), rep(0.0109, 5), rep(0.0558, 5), rep(0.1051, 5)), loc = c(0.15, 89, 3.2, 41)))
eld_base <- eld_val[1]
eld_eff <- sapply(4:20, function(x) (eld_val[x-3]/eld_base)^(1/x))
plot(x = c(4:20), y = eld_eff, type = "l", xlab = "N (Number of Observations)", ylab = "Effiency")
points(x = c(4:20), y = eld_eff, pch = 19, cex = .75)

hbd_val <- c()
hbd_val <- c(hbd_val, 
             hb_doptimal(d = c(rep(0, 5), rep(0.02, 5), rep(0.04, 5), rep(0.0991, 5)), loc = c(170, 0.04, 1.46, 40)))
hbd_base <- hbd_val[1]
hbd_eff <- sapply(4:20, function(x) (hbd_val[x-3]/hbd_base)^(1/x))
plot(x = c(4:20), y = hbd_eff, type = "l", xlab = "N (Number of Observations)", ylab = "Effiency")
points(x = c(4:20), y = hbd_eff, pch = 19, cex = .75)

hbtau_val <- c()
hbtau_val <- c(hbtau_val, 
             hb_tauoptimal(d = c(rep(0, 10), rep(0.04, 10)), loc = c(170, 0.04, 1.46, 40)))
hbtau_base <- hbtau_val[1]
hbtau_eff <- sapply(2:20, function(x) (hbtau_base/hbtau_val[x-1])^(1/x))
plot(x = c(2:20), y = hbtau_eff, type = "l", xlab = "N (Number of Observations)", ylab = "Effiency")
points(x = c(2:20), y = hbtau_eff, pch = 19, cex = .75)

eltau_val <- c()
eltau_val <- c(eltau_val, 
               exp_log_tauoptimal(d = c(rep(0, 1), rep(tauval, 1)), loc = c(0.15, 89, 3.2, 41)))
hbtau_base <- hbtau_val[1]
hbtau_eff <- sapply(2:20, function(x) (hbtau_base/hbtau_val[x-1])^(1/x))
plot(x = c(2:20), y = hbtau_eff, type = "l", xlab = "N (Number of Observations)", ylab = "Effiency")
points(x = c(2:20), y = hbtau_eff, pch = 19, cex = .75)

psoinfo <- psoinfo_setting(128, 2000)
hbh_res <- list()
for(j in 7:20){
  if (j == 11) psoinfo <- psoinfo_setting(512, 5000)
  if (j == 16) psoinfo <- psoinfo_setting(512, 20000)
  pso_results <- list()
  pso_results <- lapply(1:10, function(x) pso_results[[x]] <- hb_hoptimal_pso(j, hb_par, psoinfo, 0.15))
  
  hb_list <- list()
  hb_list$val <- sapply(1:10, function(x) pso_results[[x]]$val) 
  hb_list$design_points <- sapply(1:10, function(x) pso_results[[x]]$par) 
  hbh_res[[j]] <- hb_list
  print(hb_list)
}
hbh_result <- hb_hoptimal_pso(10, hb_par, psoinfo, 0.15)
hbh_result$par

hb_hoptimal(d = c(rep(0, 7), rep(0.0204, 9), rep(0.0469, 3), rep(0.1266, 1)), loc = hb_par)
hb_hoptimal_approx(input = c(0, 0.02, 0.04, 0.359, 0.5), loc = hb_par)
(hb_hoptimal_approx(input = c(0, 0.02, 0.04, 0.359, 0.5), loc = hb_par) / 
    hb_hoptimal(d = c(rep(0, 2), rep(0.0219, 4), rep(0.04, 2)), loc = hb_par)) ^(1/8)

psoinfo <- getPSOInfo(nSwarm = 64, maxIter = 1000, tol = 1e-6, freeRun = 0.8)
hbd_pso_res <- hb_doptimal_pso(4, hb_par, psoinfo, 0.15)
hbd_pso_res$par
hbd_pso_res$val
hbd_pso_res$cputime

library(metaheuristicOpt)
library(nlodm)
lb <- c(rep(0,3), rep(1e-4,3))
ub <- c(rep(0.15,3), rep(1,3))
range <- cbind(lb, ub)
test_res <- globpso(objFunc = hb_tauoptimal_approx, lower = lb, upper = ub, 
                    PSO_INFO = psoinfo, verbose = F, loc = hb_par)
test_res$par

test_res <- nlodm(grad_fun = hb_grad, obj = 'c_e', theta = hb_par, bound = 0.15, pts = 3, algorithm = 'PSO', 
                  swarm = 100, iter = 1000, seed = 1, exact = F, c = c(0,1,0,0))
test_res$plot

test_res_h <- nlodm(grad_fun = hb_grad, obj = 'c', theta = hb_par, bound = 0.15, pts = 4, algorithm = 'DE', 
                  swarm = 100, iter = 1000, seed = 1, exact = F, c = c(-0.04, -170, 0, 0))
test_res_h$design

eltau_test_res <- nlodm(grad_fun = exp_log_grad, obj = 'c_e', theta = hb_par, bound = 0.15, pts = 3, 
                        algorithm = 'PSO', swarm = 100, iter = 1000, seed = 1, exact = F, c = c(0,1,0,0))
test_res$plot

hb_tau_d <- test_res$par[1:3]
hb_tau_w <- test_res$par[4:6]
hbtau_mat_list <- lapply(1:3, function(x) hb_tau_w[x] * hb_mat(hb_tau_d[x], 170, 0.04, 1.46, 40))
hbtau_M <- Reduce("+", hbtau_mat_list)
diag(hbtau_M) <- diag(hbtau_M) + 1e-10
hbtau_M_inv <- solve(hbtau_M)

hbtau_get_val <- sapply(hb_grid, function(x) hbtau_get(x, hb_par, hbtau_M_inv))
plot(x = hb_grid, y = hbtau_get_val, type = "l")
abline(v = c(0, 0.04))
abline(h = 0)

psoinfo <- psoinfo_setting(128, 3000)
lb <- c(rep(0,3), rep(1e-4,3))
ub <- c(rep(0.15,3), rep(1,3))
test_res <- globpso(objFunc = hb_hoptimal_approx, lower = lb, upper = ub, 
                    PSO_INFO = psoinfo, verbose = F, loc = hb_par)
test_res$par

hbhd <- test_res$par[1:3] |> round(4)
hbhw <- test_res$par[4:6]
hbh_mat_list <- lapply(1:4, function(x) test_res_h$design$w[x] * hb_mat(test_res_h$design$x[x], 170, 0.04, 1.86, 40))
hbh_M <- Reduce("+", hbh_mat_list)
diag(hbh_M) <- diag(hbh_M) + 1e-10
hbh_M_inv <- solve(hbh_M)

hbh_get_val <- sapply(hb_grid, function(x) hbh_get(x, c(170, 0.04, 1.86, 40), hbh_M_inv))
plot(x = hb_grid, y = hbh_get_val, type = "l")
abline(v = hbhd)
abline(h = 0)

test_res <- metaOpt(test_res, optimType = "MIN", algorithm = "DE", numVar = 6, rangeVar = range)

el_hf <- c(exp_log_h(0, 0.15, 89, 3.2, 41))
el_test_res <- nlodm(grad_fun = exp_log_grad, obj = 'c', theta = el_par, bound = 0.15, pts = 4, algorithm = 'DE', 
                  swarm = 100, iter = 1000, seed = 1, exact = F, c = el_hf)
el_test_res$plot

tauval
el_bf <- c(exp_log_h(tauval, 0.15, 89, 3.2, 41))
el_tau_test_res <- nlodm(grad_fun = exp_log_grad, obj = 'c_e', theta = el_par, bound = 0.15, pts = 2, algorithm = 'PSO', 
                     swarm = 200, iter = 1000, seed = 123, exact = F, c = el_bf)
el_tau_test_res$design
el_tau_test_res$plot

exp_log_tauoptimal_approx <- function(input, loc){
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
  diag(M) <- diag(M) + 1e-5
  if (rcond(M) < 2.220446e-16){
    res = pen
  }
  else {
    M_inv <- solve(M)
    res = t(b) %*% M_inv %*% b
  }
  
  res
}

exp_log_tau_pso()

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

library(DEoptim)

hbtau_de <- DEoptim(fn = hb_tauoptimal_approx, lower = c(rep(0,6)), upper = c(rep(0.15, 3), rep(1, 3)), 
                    control = DEoptim.control(NP = 100, itermax = 1000, trace = F), loc = hb_par)
hbtau_de$optim$bestmem |> round(4)
hb_tau_d <- hbtau_de$optim$bestmem[1:3]
hb_tau_w <- hbtau_de$optim$bestmem[4:6]
hbtau_mat_list <- lapply(1:3, function(x) hb_tau_w[x] * hb_mat(hb_tau_d[x], 170, 0.04, 1.46, 40))
hbtau_M <- Reduce("+", hbtau_mat_list)
diag(hbtau_M) <- diag(hbtau_M) + 1e-5
hbtau_M_inv <- solve(hbtau_M)
hb_tauoptimal(d = c(0, 0.04), loc = hb_par)
c_vec <- matrix(c(0, 1, 0, 0))
t(c_vec) %*% ginv(hbtau_M) %*% c_vec
v = sqrt(1 / t(c_vec) %*% ginv(hbtau_M) %*% c_vec)
hbtau_elf <- lapply(1:3, function(x) (-1)^x * hb_tau_w[x] * hb_f(hb_tau_d[x], 170, 0.04, 1.46, 40))
hbtau_elf_vec <- Reduce("+", hbtau_elf)


hbtau_get_val <- sapply(hb_grid, function(x) hbtau_get(x, hb_par, hbtau_M_inv))
plot(x = hb_grid, y = hbtau_get_val, type = "l")
abline(v = c(0, 0.04))

hb_tau_d2 <- c(0, 0.04)
hb_tau_w2 <- c(0.5, 0.5)
hbtau_mat_list2 <- lapply(1:2, function(x) hb_tau_w2[x] * hb_mat(hb_tau_d2[x], 170, 0.04, 1.46, 40))
hbtau_mat2 <- Reduce("+", hbtau_mat_list2, )
hbtau_mat2_ginv <- ginv(hbtau_mat2)
sqrt(1/t(c_vec) %*% hbtau_mat2_ginv %*% c_vec)
hb_f(0, 170, 0.04, 1.46, 40) - hb_f(0.04, 170, 0.04, 1.46, 40)


hb_tauoptimal_approx(input = c(hb_tau_d, hb_tau_w), loc = hb_par)
hb_tauoptimal_approx(input = c(0, 0.04, 0.5, 0.5), loc = hb_par)

hbh_de <- DEoptim(fn = hb_hoptimal_approx, lower = c(rep(0,8)), upper = c(rep(0.15, 4), rep(1, 4)), 
                    control = DEoptim.control(NP = 100, itermax = 1000, trace = F), loc = hb_par)
hbh_de$optim$bestmem |> round(4)
hb_h_d <- hbh_de$optim$bestmem[1:4]
hb_h_w <- hbh_de$optim$bestmem[5:8]
hbh_mat_list <- lapply(1:4, function(x) hb_h_w[x] * hb_mat(hb_h_d[x], 170, 0.04, 1.46, 40))
hbh_M <- Reduce("+", hbh_mat_list)
diag(hbh_M) <- diag(hbh_M) + 1e-5
hbh_M_inv <- solve(hbh_M)

hbh_get_val <- sapply(hb_grid, function(x) hbh_get(x, hb_par, hbh_M_inv))
plot(x = hb_grid, y = hbh_get_val, type = "l")
abline(v = c(0, 0.02, 0.04))
abline(h = 0)

hbd_eff_1 <- sapply(1:17, function(x) (hbd_eff[x]^(x+3))^(1/4))
plot(x = c(4:20), y = hbd_eff_1, xlab = "Number of Observations", ylab = "D-Efficiency", 
     type = "l")
points(x = c(4:20), y = hbd_eff_1, pch = 19, cex = .75)

eld_eff_1 <- sapply(1:17, function(x) (eld_eff[x]^(x+3))^(1/4))
plot(x = c(4:20), y = eld_eff_1, xlab = "Number of Observations", ylab = "D-Efficiency", 
     type = "l")
points(x = c(4:20), y = eld_eff_1, pch = 19, cex = .75)

hbtau_val <- c(hbtau_val, hb_tauoptimal(d = c(rep(0,10), rep(0.04,10)), loc = hb_par))
hbtau_eff <- hbtau_val[1] / hbtau_val
plot(x = c(2:20), y = hbtau_eff, xlab = "Number of Observations", ylab = "tau-Efficiency", 
     type = "l")
points(x = c(2:20), y = hbtau_eff, pch = 19, cex = .75)

eltau_val <- c(eltau_val, exp_log_tauoptimal(d = c(rep(0,10), rep(tauval,10)), loc = c(el_par, tauval)))
eltau_eff <- eltau_val[1] / eltau_val
plot(x = c(2:20), y = eltau_eff, xlab = "Number of Observations", ylab = "tau-Efficiency", 
     type = "l")
points(x = c(2:20), y = eltau_eff, pch = 19, cex = .75)

hbh_val <- c(hbh_val, hb_hoptimal(d = c(rep(0,7), rep(0.0204,9), rep(0.0469,3), rep(0.1266,1)), 
                                  loc = hb_par))
hbh_eff <- hbh_base / hbh_val
plot(x = c(4:20), y = hbh_eff, xlab = "Number of Observations", ylab = "h-Efficiency", 
     type = "l")
points(x = c(4:20), y = hbh_eff, pch = 19, cex = .75)

elh_val <- c(elh_val, exp_log_hoptimal(d = c(rep(0,7), rep(0.0113,10), rep(0.0532,2), rep(0.1191,1)), 
                                       loc = el_par))
elh_val <- elh_base / elh_val
elh_eff <- elh_base / elh_val
plot(x = c(4:20), y = elh_eff, xlab = "Number of Observations", ylab = "h-Efficiency", 
     type = "l")
points(x = c(4:20), y = elh_eff, pch = 19, cex = .75)

qlog1_val <- c(qlog1_val, qlogistic_doptimal(d = c(rep(-1.238,6),rep(0,7),rep(1.238,7)), 
                                             loc = qlog_par2))
qlog1_eff <- sapply(3:20, function(x) (qlog1_val[x-2] / qlog1_val[1])^(1/3))
plot(x = c(3:20), y = qlog1_eff, xlab = "Number of Observations", ylab = "D-Efficiency", 
     type = "l")
points(x = c(3:20), y = qlog1_eff, pch = 19, cex = .75)

qlog2_val <- c(qlog2_val, qlogistic_doptimal(d = c(rep(-2.0592,6),rep(-1.3187,4),rep(1.3187,4),rep(2.0592,6)), 
                                             loc = qlog_par1))
qlog_base <- qlogistic_doptimal_approx(input = c(-2.061,-1.324,1.324,2.061,0.297,0.203,0.203), loc = qlog_par1)
qlog2_eff <- sapply(4:20, function(x) (qlog2_val[x-3] / qlog_base)^(1/3))
plot(x = c(4:20), y = qlog2_eff, xlab = "Numbers of Observations", ylab = "D-Efficiency", 
     type = "l")
points(x = c(4:20), y = qlog2_eff, pch = 19, cex = .75)

clog1_val <- c(clog1_val, clogistic_doptimal(d = c(rep(-1.6909,5), rep(-1.2436,5), 
                                                   rep(0.2388,5), rep(1.1354,5)), 
                                             loc = clog_par1))
clog1_eff <- sapply(4:12, function(x) (clog1_val[x-3] / clog1_val[1])^(1/4))
plot(x = c(4:12), y = clog1_eff, xlab = "Numbers of Observations", ylab = "D-Efficiency", 
     type = "l")
points(x = c(4:12), y = clog1_eff, pch = 19, cex = .75)

clog2_val <- c(clog2_val, clogistic_doptimal(d = c(rep(-1.1328,3), rep(-0.4231,2),
                                                   rep(0.1190,2), rep(2.9604,2), 
                                                   rep(3.1779,3)), 
                                             loc = clog_par2))
clog_approx <- clogistic_approx_pso(clog_par2, 5)
clog2_base <- clog_approx$val * -1
clog2_eff <- sapply(5:12, function(x) (clog2_val[x-4] / clog2_base)^(1/4))
plot(x = c(5:12), y = clog2_eff, xlab = "Numbers of Observations", ylab = "D-Efficiency", 
     type = "l")
points(x = c(5:12), y = clog2_eff, pch = 19, cex = .75)

elh_base / exp_log_hoptimal(d = c(0, rep(0.013,2), 0.0591, 0.1242), loc = el_par)
elh_base / exp_log_hoptimal(d = c(0, rep(0.012,2), 0.0581, 0.1232), loc = el_par)
elh_base / exp_log_hoptimal(d = c(0, rep(0.011,2), 0.0571, 0.1222), loc = el_par)
elh_base / exp_log_hoptimal(d = c(0, rep(0.010,2), 0.0561, 0.1212), loc = el_par)
elh_base / exp_log_hoptimal(d = c(0, rep(0.014,2), 0.0601, 0.1252), loc = el_par)
elh_base / exp_log_hoptimal(d = c(0, rep(0.015,2), 0.0611, 0.1262), loc = el_par)
elh_base / exp_log_hoptimal(d = c(0, rep(0.016,2), 0.0621, 0.1272), loc = el_par)

elh_base / exp_log_hoptimal(d = c(rep(0,3), rep(0.0116,5), 0.0561, 0.1242), loc = el_par)
elh_base / exp_log_hoptimal(d = c(rep(0,3), rep(0.0106,5), 0.0551, 0.1232), loc = el_par)
elh_base / exp_log_hoptimal(d = c(rep(0,3), rep(0.0096,5), 0.0541, 0.1222), loc = el_par)
elh_base / exp_log_hoptimal(d = c(rep(0,3), rep(0.0086,5), 0.0531, 0.1212), loc = el_par)
elh_base / exp_log_hoptimal(d = c(rep(0,3), rep(0.0126,5), 0.0571, 0.1252), loc = el_par)
elh_base / exp_log_hoptimal(d = c(rep(0,3), rep(0.0136,5), 0.0581, 0.1262), loc = el_par)
elh_base / exp_log_hoptimal(d = c(rep(0,3), rep(0.0146,5), 0.0591, 0.1272), loc = el_par)

elh_base / exp_log_hoptimal(d = c(rep(0,2), rep(0.0130,1), rep(0.0591,1), rep(0.1242,1)), loc = el_par)
elh_base / exp_log_hoptimal(d = c(rep(0,1), rep(0.0130,1), rep(0.0591,2), rep(0.1242,1)), loc = el_par)
elh_base / exp_log_hoptimal(d = c(rep(0,1), rep(0.0130,1), rep(0.0591,1), rep(0.1242,2)), loc = el_par)
elh_base / exp_log_hoptimal(d = c(rep(0,2), rep(0.0130,2), rep(0.0591,1), rep(0.1242,0)), loc = el_par)

elh_base / exp_log_hoptimal(d = c(rep(0,4), rep(0.0116,4), rep(0.0561,1), rep(0.1242,1)), loc = el_par)
elh_base / exp_log_hoptimal(d = c(rep(0,3), rep(0.0116,4), rep(0.0561,2), rep(0.1242,1)), loc = el_par)
elh_base / exp_log_hoptimal(d = c(rep(0,3), rep(0.0116,4), rep(0.0561,1), rep(0.1242,2)), loc = el_par)
elh_base / exp_log_hoptimal(d = c(rep(0,2), rep(0.0116,6), rep(0.0561,1), rep(0.1242,1)), loc = el_par)

hbh_base / hb_hoptimal(d = c(0, rep(0.021,2), 0.039), loc = hb_par)
hbh_base / hb_hoptimal(d = c(rep(0,4), rep(0.018,5), rep(0.039,1)), loc = hb_par)

elh_base / exp_log_hoptimal(d = c(0,0.012,0.06,0.126), loc = el_par)
elh_base / exp_log_hoptimal(d = c(rep(0,3),rep(0.012,5),0.057,0.123), loc = el_par)

(exp_log_doptimal(d = c(0,0.009,0.057,0.105), loc = el_par) / eld_base)^(1/4)
(exp_log_doptimal(d = c(rep(0,4), rep(0.012,2), rep(0.057,2), rep(0.1065,2)), loc = el_par) / eld_base)^(1/4)

(logistic_doptimal(d = c(-3.6,-3.4,-0.4), loc = c(2,-1)) / 
  logistic_doptimal(d = c(-3.5434, -0.4566), loc = c(2,-1)))^(1/2)
(logistic_doptimal(d = c(rep(-3.6,5), -0.6, rep(-0.4,4)), loc = c(2,-1)) / 
    logistic_doptimal(d = c(-3.5434, -0.4566), loc = c(2,-1)))^(1/2)

library(AlgDesign)
efficient.rounding(c(0.371, 0.501, 0.087, 0.041), 10)
efficient.rounding(c(0.297, 0.203, 0.203, 0.297), 9)
efficient.rounding(c(0.245, 0.108, 0.167, 0.238, 0.242), 10)
qlogistic_doptimal_round <- function(d_sup, loc){
  a <- loc[1]
  b1 <- loc[2]
  b2 <- loc[3]
  d <- c(rep(d_sup[1], 3), rep(d_sup[2], 3), rep(d_sup[3], 3), rep(d_sup[4], 3))
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

lb <- rep(-1 * 5, 4)
ub <- rep(5, 4)

qlog_round_pso_res <- globpso(objFunc = qlogistic_doptimal_round, lower = lb, upper = ub, 
                   PSO_INFO = psoinfo, verbose = F, loc = qlog_par)
qlog_round_pso_res$par |> round(4)
(qlog_round_pso_res$val / qlog_base) ^ (1/3)


clogistic_doptimal_round <- function(d_sup, loc){
  a <- loc[1]
  b1 <- loc[2]
  b2 <- loc[3]
  b3 <- loc[4]
  d <- c(rep(d_sup[1], 3), rep(d_sup[2], 1), rep(d_sup[3], 2), rep(d_sup[4], 2), rep(d_sup[5], 2))
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
lb <- rep(-1 * 5, 5)
ub <- rep(5, 5)

clog_round_pso_res <- globpso(objFunc = clogistic_doptimal_round, lower = lb, upper = ub, 
                              PSO_INFO = psoinfo, verbose = F, loc = clog_par2)
clog_round_pso_res$par |> round(4)
(clog_round_pso_res$val / clog2_base) ^ (1/4)

clog_pso_test <- clogistic_pso(10, clog_par2, psoinfo, 5)

hbtau_fit <- function(regressor, parm){
  n <- nrow(regressor)
  regressor[,1] <- regressor[,1] + rnorm(n, 0, 0.1)
  reg <- as.data.frame(regressor)
  colnames(reg) <- c("y", "c1", "tau", "b0", "b1")
  
  tau_lm <- lm(y~0+., reg)
  parm - tau_lm$coef
}

hbtau_sim <- function(d, parm, nsim){
  c1 = parm[1]
  tau = parm[2]
  b0 = parm[3]
  b1 = parm[4]
  n = length(d)
  
  regressor <- sapply(d, function(x) c(hunt_bowman(x, c1, tau, b0, b1), hb_f(x, c1, tau, b0, b1))) |> t()
  fitted <- sapply(1:nsim, function(x) hbtau_fit(regressor, parm)) |> t()

  fitted
}

hbtau_sim_res <- hbtau_sim(d = c(rep(0,2), rep(0.04,2)), parm = hb_par, nsim = 10000)
var2 <- var(hbtau_sim_res[,2])
#hist(hbtau_sim_res[,3])
crit2 <- hb_tauoptimal(d = c(rep(0,2), rep(0.04,2)), loc = hb_par)

hbtau_sim_res2 <- hbtau_sim(d = c(rep(0,1), rep(0.04,3)), parm = hb_par, nsim = 10000)
var1 <- var(hbtau_sim_res2[,2])
#hist(hbtau_sim_res[,3])
crit1 <- hb_tauoptimal(d = c(rep(0,1), rep(0.04,3)), loc = hb_par)

hbtau_sim_res2 <- hbtau_sim(d = c(rep(0,3), rep(0.04,1)), parm = hb_par, nsim = 10000)
var3 <- var(hbtau_sim_res2[,2])
#hist(hbtau_sim_res[,3])
crit3 <- hb_tauoptimal(d = c(rep(0,3), rep(0.04,1)), loc = hb_par)

plot(x = c(1:3), y = c(var1, var2, var3)*1000, type = "l", ylim = c(0.05, 0.3))
lines(x = c(1:3), y = c(crit1, crit2, crit3))

c(var1, var2, var3)*1000 - c(crit1, crit2, crit3)

hbtau_sim_res <- hbtau_sim(d = c(rep(0,5), rep(0.04,5)), parm = hb_par, nsim = 10000)
var2 <- var(hbtau_sim_res[,2])
#hist(hbtau_sim_res[,3])
crit2 <- hb_tauoptimal(d = c(rep(0,5), rep(0.04,5)), loc = hb_par)

hbtau_sim_res2 <- hbtau_sim(d = c(rep(0,4), rep(0.04,6)), parm = hb_par, nsim = 10000)
var1 <- var(hbtau_sim_res2[,2])
#hist(hbtau_sim_res[,3])
crit1 <- hb_tauoptimal(d = c(rep(0,4), rep(0.04,6)), loc = hb_par)

hbtau_sim_res4 <- hbtau_sim(d = c(rep(0,3), rep(0.04,7)), parm = hb_par, nsim = 10000)
var4 <- var(hbtau_sim_res4[,2])
#hist(hbtau_sim_res[,3])
crit4 <- hb_tauoptimal(d = c(rep(0,3), rep(0.04,7)), loc = hb_par)

hbtau_sim_res2 <- hbtau_sim(d = c(rep(0,6), rep(0.04,4)), parm = hb_par, nsim = 10000)
var3 <- var(hbtau_sim_res2[,2])
#hist(hbtau_sim_res[,3])
crit3 <- hb_tauoptimal(d = c(rep(0,6), rep(0.04,4)), loc = hb_par)

hbtau_sim_res5 <- hbtau_sim(d = c(rep(0,7), rep(0.04,3)), parm = hb_par, nsim = 10000)
var5 <- var(hbtau_sim_res5[,2])
#hist(hbtau_sim_res[,3])
crit5 <- hb_tauoptimal(d = c(rep(0,7), rep(0.04,3)), loc = hb_par)

plot(x = c(1:5), y = c(var4, var1, var2, var3, var5)*850, type = "l", ylim = c(0.07, 0.11), 
     xaxt = "n", xlab = "Number of Replicates on each support (0, 0.04)", ylab = "")
axis(1, at = 1:5, labels = c("(3,7)", "(4,6)", "(5,5)", "(6,4)", "(7,3)"))
lines(x = c(1:5), y = c(crit4, crit1, crit2, crit3, crit5), lty = "dashed")
legend("bottomright", c("Simulated Variance", "Criterion Value"), lty = c(1, 3))

c(var1, var2, var3)*1000 - c(crit1, crit2, crit3)

hbh_fit <- function(regressor, parm){
  n <- nrow(regressor)
  regressor[,1] <- regressor[,1] + rnorm(n, 0, 0.1)
  reg <- as.data.frame(regressor)
  colnames(reg) <- c("y", "c1", "tau", "b0", "b1")
  
  h_lm <- lm(y~0+., reg)
  coeff <- parm - h_lm$coef
  -coeff[1] * coeff[2]
}

hbh_sim <- function(d, parm, nsim){
  c1 = parm[1]
  tau = parm[2]
  b0 = parm[3]
  b1 = parm[4]
  n = length(d)
  
  regressor <- sapply(d, function(x) c(hunt_bowman(x, c1, tau, b0, b1), hb_f(x, c1, tau, b0, b1))) |> t()
  fitted <- sapply(1:nsim, function(x) hbh_fit(regressor, parm))
  
  fitted
}

var2 <- hbh_sim(d = c(rep(0,2), rep(0.02, 4), rep(0.04,2)), parm = hb_par, nsim = 10000) |> var()
crit2 <- hb_hoptimal(d = c(rep(0,2), rep(0.02, 4), rep(0.04,2)), loc = hb_par)

var1 <- hbh_sim(d = c(rep(0,1), rep(0.02, 5), rep(0.04,2)), parm = hb_par, nsim = 10000) |> var()
crit1 <- hb_hoptimal(d = c(rep(0,1), rep(0.02, 5), rep(0.04,2)), loc = hb_par)

var3 <- hbh_sim(d = c(rep(0,3), rep(0.02, 3), rep(0.04,2)), parm = hb_par, nsim = 10000) |> var()
crit3 <- hb_hoptimal(d = c(rep(0,3), rep(0.02, 3), rep(0.04,2)), loc = hb_par)

plot(x = c(1:3), y = c(var1, var2, var3)*600, type = "l", ylim = c(33000, 50000))
lines(x = c(1:3), y = c(crit1, crit2, crit3))

eltau_fit <- function(regressor, parm){
  n <- nrow(regressor)
  regressor[,1] <- regressor[,1] + rnorm(n, 0, 0.1)
  reg <- as.data.frame(regressor)
  colnames(reg) <- c("y", "c0", "c1", "b0", "b1")
  #print(reg)
  
  tau_lm <- lm(y~0+., reg)
  coeff <- parm - tau_lm$coef
  coeff[2] = parm[2]
  coeff[4] = parm[4]
  print(coeff)
  tau = tryCatch({
    uniroot(tau_func, c(1e-10, 0.15), tol = 1e-10, 
          c0 = coeff[1], c1 = coeff[2], b0 = coeff[3], b1 = coeff[4], extendInt = "yes")$root
  })
  print(tau)
}

eltau_sim <- function(d, parm, nsim){
  c0 = parm[1]
  c1 = parm[2]
  b0 = parm[3]
  b1 = parm[4]
  n = length(d)
  
  regressor <- sapply(d, function(x) c(exp_log(x, c0, c1, b0, b1), hb_f(x, c0, c1, b0, b1))) |> t()
  fitted <- sapply(1:nsim, function(x) eltau_fit(regressor, parm)) |> t()
  
  fitted
}

eltau_sim(d = c(rep(0,5), rep(tauval, 5)), parm = el_par, nsim = 10000)

qpsoinfo <- getPSOInfo(nSwarm = 64, maxIter = 500, psoType = "quantum")
hbh_lb <- rep(0, 4)
hbh_ub <- rep(0.15, 4)
qpso_hbh <- globpso(hb_hoptimal, hbh_lb, hbh_ub, PSO_INFO = qpsoinfo, loc = hb_par)
qpso_hbh$par |> sort() |> round(4)
qpso_hbh$val

hbh_qpso <- list()

for(j in 4:20){
  if(j <= 10){
    qpsoinfo <- getPSOInfo(nSwarm = 128, maxIter = 2000, psoType = "quantum")
  } 
  else if(j >= 11 & j <= 15){
    qpsoinfo <- getPSOInfo(nSwarm = 1024, maxIter = 10000, psoType = "quantum")
  }
  else if(j >= 16){
    qpsoinfo <- getPSOInfo(nSwarm = 1024, maxIter = 30000, psoType = "quantum")
  }
  
  hbh_lb <- rep(0, j)
  hbh_ub <- rep(0.15, j)
  
  qpso_list <- list()
  qpso_list <- lapply(1:10, function(x) qpso_list[[x]] <- globpso(hb_hoptimal, hbh_lb, hbh_ub, PSO_INFO = qpsoinfo, loc = hb_par))
  
  hbh_qpso[[j]] <- list(val = c(), design_points = c())
  hbh_qpso[[j]]$val <- sapply(1:10, function(x) qpso_list[[x]]$val) 
  hbh_qpso[[j]]$design_points <- sapply(1:10, function(x) qpso_list[[x]]$par |> sort() |> round(4)) 
  print(hbh_qpso[[j]])
}

hbh_qpso_res <- list(eff = c(), dp = list())
for(j in 4:20){
  idx <- which.min(hbh_qpso[[j]]$val)
  hbh_qpso_res$eff <- c(hbh_qpso_res$eff, (hbh_base/hbh_qpso[[j]]$val[idx]))
  hbh_qpso_res$dp[[j]] <- hbh_qpso[[j]]$design_points[,idx]
}
hbh_qpso_res$eff
hbh_qpso_res$dp[[6]]

hb_hoptimal_approx(input = c(0, 0.02, 0.04, 0.359, 0.500, 0.141), loc = hb_par)
hbh_base/hb_hoptimal(d = c(rep(0, 1), rep(0.0219, 2), rep(0.04, 1)), loc = hb_par)
hbh_base/hb_hoptimal(d = c(rep(0, 1), rep(0.0218, 1), rep(0.0219, 1), rep(0.04, 1)), loc = hb_par)

hbh_lb <- rep(0, 10)
hbh_ub <- rep(0.15, 10)
qpso_test <- globpso(exp_log_hoptimal, hbh_lb, hbh_ub, PSO_INFO = qpsoinfo, loc = el_par)
(elh_base / qpso_test$val)
qpso_test$par |> round(4) |> sort()

qpsoinfo <- getPSOInfo(nSwarm = 64, maxIter = 1000, psoType = "quantum")

log_lb <- rep(-5, 4)
log_ub <- rep(5, 4)
qpso_test <- globpso(clogistic_doptimal, log_lb, log_ub, PSO_INFO = qpsoinfo, loc = clog_par)
(qpso_test$val / clog_base)^(1/4)
qpso_test$par |> round(4) |> sort()

qpsoinfo <- getPSOInfo(nSwarm = 128, maxIter = 2000, psoType = "quantum")
sapply(1:10, function(x) hbh_base/globpso(hb_hoptimal, hbh_lb, hbh_ub, PSO_INFO = qpsoinfo, loc = hb_par)$val) |> sort() |> round(4)
sapply(1:10, function(x) elh_base/globpso(exp_log_hoptimal, hbh_lb, hbh_ub, PSO_INFO = qpsoinfo, loc = el_par)$val) |> sort() |> round(4)

log_lb <- rep(-5, 10)
log_ub <- rep(5, 10)
qpsoinfo <- getPSOInfo(nSwarm = 128, maxIter = 2000, psoType = "quantum")
sapply(1:10, function(x) (globpso(logistic_doptimal, log_lb, log_ub, PSO_INFO = qpsoinfo, loc = log_par)$val/log_base)^(1/2)) |> sort() |> round(4)
sapply(1:10, function(x) (globpso(qlogistic_doptimal, log_lb, log_ub, PSO_INFO = qpsoinfo, loc = qlog_par)$val/qlog_base)^(1/3)) |> sort() |> round(4)
sapply(1:10, function(x) (globpso(clogistic_doptimal, log_lb, log_ub, PSO_INFO = qpsoinfo, loc = clog_par)$val/clog_base)^(1/4)) |> sort() |> round(4)


library(DEoptim)

hbh_de <- DEoptim(fn = hb_hoptimal, lower = c(rep(0,10)), upper = c(rep(0.15, 10)), 
                    control = DEoptim.control(NP = 128, itermax = 2000, trace = F), loc = hb_par)
hbh_de$optim$bestmem |> round(4) |> sort()
hbh_de$optim$bestval
hbh_base/hb_hoptimal(d = c(rep(0, 4), rep(0.0195, 5), rep(0.04, 1)), loc = hb_par)

hbh_de$optim

hbh_de <- list()

for(j in 19:20){
  if(j <= 10){
    deinfo <- DEoptim.control(NP = 128, itermax = 2000, trace = F)
  } 
  else if(j >= 11 & j <= 15){
    deinfo <- DEoptim.control(NP = 1024, itermax = 10000, trace = F)
  }
  else if(j >= 16){
    deinfo <- DEoptim.control(NP = 1024, itermax = 30000, trace = F)
  }
  
  hbh_lb <- rep(0, j)
  hbh_ub <- rep(0.15, j)
  
  de_list <- list()
  de_list <- lapply(1:3, function(x) de_list[[x]] <- DEoptim(fn = hb_hoptimal, lower = hbh_lb, upper = hbh_ub, 
                                                              control =deinfo, loc = hb_par))
  
  hbh_de[[j]] <- list(val = c(), design_points = c())  
  hbh_de[[j]]$val <- sapply(1:3, function(x) de_list[[x]]$optim$bestval) 
  hbh_de[[j]]$design_points <- sapply(1:3, function(x) de_list[[x]]$optim$bestmem |> sort() |> round(4)) 
  print(hbh_de[[j]])
}

deinfo <- DEoptim.control(NP = 256, itermax = 2000, trace = F)
hbh_lb <- rep(0, 8)
hbh_ub <- rep(0.15, 8)
hbh_detest <- DEoptim(fn = hb_hoptimal, lower = hbh_lb, upper = hbh_ub, 
        control =deinfo, loc = hb_par)
hbh_detest$optim$bestmem |> sort() |> round(4)

st_time <- Sys.time()
ed_time <- Sys.time()

hbh_de_res <- list(eff = c(), dp = list())
for(j in 4:20){
  idx <- which.min(hbh_de[[j]]$val)
  hbh_de_res$eff <- c(hbh_de_res$eff, (hbh_base/hbh_de[[j]]$val[idx]))
  hbh_de_res$dp[[j]] <- hbh_de[[j]]$design_points[,idx]
}
hbh_de_res$eff
hbh_qpso_res$dp[[10]]

hbh_base/hb_hoptimal(d = c(rep(0, 6), rep(0.020, 9), rep(0.0482, 2), 0.1289), loc = hb_par)

hbh_lb <- rep(-5, 4)
hbh_ub <- rep(5, 4)
deinfo <- DEoptim.control(NP = 128, itermax = 2000, trace = F)
st_time <- Sys.time()
detest <- DEoptim(fn = clogistic_doptimal, lower = hbh_lb, upper = hbh_ub, 
                      control = deinfo, loc = clog_par)
ed_time <- Sys.time()
ed_time - st_time
(detest$optim$bestval/clog_base)^(1/4)
detest$optim$bestmem |> sort() |> round(4)

eodh <- function(model, criterion, Nrun, ){
  
}

library(AlgDesign)

hbh_effr <- sapply(4:20, function(x) efficient.rounding(c(0.359, 0.5, 0.141), x))
hbh_effr_val <- sapply(1:17, function(x) hb_hoptimal(d = c(rep(0, hbh_effr[1,x]), rep(0.02, hbh_effr[2,x]), rep(0.04, hbh_effr[3,x])), loc = hb_par))
hb_hoptimal(d = c(0, 0.02, 0.02, 0.04), loc = hb_par)

round(((hbh_base/hbh_effr_val)/hbh_eff), 4)


elh_effr <- sapply(4:20, function(x) efficient.rounding(c(0.371, 0.501, 0.087, 0.041), x))
elh_effr_val <- 
  sapply(1:17, function(x) exp_log_hoptimal(d = c(rep(0, elh_effr[1,x]), rep(0.0108, elh_effr[2,x]), rep(0.0526, elh_effr[3,x]), rep(0.1187, elh_effr[4,x])), loc = el_par))

(elh_val / elh_effr_val) |> round(4)


qlog_effr <- sapply(4:20, function(x) efficient.rounding(c(0.297, 0.203, 0.203, 0.297), x, random = F))
qlog_effr_val <- 
  sapply(1:17, function(x) qlogistic_doptimal(d = c(rep(-2.061, qlog_effr[1,x]), rep(-1.324, qlog_effr[2,x]), rep(1.324, qlog_effr[3,x]), rep(2.061, qlog_effr[4,x])), loc = qlog_par))

(qlog_effr_val / qlog2_val)^(1/3) |> round(4)

clog_effr <- sapply(5:20, function(x) efficient.rounding(c(0.2454, 0.1082, 0.1669, 0.2377, 0.2418), x, random = F))
clog_effr_val <- 
  sapply(1:16, function(x) clogistic_doptimal(d = c(rep(-1.1703, clog_effr[1,x]), rep(-0.3497, clog_effr[2,x]), rep(0.0005, clog_effr[3,x]), rep(2.9576, clog_effr[4,x]), rep(3.1874, clog_effr[5,x])), loc = clog_par))
((clog_effr_val/clog2_base) / clog1_eff)^(1/4) |> round(4)

psoinfo <- psoinfo_setting(256, 2000)
el_par <- el_params(c0 = 0.15, c1 = 89, b0 = 3.2, b1 = 41)
exp_log_plot(parms = el_par, upper_bound = 0.15)

#c0: control the effect size of hormesis, no effect to the design when changing c0 and fixing the other parameters
# c0 = 0.15
el_par <- el_params(c0 = 0.15, c1 = 89, b0 = 3.2, b1 = 41)
exp_log_plot(parms = el_par, upper_bound = 0.15)
el_c0_1 <- hormesis_exact(model = "ExpLog", criterion = "h", nPoints = 10, parms = el_par, 
                          psoinfo = psoinfo, upper = 0.15, lower = 0)
el_c0_1$exact

# c0 = 0.05
el_par <- el_params(c0 = 0.05, c1 = 89, b0 = 3.2, b1 = 41)
exp_log_plot(parms = el_par, upper_bound = 0.15)
el_c0_2 <- hormesis_exact(model = "ExpLog", criterion = "h", nPoints = 10, parms = el_par, 
                          psoinfo = psoinfo, upper = 0.15, lower = 0)
el_c0_2$exact

# c0 = 0.25
el_par <- el_params(c0 = 0.25, c1 = 89, b0 = 3.2, b1 = 41)
exp_log_plot(parms = el_par, upper_bound = 0.15)
el_c0_3 <- hormesis_exact(model = "ExpLog", criterion = "h", nPoints = 10, parms = el_par, 
                          psoinfo = psoinfo, upper = 0.15, lower = 0)
el_c0_3$exact


# c1
# c1 = 89
el_par <- el_params(c0 = 0.15, c1 = 89, b0 = 3.2, b1 = 41)
exp_log_plot(parms = el_par, upper_bound = 0.15)
el_c1_1 <- hormesis_exact(model = "ExpLog", criterion = "h", nPoints = 10, parms = el_par, 
                          psoinfo = psoinfo, upper = 0.15, lower = 0)
el_c1_1$exact

# c1 = 60.
el_par <- el_params(c0 = 0.15, c1 = 60, b0 = 3.2, b1 = 41)
exp_log_plot(parms = el_par, upper_bound = 0.15)
el_c1_2 <- hormesis_exact(model = "ExpLog", criterion = "h", nPoints = 10, parms = el_par, 
                          psoinfo = psoinfo, upper = 0.15, lower = 0)
el_c1_2$exact

# c1 = 120
el_par <- el_params(c0 = 0.15, c1 = 120, b0 = 3.2, b1 = 41)
exp_log_plot(parms = el_par, upper_bound = 0.15)
el_c1_3 <- hormesis_exact(model = "ExpLog", criterion = "h", nPoints = 10, parms = el_par, 
                          psoinfo = psoinfo, upper = 0.15, lower = 0)
el_c1_3$exact


# b0
# b0 = 3.2
el_par <- el_params(c0 = 0.15, c1 = 89, b0 = 3.2, b1 = 41)
exp_log_plot(parms = el_par, upper_bound = 0.15)
el_b0_1 <- hormesis_exact(model = "ExpLog", criterion = "h", nPoints = 10, parms = el_par, 
                          psoinfo = psoinfo, upper = 0.15, lower = 0)
el_b0_1$exact

# b0 = 2.0
el_par <- el_params(c0 = 0.15, c1 = 89, b0 = 2, b1 = 41)
exp_log_plot(parms = el_par, upper_bound = 0.15)
el_b0_2 <- hormesis_exact(model = "ExpLog", criterion = "h", nPoints = 10, parms = el_par, 
                          psoinfo = psoinfo, upper = 0.15, lower = 0)
el_b0_2$exact
el_b0_2$val

# b0 = 4.5
el_par <- el_params(c0 = 0.15, c1 = 89, b0 = 4.5, b1 = 41)
exp_log_plot(parms = el_par, upper_bound = 0.15)
el_b0_3 <- hormesis_exact(model = "ExpLog", criterion = "h", nPoints = 10, parms = el_par, 
                          psoinfo = psoinfo, upper = 0.15, lower = 0)
el_b0_3$exact

# b1
# b1 = 41
el_par <- el_params(c0 = 0.15, c1 = 89, b0 = 3.2, b1 = 41)
exp_log_plot(parms = el_par, upper_bound = 0.15)
el_b1_1 <- hormesis_exact(model = "ExpLog", criterion = "h", nPoints = 10, parms = el_par, 
                          psoinfo = psoinfo, upper = 0.15, lower = 0)
el_b1_1$exact

# b1 = 60.
el_par <- el_params(c0 = 0.15, c1 = 89, b0 = 3.2, b1 = 60)
exp_log_plot(parms = el_par, upper_bound = 0.15)
el_b1_2 <- hormesis_exact(model = "ExpLog", criterion = "h", nPoints = 10, parms = el_par, 
                          psoinfo = psoinfo, upper = 0.15, lower = 0)
el_b1_2$exact

# b1 = 20
el_par <- el_params(c0 = 0.15, c1 = 89, b0 = 3.2, b1 = 20)
exp_log_plot(parms = el_par, upper_bound = 0.15)
el_b1_3 <- hormesis_exact(model = "ExpLog", criterion = "h", nPoints = 10, parms = el_par, 
                          psoinfo = psoinfo, upper = 0.15, lower = 0)
el_b1_3$exact




# Hunt-Bowman origin D-optimal design
hb_par <- hb_parms(c1 = 170, tau = 0.04, b0 = 1.46, b1 = 40)
hunt_bowman_plot(hb_par, 0.15)
hb_d_origin <- hormesis_exact(model = "HuntBowman", criterion = "D", nPoints = 10, parms = hb_par, 
                              psoinfo = psoinfo, upper = 0.15, lower = 0)
hb_d_origin$exact

hb_par <- hb_parms(c1 = 6.8, tau = 0.2, b0 = 1.46, b1 = 8) # c1: /25, tau: *5, b1: /5
hunt_bowman_plot(hb_par, 0.75)

hb_d_5t <- hormesis_exact(model = "HuntBowman", criterion = "D", nPoints = 10, parms = hb_par, 
                          psoinfo = psoinfo, upper = 0.75, lower = 0)
hb_d_5t$exact
hb_d_5t$val

hb_par <- hb_parms(c1 = 1.7, tau = 0.4, b0 = 1.46, b1 = 4) # c1: /100, tau: *10, b1: /10
hunt_bowman_plot(hb_par, 1.5)
hb_d_10t <- hormesis_exact(model = "HuntBowman", criterion = "D", nPoints = 10, parms = hb_par, 
                          psoinfo = psoinfo, upper = 1.5, lower = 0)
hb_d_10t$exact

hb_par <- hb_parms(c1 = 0.425, tau = 0.8, b0 = 1.46, b1 = 2) # c1: /400, tau: *20, b1: /20
hunt_bowman_plot(hb_par, 3)
hb_d_20t <- hormesis_exact(model = "HuntBowman", criterion = "D", nPoints = 10, parms = hb_par, 
                          psoinfo = psoinfo, upper = 3, lower = 0)
hb_d_20t$exact

psoinfo <- psoinfo_setting(256, 2000)
hb_dop_5 <- globpso(hb_doptimal_approx, lower = rep(0, 10), upper = c(rep(0.15, 5), rep(1, 5)), loc = hb_par, PSO_INFO = psoinfo)
hb_dop_5$par |> round(4)
hb_dop_5$val


hb_par <- hb_parms(c1 = 170, tau = 0.04, b0 = 1.46, b1 = 40)
hb_hop <- hormesis_exact(model = "HuntBowman", criterion = "h", nPoints = 10, parms = hb_par, upper = 0.15, lower = 0, psoinfo = psoinfo)
hb_hop$exact

hb_par <- hb_parms(c1 = 1.7, tau = 0.4, b0 = 1.46, b1 = 4)
hb_hop_10 <- hormesis_exact(model = "HuntBowman", criterion = "h", nPoints = 10, parms = hb_par, upper = 1.5, lower = 0, psoinfo = psoinfo)
hb_hop_10$exact

# Hunt-Bowman D-optimal c1 decrease
hb_par <- hb_parms(c1 = 85, tau = 0.04, b0 = 1.46, b1 = 40)
hunt_bowman_plot(hb_par, 0.15)
hb_d_c1d <- hormesis_exact(model = "HuntBowman", criterion = "D", nPoints = 10, parms = hb_par, 
                              psoinfo = psoinfo, upper = 0.15, lower = 0)
hb_d_c1d$exact

# Hunt-Bowman D-optimal c1 increase
hb_par <- hb_parms(c1 = 255, tau = 0.04, b0 = 1.46, b1 = 40)
hunt_bowman_plot(hb_par, 0.15)
hb_d_c1i <- hormesis_exact(model = "HuntBowman", criterion = "D", nPoints = 10, parms = hb_par, 
                           psoinfo = psoinfo, upper = 0.15, lower = 0)
hb_d_c1i$exact

# Hunt-Bowman D-optimal tau decrease
hb_par <- hb_parms(c1 = 170, tau = 0.02, b0 = 1.46, b1 = 40)
hunt_bowman_plot(hb_par, 0.15)
hb_d_taud <- hormesis_exact(model = "HuntBowman", criterion = "D", nPoints = 10, parms = hb_par, 
                              psoinfo = psoinfo, upper = 0.15, lower = 0)
hb_d_taud$exact

# Hunt-Bowman D-optimal tau increase
hb_par <- hb_parms(c1 = 170, tau = 0.06, b0 = 1.46, b1 = 40)
hunt_bowman_plot(hb_par, 0.15)
hb_d_taui <- hormesis_exact(model = "HuntBowman", criterion = "D", nPoints = 10, parms = hb_par, 
                              psoinfo = psoinfo, upper = 0.15, lower = 0)
hb_d_taui$exact

# Hunt-Bowman D-optimal b0 decrease
hb_par <- hb_parms(c1 = 170, tau = 0.04, b0 = 0.73, b1 = 40)
hunt_bowman_plot(hb_par, 0.15)
hb_d_b0d <- hormesis_exact(model = "HuntBowman", criterion = "D", nPoints = 10, parms = hb_par, 
                              psoinfo = psoinfo, upper = 0.15, lower = 0)
hb_d_b0d$exact

# Hunt-Bowman D-optimal b0 increase
hb_par <- hb_parms(c1 = 170, tau = 0.04, b0 = 2.19, b1 = 40)
hunt_bowman_plot(hb_par, 0.15)
hb_d_b0i <- hormesis_exact(model = "HuntBowman", criterion = "D", nPoints = 10, parms = hb_par, 
                              psoinfo = psoinfo, upper = 0.15, lower = 0)
hb_d_b0i$exact

hbd_b0i_approx <- hormesis_approx(model = "HuntBowman", criterion = "D", parms = hb_par, 
                                  psoinfo = psoinfo, upper = 0.15, lower = 0)
hbd_b0i_approx$approx

# Hunt-Bowman D-optimal b1 decrease
hb_par <- hb_parms(c1 = 170, tau = 0.04, b0 = 1.46, b1 = 20)
hunt_bowman_plot(hb_par, 0.15)
hb_d_b1d <- hormesis_exact(model = "HuntBowman", criterion = "D", nPoints = 10, parms = hb_par, 
                              psoinfo = psoinfo, upper = 0.15, lower = 0)
hb_d_b1d$exact
hb_d_b1d$val

# Hunt-Bowman D-optimal b1 increase
hb_par <- hb_parms(c1 = 170, tau = 0.04, b0 = 1.46, b1 = 60)
hunt_bowman_plot(hb_par, 0.15)
hb_d_b1i <- hormesis_exact(model = "HuntBowman", criterion = "D", nPoints = 10, parms = hb_par, 
                              psoinfo = psoinfo, upper = 0.15, lower = 0)
hb_d_b1i$exact

log_par <- c(2, 1)
logistic_plot(parms = log_par, lb = -5, ub = 5)

log_par <- c(1, 1/2)
logistic_plot(parms = log_par, lb = -10, ub = 10)

log_par <- c(3/4, 5/4)
logistic_plot(parms = log_par, lb = -3, ub = 5)

library(AlgDesign)

psoinfo <- psoinfo_setting(256, 2000)
hb_approx <- hormesis_approx(model = "HuntBowman", criterion = "h", parms = hb_par, 
                             psoinfo = psoinfo, upper = 0.15, lower = 0)
hb_approx$approx
hb_effr <- efficient.rounding(hb_approx$approx$weight, 10)
mu_pso <- sapply(1:length(hb_effr), function(x) rep(hb_approx$approx$support[x], hb_effr[x])) |> unlist()
mvnorm_pso <- mvrnorm(n = 15, mu = mu_pso, Sigma = diag(10) * 0.01)
mvnorm_pso[mvnorm_pso < 0] <- 0
mvnorm_pso[mvnorm_pso > 0.15] <- 0.15
mvnorm_pso
mvnorm_pso <- rbind(mvnorm_pso, mu_pso)

psoinfo <- psoinfo_setting(16, 1000)
hbh_pso_mvnorm <- globpso(objFunc = hb_hoptimal, lower = rep(0, 10), upper = rep(0.15, 10), 
                          init = mvnorm_pso, PSO_INFO = psoinfo, loc = hb_par)
hbh_pso_mvnorm$par |> round(4) |> sort()
hbh_pso_mvnorm$val


hormesis_modified_hbh <- function(model, criterion, parms, psoinfo, upper, lower, nPoints){
  
  psoinfo_approx <- psoinfo_setting(128, 1000)
  approx_lb <- rep(0, 7)
  approx_ub <- c(rep(0.15, 4), rep(1, 3))
  approx_result <- globpso(objFunc = hb_hoptimal_approx, lower = approx_lb, upper = approx_ub, PSO_INFO = psoinfo_approx, 
                           loc = hb_par)
  approx_d <- approx_result$par[1:4]
  approx_w <- approx_result$par[5:7]
  approx_w <- c(approx_w, 1 - sum(approx_w))
  print(round(approx_w, 4))
  approx_val <- approx_result$val
  
  effr <- efficient.rounding(approx_w, nPoints)
  mu_pso <- sapply(1:length(effr), function(x) rep(approx_d[x], effr[x])) |> unlist()
  mvnorm_pso <- mvrnorm(n = 31, mu = mu_pso, Sigma = diag(10) * 0.01)
  mvnorm_pso[mvnorm_pso < 0] <- 0
  mvnorm_pso[mvnorm_pso > 0.15] <- 0.15
  mvnorm_pso
  mvnorm_pso <- rbind(mvnorm_pso, mu_pso)
  
  psoinfo <- psoinfo_setting(32, 1000)
  pso_mvnorm <- globpso(objFunc = hb_hoptimal, lower = rep(0, nPoints), upper = rep(0.15, nPoints), 
                        init = mvnorm_pso, PSO_INFO = psoinfo, loc = hb_par)
  pso_mvnorm$par |> round(4) |> sort()
  pso_mvnorm
  
}

hbh_m_4p <- hormesis_modified_hbh(model = "HuntBowman", criterion = "h", parms = hb_par, psoinfo = psoinfo, 
                                  upper = 0.15, lower = 0, nPoints = 10)
hbh_m_4p$par |> round(4) |> sort()

psoinfo_exact <- psoinfo_setting(32, 500)
psoinfo_exact2 <- psoinfo_setting(64, 1000)
psoinfo_approx <- psoinfo_setting(64, 1000)
hb_par <- hb_parms(c1 = 170, tau = 0.04, b0 = 1.46, b1 = 40)
hb_d_mvnorm <- list()
for(j in 4:20){
  hb_d_mvnorm[[j]] <- hormesis_exact_mvnorm(model = "HuntBowman", criterion = "D", parms = hb_par, 
                                            upper = 0.15, lower = 0, nPoints = j, psoinfo_exact = psoinfo_exact, 
                                            psoinfo_approx = psoinfo_approx)
  #print(hb_d_mvnorm[[j]])
}
hb_d_mvnorm


hb_tau_mvnorm <- list()
for(j in 2:20){
  hb_tau_mvnorm[[j]] <- hormesis_exact_mvnorm(model = "HuntBowman", criterion = "tau", parms = hb_par, 
                                            upper = 0.15, lower = 0, nPoints = j, psoinfo_exact = psoinfo_exact, 
                                            psoinfo_approx = psoinfo_approx)
  #print(hb_tau_mvnorm[[j]])
}
hb_tau_mvnorm

hb_h_mvnorm <- list()
for(j in 4:20){
  hb_h_mvnorm[[j]] <- hormesis_exact_mvnorm(model = "HuntBowman", criterion = "h", parms = hb_par, 
                                              upper = 0.15, lower = 0, nPoints = j, psoinfo_exact = psoinfo_exact, 
                                            psoinfo_approx = psoinfo_approx, nRep = 10)
  #print(hb_h_mvnorm[[j]])
}
for(j in c(10, 17:20)){
  hb_h_mvnorm[[j]] <- hormesis_exact_mvnorm(model = "HuntBowman", criterion = "h", parms = hb_par, 
                                            upper = 0.15, lower = 0, nPoints = j, psoinfo_exact = psoinfo_exact2, 
                                            psoinfo_approx = psoinfo_approx, nRep = 10)
}
hb_h_mvnorm


el_par <- el_parms(c0 = 0.15, c1 = 89, b0 = 3.2, b1 = 41)
el_d_mvnorm <- list()
for(j in 4:20){
  el_d_mvnorm[[j]] <- hormesis_exact_mvnorm(model = "ExpLog", criterion = "D", parms = el_par, 
                                            upper = 0.15, lower = 0, nPoints = j, psoinfo_exact = psoinfo_exact, 
                                            psoinfo_approx = psoinfo_approx)
  #print(el_d_mvnorm[[j]])
}
el_d_mvnorm


el_tau_mvnorm <- list()
for(j in 2:20){
  el_tau_mvnorm[[j]] <- hormesis_exact_mvnorm(model = "ExpLog", criterion = "tau", parms = el_par, 
                                              upper = 0.15, lower = 0, nPoints = j, psoinfo_exact = psoinfo_exact, 
                                              psoinfo_approx = psoinfo_approx)
  #print(el_tau_mvnorm[[j]])
}
el_tau_mvnorm

el_h_mvnorm <- list()
for(j in 4:20){
  el_h_mvnorm[[j]] <- hormesis_exact_mvnorm(model = "ExpLog", criterion = "h", parms = el_par, 
                                            upper = 0.15, lower = 0, nPoints = j, psoinfo_exact = psoinfo_exact, 
                                            psoinfo_approx = psoinfo_approx)
  #print(el_h_mvnorm[[j]])
}
el_h_mvnorm


log_par <- logistic_parms(alpha = 2, beta = 1)
log_mvnorm <- list()
for(j in 2:20){
  log_mvnorm[[j]] <- hormesis_exact_mvnorm(model = "logistic", criterion = "D", parms = log_par, 
                                            upper = 5, lower = -5, nPoints = j, psoinfo_exact = psoinfo_exact, 
                                           psoinfo_approx = psoinfo_approx)
  #print(log_mvnorm[[j]])
}
log_mvnorm

qlog_par1 <- qlogistic_parms(alpha = 3, beta1 = 0, beta2 = -1)
qlog_mvnorm1 <- list()
for(j in 4:20){
  qlog_mvnorm1[[j]] <- hormesis_exact_mvnorm(model = "qlogistic", criterion = "D", parms = qlog_par1, 
                                             upper = 5, lower = -5, nPoints = j, psoinfo_exact = psoinfo_exact2, 
                                             psoinfo_approx = psoinfo_approx)
  #print(qlog_mvnorm1[[j]])
}
qlog_mvnorm1

qlog_par2 <- qlogistic_parms(alpha = -3, beta1 = 0, beta2 = -1)
qlog_mvnorm2 <- list()
for(j in 4:20){
  qlog_mvnorm2[[j]] <- hormesis_exact_mvnorm(model = "qlogistic", criterion = "D", parms = qlog_par2, 
                                             upper = 5, lower = -5, nPoints = j, psoinfo_exact = psoinfo_exact, 
                                             psoinfo_approx = psoinfo_approx)
  #print(qlog_mvnorm2[[j]])
}
qlog_mvnorm2

clog_par1 <- clogistic_parms(alpha = 1, beta1 = 3, beta2 = 2, beta3 = -1)
clog_mvnorm1 <- list()
for(j in 5:20){
  clog_mvnorm1[[j]] <- hormesis_exact_mvnorm(model = "clogistic", criterion = "D", parms = clog_par1, 
                                             upper = 5, lower = -5, nPoints = j, psoinfo_exact = psoinfo_exact2, 
                                             psoinfo_approx = psoinfo_approx, nRep = 1)
  #print(clog_mvnorm1[[j]])
}
clog_mvnorm1

clog_par2 <- clogistic_parms(alpha = -3, beta1 = 0, beta2 = 0, beta3 = -1)
clog_mvnorm2 <- list()
for(j in 5:20){
  clog_mvnorm2[[j]] <- hormesis_exact_mvnorm(model = "clogistic", criterion = "D", parms = clog_par2, 
                                             upper = 5, lower = -5, nPoints = j, psoinfo_exact = psoinfo_exact, 
                                             psoinfo_approx = psoinfo_approx)
  #print(clog_mvnorm2[[j]])
}
clog_mvnorm2

hbh_eff <- c()
hbh_eff <- c(hbh_eff, hb_hbase/hb_hoptimal(d = c(rep(0,7), rep(0.0202, 10), rep(0.04, 3)), loc = hb_par))
round(hbh_eff / 1.0026201, 4)
plot(x = c(4:20), y = hbh_eff / 1.0026201, xlab = "Numbers of Observations", ylab = "h-Efficiency", 
     type = "l")
points(x = c(4:20), y = hbh_eff / 1.0026201, pch = 19, cex = .75)

plot(x = c(2:20), y = eltau_eff, xlab = "Numbers of Observations", ylab = "tau-Efficiency", 
     type = "l")
points(x = c(2:20), y = eltau_eff, pch = 19, cex = .75)

clog1_eff <- sapply(5:20, function(x) clog_mvnorm1[[x]]$efficiency)
plot(x = c(5:20), y = clog1_eff, xlab = "Numbers of Observations", ylab = "D-Efficiency", 
     type = "l")
points(x = c(5:20), y = clog1_eff, pch = 19, cex = .75)

deinfo_exact <- getDEInfo(nPop = 32, maxIter = 500, deType = "rand-to-best-1")
deinfo_approx <- getDEInfo(nPop = 64, maxIter = 1000, deType = "rand-to-best-1")

hormesis_exact_de(model = "HuntBowman", criterion = "D", parms = hb_par, 
                      upper = 0.15, lower = 0, nPoints = 10, deinfo_exact = deinfo_exact, 
                      deinfo_approx = deinfo_approx)
(hb_doptimal(d = c(rep(0, 3), rep(0.02, 2), rep(0.04, 3), rep(0.0989, 2)), loc = hb_par) / hbd_base)^(1/4)

hormesis_exact_de(model = "HuntBowman", criterion = "h", parms = hb_par, 
                  upper = 0.15, lower = 0, nPoints = 10, deinfo_exact = deinfo_exact, 
                  deinfo_approx = deinfo_approx)
hb_hoptimal(d = c(rep(0, 3), rep(0.02, 4), 0.04), loc = hb_par)/hb_hoptimal(d = c(rep(0, 4), 0.0183, 0.0194, 0.0195, 0.0195, 0.0196, 0.04), loc = hb_par)

hormesis_exact_de(model = "ExpLog", criterion = "h", parms = el_par, 
                  upper = 0.15, lower = 0, nPoints = 10, deinfo_exact = deinfo_exact, 
                  deinfo_approx = deinfo_approx)

hormesis_exact_de(model = "logistic", criterion = "D", parms = log_par, 
                  upper = 5, lower = -5, nPoints = 4, deinfo_exact = deinfo_exact, 
                  deinfo_approx = deinfo_approx)

hormesis_exact_de(model = "qlogistic", criterion = "D", parms = qlog_par1, 
                  upper = 5, lower = -5, nPoints = 10, deinfo_exact = deinfo_exact, 
                  deinfo_approx = deinfo_approx)

hormesis_exact_de(model = "clogistic", criterion = "D", parms = clog_par2, 
                  upper = 5, lower = -5, nPoints = 4, deinfo_exact = deinfo_exact, 
                  deinfo_approx = deinfo_approx)
(clogit_doptimal(d = c(-1.6943, -1.6887, -1.2471, -1.2425, 0.2331, 0.2386, 0.2434, 1.1327, 1.1327, 1.1356), loc = clog_par2) / 
    clogit_doptimal(d = c(-1.6909, -1.2436, 0.2388, 1.1354), loc = clog_par2))^(1/4)

hbd_approx_result <- globpso(objFunc = hb_doptimal_approx, lower = rep(0, 7), upper = c(rep(0.15, 4), rep(1, 3)), PSO_INFO = psoinfo_approx, 
                         loc = hb_par, verbose = T)
hbd_approx_result$par |> round(4)
plot(x = c(1:1001), y = hbd_approx_result$history / hbd_base, type = "l", xlab = "Iteration", ylab = "D-efficiency")

hbh_approx_result <- globpso(objFunc = hb_hoptimal_approx, lower = rep(0, 7), upper = c(rep(0.15, 4), rep(1, 3)), PSO_INFO = psoinfo_approx, 
                         loc = hb_par, verbose = T)
hbh_approx_result$par
hbh_base <- hb_hoptimal(d = c(rep(0, 3), rep(0.02, 4), 0.04), loc = hb_par) |> as.numeric()
plot(x = c(1:1001), y = hbh_base / hbh_approx_result$history, type = "l", xlab = "Iteration", ylab = "h-efficiency")

hbd_support <- c(0, 0.02, 0.04, 0.0991)
hbd_5_effr <- efficient.rounding(rep(0.25, 4), 5)
hbd_5_mu_pso <- lapply(1:length(hbd_5_effr), function(x) rep(hbd_support[x], hbd_5_effr[x])) |> unlist()
mvnorm_sd = 0.15/4
nswarm <- psoinfo_exact$nSwarm/2
hbd_5_mvnorm_pso <- mvrnorm(n = (nswarm-1), mu = hbd_5_mu_pso, Sigma = diag(5) * mvnorm_sd)
hbd_5_mvnorm_pso[hbd_5_mvnorm_pso < 0] <- 0
hbd_5_mvnorm_pso[hbd_5_mvnorm_pso > 0.15] <- 0.15
hbd_5_mvnorm_pso <- rbind(hbd_5_mvnorm_pso, hbd_5_mu_pso)
hbd_5_result <- globpso(objFunc = hb_doptimal, lower = rep(0, 5), upper = rep(0.15, 5), PSO_INFO = psoinfo_exact, 
                             loc = hb_par, verbose = T, init = hbd_5_mvnorm_pso)
plot(x = c(1:501), y = hbd_5_result$history / hbd_base, type = "l", xlab = "Iteration", ylab = "D-efficiency", ylim = c(0, 1))

hbd_10_effr <- efficient.rounding(rep(0.25, 4), 10)
hbd_10_mu_pso <- lapply(1:length(hbd_10_effr), function(x) rep(hbd_support[x], hbd_10_effr[x])) |> unlist()
hbd_10_mvnorm_pso <- mvrnorm(n = 15, mu = hbd_10_mu_pso, Sigma = diag(10) * mvnorm_sd)
hbd_10_mvnorm_pso[hbd_10_mvnorm_pso < 0] <- 0
hbd_10_mvnorm_pso[hbd_10_mvnorm_pso > 0.15] <- 0.15
hbd_10_mvnorm_pso <- rbind(hbd_10_mvnorm_pso, hbd_10_mu_pso)
hbd_10_result <- globpso(objFunc = hb_doptimal, lower = rep(0, 10), upper = rep(0.15, 10), PSO_INFO = psoinfo_exact, 
                        loc = hb_par, verbose = T, init = hbd_10_mvnorm_pso)
plot(x = c(1:501), y = hbd_10_result$history / hbd_base, type = "l", xlab = "Iteration", ylab = "D-efficiency", ylim = c(0, 1))

hbd_15_effr <- efficient.rounding(rep(0.25, 4), 15)
hbd_15_mu_pso <- lapply(1:length(hbd_15_effr), function(x) rep(hbd_support[x], hbd_15_effr[x])) |> unlist()
hbd_15_mvnorm_pso <- mvrnorm(n = 15, mu = hbd_15_mu_pso, Sigma = diag(15) * mvnorm_sd)
hbd_15_mvnorm_pso[hbd_15_mvnorm_pso < 0] <- 0
hbd_15_mvnorm_pso[hbd_15_mvnorm_pso > 0.15] <- 0.15
hbd_15_mvnorm_pso <- rbind(hbd_15_mvnorm_pso, hbd_15_mu_pso)
hbd_15_result <- globpso(objFunc = hb_doptimal, lower = rep(0, 15), upper = rep(0.15, 15), PSO_INFO = psoinfo_exact, 
                         loc = hb_par, verbose = T, init = hbd_15_mvnorm_pso)
plot(x = c(1:501), y = hbd_15_result$history / hbd_base, type = "l", xlab = "Iteration", ylab = "D-efficiency", ylim = c(0, 1))

hbd_20_effr <- efficient.rounding(rep(0.25, 4), 20)
hbd_20_mu_pso <- lapply(1:length(hbd_20_effr), function(x) rep(hbd_support[x], hbd_20_effr[x])) |> unlist()
hbd_20_mvnorm_pso <- mvrnorm(n = 31, mu = hbd_20_mu_pso, Sigma = diag(20) * mvnorm_sd)
hbd_20_mvnorm_pso[hbd_20_mvnorm_pso < 0] <- 0
hbd_20_mvnorm_pso[hbd_20_mvnorm_pso > 0.15] <- 0.15
hbd_20_mvnorm_pso <- rbind(hbd_20_mvnorm_pso, hbd_20_mu_pso)
hbd_20_result <- globpso(objFunc = hb_doptimal, lower = rep(0, 20), upper = rep(0.15, 20), PSO_INFO = psoinfo_exact2, 
                         loc = hb_par, verbose = T, init = hbd_20_mvnorm_pso)
plot(x = c(1:1001), y = hbd_20_result$history / hbd_base, type = "l", xlab = "Iteration", ylab = "D-efficiency", ylim = c(0, 1))

hbh_support <- c(0, 0.02, 0.04)
hbh_weight <- c(0.359, 0.5, 0.141)
hbh_5_effr <- efficient.rounding(hbh_weight, 5)
hbh_5_mu_pso <- lapply(1:length(hbh_5_effr), function(x) rep(hbh_support[x], hbh_5_effr[x])) |> unlist()
hbh_5_mvnorm_pso <- mvrnorm(n = 15, mu = hbh_5_mu_pso, Sigma = diag(5) * mvnorm_sd)
hbh_5_mvnorm_pso[hbh_5_mvnorm_pso < 0] <- 0
hbh_5_mvnorm_pso[hbh_5_mvnorm_pso > 0.15] <- 0.15
hbh_5_mvnorm_pso <- rbind(hbh_5_mvnorm_pso, hbh_5_mu_pso)
hbh_5_result <- globpso(objFunc = hb_hoptimal, lower = rep(0, 5), upper = rep(0.15, 5), PSO_INFO = psoinfo_exact, 
                        loc = hb_par, verbose = T, init = hbh_5_mvnorm_pso)
hbh_5_result$par|>round(4)|>table()
plot(x = c(1:501), y = hbh_base/hbh_5_result$history, type = "l", xlab = "Iteration", ylab = "h-efficiency", ylim = c(0.9, 1))

hbh_10_effr <- efficient.rounding(hbh_weight, 10)
hbh_10_mu_pso <- lapply(1:length(hbh_10_effr), function(x) rep(hbh_support[x], hbh_10_effr[x])) |> unlist()
hbh_10_mvnorm_pso <- mvrnorm(n = 15, mu = hbh_10_mu_pso, Sigma = diag(10) * mvnorm_sd)
hbh_10_mvnorm_pso[hbh_10_mvnorm_pso < 0] <- 0
hbh_10_mvnorm_pso[hbh_10_mvnorm_pso > 0.15] <- 0.15
hbh_10_mvnorm_pso <- rbind(hbh_10_mvnorm_pso, hbh_10_mu_pso)
hbh_10_result <- globpso(objFunc = hb_hoptimal, lower = rep(0, 10), upper = rep(0.15, 10), PSO_INFO = psoinfo_exact, 
                        loc = hb_par, verbose = T, init = hbh_10_mvnorm_pso)
hbh_10_result$par|>round(4)|>table()
plot(x = c(1:501), y = hbh_base/hbh_10_result$history, type = "l", xlab = "Iteration", ylab = "h-efficiency", ylim = c(0.9, 1))

hbh_15_effr <- efficient.rounding(hbh_weight, 15)
hbh_15_mu_pso <- lapply(1:length(hbh_15_effr), function(x) rep(hbh_support[x], hbh_15_effr[x])) |> unlist()
hbh_15_mvnorm_pso <- mvrnorm(n = 15, mu = hbh_15_mu_pso, Sigma = diag(15) * mvnorm_sd)
hbh_15_mvnorm_pso[hbh_15_mvnorm_pso < 0] <- 0
hbh_15_mvnorm_pso[hbh_15_mvnorm_pso > 0.15] <- 0.15
hbh_15_mvnorm_pso <- rbind(hbh_15_mvnorm_pso, hbh_15_mu_pso)
hbh_15_result <- globpso(objFunc = hb_hoptimal, lower = rep(0, 15), upper = rep(0.15, 15), PSO_INFO = psoinfo_exact, 
                         loc = hb_par, verbose = T, init = hbh_15_mvnorm_pso)
hbh_15_result$par|>round(4)|>table()
plot(x = c(1:501), y = hbh_base/hbh_15_result$history, type = "l", xlab = "Iteration", ylab = "h-efficiency", ylim = c(0.9, 1))

hbh_20_effr <- efficient.rounding(hbh_weight, 20)
hbh_20_mu_pso <- lapply(1:length(hbh_20_effr), function(x) rep(hbh_support[x], hbh_20_effr[x])) |> unlist()
hbh_20_mvnorm_pso <- mvrnorm(n = 31, mu = hbh_20_mu_pso, Sigma = diag(20) * mvnorm_sd)
hbh_20_mvnorm_pso[hbh_20_mvnorm_pso < 0] <- 0
hbh_20_mvnorm_pso[hbh_20_mvnorm_pso > 0.15] <- 0.15
hbh_20_mvnorm_pso <- rbind(hbh_20_mvnorm_pso, hbh_20_mu_pso)
hbh_20_result <- globpso(objFunc = hb_hoptimal, lower = rep(0, 20), upper = rep(0.15, 20), PSO_INFO = psoinfo_exact2, 
                         loc = hb_par, verbose = T, init = hbh_20_mvnorm_pso)
hbh_20_result$par|>round(4)|>table()
plot(x = c(1:1001), y = hbh_base/hbh_20_result$history, type = "l", xlab = "Iteration", ylab = "h-efficiency", ylim = c(0.9, 1))

logistic_plot(log_par, lb = -5, ub = 5)
logistic_plot(c(2, 5), lb = -1, ub = 1)
logistic_plot(c(-3, 10), lb = 0, ub = 1)

logistic_plot(c(2, 2), lb = -5, ub = 5)
logistic_plot(c(-3, 20), lb = 0, ub = 1)

globpso(objFunc = logit_doptimal, lower = rep(-1, 2), upper = rep(1, 2), PSO_INFO = psoinfo_exact, 
        loc = c(2, 5), verbose = T)$par |> round(4)

globpso(objFunc = logit_doptimal, lower = rep(0, 2), upper = rep(1, 2), PSO_INFO = psoinfo_exact, 
        loc = c(-3, 10), verbose = T)$par |> round(4)


hb_doptimal_round <- function(input, loc){
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  nobs <- loc[5]
  
  n <- length(input)/2
  d <- input[1:n]
  w <- input[(n+1):(2*n)] |> round()
  #print(w)
  
  if(sum(w) !=  nobs) res = 1e+10
  else{
    mat_list <- lapply(1:n, function(x) w[x]/n * hb_mat(d[x], c1, tau, b0, b1))
    inf_mat <- Reduce("+", mat_list)
    res = -det(inf_mat)
  }
  res
}

hb_d_round <- globpso(hb_doptimal_round, lower = c(rep(0,8)), upper = c(rep(0.15, 4), rep(10, 4)), PSO_INFO = psoinfo, 
                      loc = c(hb_par, 10))
c(round(hb_d_round$par[1,1:4], 4), round(hb_d_round$par[1,5:8]))


hb_doptimal_round <- function(input, loc){
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  nobs <- loc[5]
  
  n <- length(input)/2
  d <- input[1:n]
  w <- input[(n+1):(2*n)] |> round()
  #print(w)
  
  if(sum(w) !=  nobs) res = 1e+10
  else{
    mat_list <- lapply(1:n, function(x) w[x]/n * hb_mat(d[x], c1, tau, b0, b1))
    inf_mat <- Reduce("+", mat_list)
    res = -det(inf_mat)
  }
  res
}

hb_d_round <- globpso(hb_doptimal_round, lower = c(rep(0,8)), upper = c(rep(0.15, 4), rep(20, 4)), PSO_INFO = psoinfo, 
                      loc = c(hb_par, 20))
c(round(hb_d_round$par[1,1:4], 4), round(hb_d_round$par[1,5:8]))

hb_hoptimal_round <- function(input, loc){
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  nobs <- loc[5]
  
  n <- length(input)/2
  d <- input[1:n]
  w <- input[(n+1):(2*n)] |> round()
  #print(w)
  
  mat_list <- lapply(1:n, function(x) w[x]/n * hb_mat(d[x], c1, tau, b0, b1))
  M <- Reduce("+", mat_list)
  diag(M) <- diag(M) + 1e-10
  h <- matrix(c(-tau, -c1, 0, 0))
  
  if (rcond(M) < 2.220446e-16 || sum(w) !=  nobs){
    res = 9e+10
  } 
  else{
    M_inv <- solve(M)
    res <- t(h) %*% M_inv %*% h
  }
  res
}

hb_h_round <- globpso(hb_hoptimal_round, lower = c(rep(0,8)), upper = c(rep(0.15, 4), rep(8, 4)), PSO_INFO = psoinfo, 
                      loc = c(hb_par, 8))
c(round(hb_h_round$par[1,1:4], 4), round(hb_h_round$par[1,5:8]))
