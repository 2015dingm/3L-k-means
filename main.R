rm(list=ls())
source("algorithm.R")
source("function.R")
# library('LICORS')
library('wskm')

x <- read.csv(file="../Data/new/User_del0.csv", header = T)
e <- read.csv(file="../Data/new/Energy.csv", header = T)
x = as.matrix(x)
e <- as.matrix(e)

print('Read finish!')

nr = nrow(x)
nc = ncol(x)
ng = 33
k = 10
lambda = 50000
eta = 10000

x_scale <- scale(x)
scale_center = attr(x_scale, "scaled:center")
scale_std = attr(x_scale, "scaled:scale")

e_scale <- e
e_mean = mean(e)
e_sd = sd(e)
e_min = min(e)
e_max = e_mean + 3 * e_sd
for (i in seq(1, nr)){
  e_scale[i] = (e[i]-e_min)/(e_max-e_min)
}

dum = matrix(nrow = nr, ncol = 1)
for (i in seq(1, nr)){
  dum[i] = 1
}

group_p <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,4,5,5,5,6,6,6,6,6,6,6,6,6,6,6,
             7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,
             11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,
             15,15,15,15,15,15,15,15,15,15,15,15,16,16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,17,17,17,18,18,18,18,18,18,18,18,18,18,18,
             19,19,19,19,19,19,19,19,19,19,19,19,20,20,20,20,20,20,20,20,20,20,20,20,21,21,21,21,21,21,21,21,21,21,21,21,22,22,22,22,22,22,22,22,22,22,22,22,
             23,23,23,23,23,23,23,23,23,23,23,23,24,24,24,24,24,24,24,24,24,24,24,24,25,25,25,25,25,25,25,25,25,25,25,26,26,26,26,26,26,26,26,26,26,26,26,
             27,27,27,27,27,27,27,27,28,28,28,28,28,28,28,28,28,29,29,29,29,29,29,29,29,29,29,29,29,30,30,30,30,30,30,30,30,30,30,30,30,
             31,31,31,31,31,31,31,31,31,31,31,31,32,32,32,32,32,32,32,32,32,32,32,32)

group_hw <- c(0,1,2,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,0,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,9,7,8,9,0,1,3,4,5,6,7,8,9,10,11,
              0,1,2,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,0,1,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,
              0,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,0,1,3,4,5,6,7,8,9,10,11,
              0,1,2,3,4,5,6,7,8,9,10,11,0,1,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,
              0,1,2,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,0,1,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,
              0,1,2,3,5,6,8,11,0,3,4,5,6,7,8,9,11,0,1,2,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,
              0,1,2,3,4,5,6,7,8,9,10,11)

x_p = aggregate_function(x, group_p)
x_hw = aggregate_function(x, group_hw)

# init_centers <- kmeanspp(x, k, start = "normal", nstart = 1)
# write.csv(init_centers[["inicial.centers"]], file="../Result/new/inicial.csv")

# mymfgkm_ini <- mfgkm(x_scale, e_scale, init_centers$inicial.centers, group_p, lambda = lambda, eta = eta, maxiter = 500, delta = 1e-15, seed = 22)
# mymfgkm <- mfgkm(x_scale, e_scale, k, group_p, lambda = lambda, eta = eta, maxiter = 500, delta = 1e-15, seed = 22)
# mymtwkm <- mtwkm(x_scale, e_scale, k, group_p, lambda = lambda, eta = eta, maxiter = 500, delta = 1e-15, seed = 22)
# myfgkm <- fgkm(x_scale, k, group_p, lambda = lambda, eta = eta, maxiter = 500, delta = 1e-15, seed = 22)
# mytwkm <- twkm(x_scale, k, group_p, lambda = lambda, eta = eta, maxiter = 500, delta = 1e-15, seed = 22)
# mykm <- kmeans(x_scale, k)
# 
# output_function(mymfgkm_ini, "mfgkm_ini(p)")
# output_function(mymfgkm, "mfgkm(p)")
# output_function(mymtwkm, "mtwkm(p)")
# output_function(myfgkm, "fgkm(p)")
# output_function(mytwkm, "twkm(p)")
# output_function(mykm, "km(p)")

# evaluation_function(mymtwkm, x, x_p, x_hw, e, 50, nr)

# cost_mtwkm = cost_function_twkm(mymtwkm, x_scale, e_scale, group_p, nr, nc, ng, k, lambda, eta)
# cost_twkm = cost_function_twkm(mytwkm, x_scale, e_scale, group_p, nr, nc, ng, k, lambda, eta)

# lambda_range = c(100, 500, 1000, 5000, 10000, 50000, 100000)
# eta_range = c(10, 50, 100, 500, 1000, 5000, 10000)

lambda_range = c(10000, 50000, 100000)
eta_range = c(50000, 100000)
k_range = c(15, 25, 35, 45)
tuning_function(x, x_scale, x_p, x_hw, e, e_scale, group_p, nr, nc, ng, k_range, lambda_range, eta_range, 50, "mtwkm(p)")
  
