rm(list=ls())
dyn.load("mfgkm.so")
dyn.load("mtwkm.so")
library('LICORS')
library('wskm')

mfgkm <- function(x, e, centers, groups, lambda, eta, maxiter=100, delta=0.000001, maxrestart=10,seed=-1) 
{
  if (missing(centers))
    stop("the number or initial clusters 'centers' must be provided")

  if(seed<=0){
    seed <-runif(1,0,10000000)[1]
  }

  vars <- colnames(x)
  
  nr <-nrow(x) # nrow() return a integer type
  nc <-ncol(x) # integer

  if (is.data.frame(centers) || is.matrix(centers))
  {
    init <- TRUE
    k <- nrow(centers)
  }
  else
  {
    init <- FALSE
    k <- centers
    centers <- double(k * nc)
  }
  
  # get the setting of feature group
  if (is.character(groups) && length(groups) == 1) {
    G <- .C("parseGroup",as.character(groups),numGroups=integer(1), groupInfo=integer(nc),PACKAGE="wskm")
  } else if (is.vector(groups) && length(groups) == nc) {
    G <- list()
    grps <- as.factor(groups)
    groupNames <- levels(grps)
    G$numGroups <- nlevels(grps)
    G$groupInfo <- as.integer(as.integer(grps) - 1)
  }

  set.seed(seed)
  Z <- .C("mfgkm",
          x = as.double(as.matrix(x)),
          e = as.double(as.matrix(e)),
          nr,
          nc,
          k = as.integer(k),
          lambda = as.double(lambda),
          eta = as.double(eta),
          G$numGroups,
          G$groupInfo,
          delta = as.double(delta),
          maxIterations = as.integer(maxiter),
          maxRestarts = as.integer(maxrestart),
          as.logical(init),
#          seed,
          cluster = integer(nr),
          centers=as.double(as.matrix(centers)),
          featureWeight = double(k * nc),
          groupWeight = double(k * G$numGroups),
          iterations = integer(1),
          restarts = integer(1),
          totiters = integer(1),
          totalCost = double(1),
          totss = double(1),
		      withiness = double(k)
          # PACKAGE="wskm"
          )
       
  centers <- matrix( Z$centers)
  dim(centers) <- c(k, nc)
  colnames(centers) <- vars
  
  featureWeight <- matrix(Z$featureWeight)
  dim(featureWeight) <- c(k, nc)
  colnames(featureWeight) <- vars
  
  groupWeight <- matrix(Z$groupWeight)
  dim(groupWeight) <- c(k, G$numGroups )
  colnames(groupWeight) <- 1:ncol(groupWeight)
  
  ignore <- which(rowSums(centers==0) == ncol(centers))
  if (length(ignore)) {
    centers <- centers[-ignore,, drop=FALSE]
    featureWeight <- featureWeight[-ignore,, drop=FALSE]
  }
  
  rownames(centers) <- 1:nrow(centers)
  rownames(featureWeight) <- 1:nrow(featureWeight)
  rownames(groupWeight) <- 1:nrow(groupWeight)
  
  cluster <- Z$cluster + 1
  
  size <- aggregate(cluster, list(cluster=cluster), length)[[2]]
  
  result <- list(cluster = cluster,
                 centers = Z$centers,
                 totss = Z$totss, 
                 withinss = Z$withinss, 
                 tot.withinss = sum(Z$withiness), 
                 betweenss = Z$totss-sum(Z$withinss),
                 size = size,
                 iterations = Z$iterations,
                 restarts = Z$restarts,
                 totiters=Z$totiters,
                 featureWeight = Z$featureWeight,
                 groupWeight = Z$groupWeight)
  
  dim(result$centers) <- c(k, nc)
  dim(result$featureWeight) <- c(k, nc)
  dim(result$groupWeight) <- c(k, G$numGroups)
  
  class(result) <- c("kmeans", "mfgkm")
  return(result)
}

mtwkm <- function(x, e, centers, groups, lambda, eta, maxiter=100, delta=0.000001, maxrestart=10,seed=-1) 
{
  if (missing(centers))
    stop("the number or initial clusters 'centers' must be provided")
  
  if(seed<=0){
    seed <-runif(1,0,10000000)[1]
  }
  
  vars <- colnames(x)
  
  nr <-nrow(x) # nrow() return a integer type
  nc <-ncol(x) # integer
  
  if (is.data.frame(centers) || is.matrix(centers))
  {
    init <- TRUE
    k <- nrow(centers)
  }
  else
  {
    init <- FALSE
    k <- centers
    centers <- double(k * nc)
  }
  
  # get the setting of feature group
  if (is.character(groups) && length(groups) == 1) {
    G <- .C("parseGroup",as.character(groups),numGroups=integer(1), groupInfo=integer(nc),PACKAGE="wskm")
  } else if (is.vector(groups) && length(groups) == nc) {
    G <- list()
    grps <- as.factor(groups)
    groupNames <- levels(grps)
    G$numGroups <- nlevels(grps)
    G$groupInfo <- as.integer(as.integer(grps) - 1)
  }
  
  set.seed(seed)
  Z <- .C("mtwkm",
          x = as.double(as.matrix(x)),
          e = as.double(as.matrix(e)),
          nr,
          nc,
          k = as.integer(k),
          lambda = as.double(lambda),
          eta = as.double(eta),
          G$numGroups,
          G$groupInfo,
          delta = as.double(delta),
          maxIterations = as.integer(maxiter),
          maxRestarts = as.integer(maxrestart),
          as.logical(init),
          #          seed,
          cluster = integer(nr),
          centers=as.double(as.matrix(centers)),
          featureWeight = double( nc),
          groupWeight = double( G$numGroups),
          iterations = integer(1),
          restarts = integer(1),
          totiters = integer(1),
          totalCost = double(1),
          totss = double(1),
          withiness = double(k)
#          PACKAGE="wskm"
  )
  
  centers <- matrix( Z$centers)
  dim(centers) <- c(k, nc)
  colnames(centers) <- vars
  
  featureWeight <- Z$featureWeight
  
  groupWeight <- Z$groupWeight
  
  ignore <- which(rowSums(centers==0) == ncol(centers))
  if (length(ignore)) {
    centers <- centers[-ignore,, drop=FALSE]
    featureWeight <- featureWeight[-ignore,, drop=FALSE]
  }
  
  rownames(centers) <- 1:nrow(centers)
  
  cluster <- Z$cluster + 1
  
  size <- aggregate(cluster, list(cluster=cluster), length)[[2]]
  
  result <- list(cluster = cluster,
                 centers = Z$centers,
                 totss = Z$totss, 
                 withinss = Z$withinss, 
                 tot.withinss = sum(Z$withiness), 
                 betweenss = Z$totss-sum(Z$withinss),
                 size = size,
                 iterations = Z$iterations,
                 restarts = Z$restarts,
                 totiters=Z$totiters,
                 featureWeight = Z$featureWeight,
                 groupWeight = Z$groupWeight)
  
  dim(result$centers) <- c(k, nc)
  
  class(result) <- c("kmeans", "mtwkm")
  return(result)
}

#########################################################################

# x <- read.csv(file="../Data/new/User_del0_del_cut.csv", header = T)
# e <- read.csv(file="../Data/new/Energy_cut.csv", header = T)
x <- read.csv(file="../Data/new/User_del0.csv", header = T)
e <- read.csv(file="../Data/new/Energy.csv", header = T)

nr = nrow(x)
nc = ncol(x)
ng = 33
k = 10
# lambda = 12
# eta = 2
lambda = 50000
eta = 10000
x = as.matrix(x)
e <- as.matrix(e)

x <- scale(x)
scale_center = attr(x, "scaled:center")
scale_std = attr(x, "scaled:scale")

e <- as.matrix(e)
e_mean = mean(e)
e_sd = sd(e)
e_min = min(e)
e_max = e_mean + 3 * e_sd
for (i in seq(1, nr)){
  e[i] = (e[i]-e_min)/(e_max-e_min)
}


dum = matrix(nrow = nr, ncol = 1)
for (i in seq(1, nr)){
  dum[i] = 1
}


inv_scale_function <- function(x, scale_center, scale_std){
  inv_scale_x <- matrix(nrow = dim(x)[1], ncol = dim(x)[2])
  for (c in seq(1, dim(x)[2])){
    for (r in seq(1, dim(x)[1])){
      inv_scale_x[r, c] = x[r, c] * scale_std[c] + scale_center[c]
    }
  }
  return(inv_scale_x)
}


aggregate_function <- function(x, group, ng){
  inv_scale_x = t(inv_scale_function(x, scale_center, scale_std))
  group = as.matrix(group)
  temp = cbind(inv_scale_x, group)
  x_g = t(aggregate(temp[,seq(1, dim(temp)[2]-1)],list(temp[,dim(temp)[2]]),sum))
  
  return(x_g)
}


output_function <- function(result, file_name){
  write.csv(result$featureWeight, file=paste("../Result/new/", file_name, "_fw.csv", sep = ''))
  write.csv(result$groupWeight, file=paste("../Result/new/", file_name, "_gw.csv", sep = ''))
  write.csv(result$cluster, file=paste("../Result/new/", file_name, "_cluster.csv", sep = ''))
  # write.csv(result$centers, file=paste("../Result/new/", file_name, "_center.csv", sep = ''))
  centers_scaled = result$centers
  centers_noscaled <- matrix(nrow = k, ncol = nc)
  for (c in seq(1, nc)){
    for (r in seq(1, k)){
      centers_noscaled[r, c] = centers_scaled[r, c] * scale_std[c] + scale_center[c]
    }
  }
  write.csv(centers_noscaled, file=paste("../Result/new/", file_name, "_center.csv", sep = ''))
}

cost_function <- function(result, x, e=dum, group, nr, nc, ng, k, lambda, eta){
  sum1 = 0
  sum2 = 0
  sum3 = 0
  for (i in seq(1,nr)){
    for (j in seq(1,nc)){
      sum1 = sum1 + e[i] * result$groupWeight[group[j] * k + result$cluster[i]] * result$featureWeight[(j-1) * k + result$cluster[i]] * (x[(j-1) * nr + i] - result$centers[(j-1) * k + result$cluster[i]])^2
    }
  }
  
  for (l in seq(1,k)){
    for (t in seq(1,ng)){
      sum2 = sum2 + result$groupWeight[(t-1) * k + l] * log(result$groupWeight[(t-1) * k + l])
    }
    
    for (j in seq(1,nc)){
      sum3 = sum3 + result$featureWeight[(j-1) * k + l] * log(result$featureWeight[(j-1) * k + l])
    }
  }
  sum2 = sum2 * lambda
  sum3 = sum3 * eta
  
  cost = c(sum1, sum2, sum3)
  
  return(cost)
}

distortion_function <- function(result, center, x, e, nr, scale_center, scale_std){
  user = inv_scale_function(x, scale_center, scale_std)
  center = inv_scale_function(center, scale_center, scale_std)
  sum = 0
  tot = 0
  for (i in seq(1, nr)){
    distortion = dist(rbind(user[i+1,], center[result$cluster[i]+1,]), p=2)
    distortion = distortion * e[i,1]
    sum = sum + distortion
    tot = tot + e[i,1]
  }
  md = sum/tot
  
  return(md)
}

evaluation_function <- function(result, x, x_p, x_hw, e, ceiling, nr, scale_center, scale_std){
  k_real = length(result$size)
  
  # entropy
  max = 0
  for (i in seq(1, k_real)){
    if (result$size[i] > max) {
      max <- result$size[i]
    }
  }
  entropy = (nr-max)/nr
  
  # complexity
  complexity = (ceiling - k_real)/ceiling
  
  # distortion
  center = result$centers
  center_p = aggregate_function(result$centers, group_p, 33)
  center_hw = aggregate_function(result$centers, group_hw, 12)
  
  mdph = distortion_function(result, center, x, e, nr, scale_center, scale_std)
  mdh = distortion_function(result, center_hw, x_hw, e, nr, scale_center, scale_std)
  mdp = distortion_function(result, center_p, x_p, e, nr, scale_center, scale_std)
  
  evaluation = c(entropy, complexity, mdph, mdh, mdp)
  
  return(evaluation)
}

tuning_function <- function(result, x, x_p, x_hw, e, group, nr, nc, ng, k_range, lambda_range, eta_range, ceiling, scale_center, scale_std, file_name){
  df <- data.frame(lambda = numeric(0), eta = numeric(0), k = numeric(0), k_real = numeric(0), size = character(0), 
                   var_fw = numeric(0), var_gw = numeric(0), sum1 = numeric(0), sum2 = numeric(0), sum3 = numeric(0), 
                   entropy = numeric(0), complexity = numeric(0), mdph = numeric(0), mdh = numeric(0), mdp = numeric(0), 
                   stringsAsFactors=FALSE)
  index = 1
  # for (lam in c(1e-3, 1e-2, 1e-1, 1, 2, 4, 8, 12, 16, 24, 32, 48, 64, 80, 100)){
  #   for (et in c(1e-3, 1e-2, 1e-1, 1, 2, 4, 8, 12, 16, 24, 32, 48, 64, 80, 100, 200)){
  #     for (k in seq(5,20)){
  for (lam in lambda_range){
    for (et in eta_range){
      for (k in k_range){
        # init_centers <- kmeanspp(x, k, start = "normal")
        # mymfgkm <- mfgkm(x, e, init_centers$inicial.centers, group_p, lambda = lam, eta = et, maxiter = 500, delta = 1e-10, seed = 22)
        result <- mtwkm(x, e, k, group, lambda = lam, eta = et, maxiter = 500, delta = 1e-15, seed = 22)
        fw = as.vector(result$featureWeight)
        gw = as.vector(result$groupWeight)
        k_real = length(result$size)
        size = ''
        for (i in seq(k_real)){
          size = paste(size, as.character(result$size[i]), sep = ',')
        }
        cost = cost_function(result, x, e, group, nr, nc, ng, k, lam, et)
        evaluation = evaluation_function(result, x, x_p, x_hw, e, ceiling, nr, scale_center, scale_std)
        
        df[index,] = c(lam, et, k, k_real, size, var(fw), var(gw), 
                       cost[1], cost[2], cost[3], evaluation[1], evaluation[2], evaluation[3], evaluation[4], evaluation[5])
        print(c('lambda=', lam, 'eta=', et, 'k=', k, 'k_real=', k_real, 'size=', size, 
                'var(fg)=', var(fw), 'var(gw)=', var(gw), 'sum1=', cost[1], 'sum2=', cost[2], 'sum3=', cost[3], 
                'entropy=', evaluation[1], 'complexity=', evaluation[2], 'mdph=', evaluation[3], 'mdh=', evaluation[4], 'mdp=', evaluation[5]))
        
        index = index + 1
      }
    }
  }
  write.csv(df, file=paste("../Result/new/tuning1_", file_name, ".csv", sep = ''))
}



# group <- c(0,1,3,4,5,6,7,8,9,11,0,1,2,3,4,5,6,8,10,11,0,1,2,3,4,5,6,8,10,11,0,4,5,6,7,8,9,11,5,8,0,4,5,6,7,8,9,11,4,5,6,7,8,9,11,0,1,4,5,6,7,8,9,10,11,0,5,6,7,8,9,11,0,4,5,6,7,8,9,11,0,1,5,6,7,8,9,10,11,9,0,7,8,9,0,1,3,4,5,6,7,8,9,10,11,0,1,3,4,5,6,7,8,9,10,11,0,1,3,4,5,6,7,8,9,10,11,0,4,5,6,7,8,9,11,0,3,4,5,6,7,8,9,10,11,0,4,5,6,7,8,9,10,11,0,4,5,6,7,8,9,10,11,0,3,6,7,8,9,11,0,1,3,4,5,6,7,8,9,10,11,0,3,4,5,6,7,8,9,11,0,1,3,4,5,6,7,8,9,10,11,0,3,4,5,6,7,8,9,11,0,3,5,6,7,8,9,11,0,1,2,3,4,5,6,7,8,9,10,11,0,3,4,5,6,7,8,9,10,11,0,1,3,4,5,6,7,8,9,10,11,0,1,3,4,5,6,7,8,9,10,11,0,3,5,6,7,8,11,0,1,4,5,6,7,8,9,11,0,1,2,4,5,6,7,8,9,10,11,0,1,3,5,6,7,8,9,10,11,0,1,3,4,5,6,7,8,9,10,11,0,1,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,0,4,5,6,7,8,9,10,11,0,5,6,7,8,9,10,11,0,4,5,6,7,8,9,10,11,0,6,7,8,9,11,0,1,2,3,4,5,6,7,8,9,10,11,0,5,6,7,8,9,11,0,6,7,8,9,11,0,5,6,8,9,11,0,1,4,5,6,7,8,9,11,0,3,6,7,8,9,11,0,1,2,3,5,6,7,8,11,0,1,2,3,6,8,11,0,3,4,5,6,7,8,9,11,0,1,2,3,4,5,6,7,8,9,10,11)
#### all
# group <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,11,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,15,16,16,16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,17,17,17,18,18,18,18,18,18,18,18,18,18,18,18,19,19,19,19,19,19,19,19,19,19,19,20,20,20,20,20,20,20,20,20,21,21,21,21,21,21,21,21,21,21,21,21,22,22,22,22,22,22,22,22,22,22,22,22,23,23,23,23,23,23,23,23,23,23,23,23,24,24,24,24,24,24,24,24,24,24,24,25,25,25,25,25,25,25,25,25,26,26,26,26,26,26,26,26,26,26,26,26,27,27,27,27,27,27,27,27,27,27,27,27,28,28,28,28,28,28,28,28,28,28,28,28,29,29,29,29,29,29,29,29,29,29,29,29,30,30,30,30,30,30,30,30,30,30,30,31,31,31,31,31,31,31,31,31,31,31,31,32,32,32,32,32,32,32,32,32,32,32,32,33,33,33,33,33,33,33,33,33,33,33,34,34,34,34,34,34,34,34,34,34,34,34,35,35,35,35,35,35,35,35,35,35,35,35,36,36,36,36,36,36,36,36,36,36,36,36,37,37,37,37,37,37,37,37,37,37,37,38,38,38,38,38,38,38,38,38,38,38,38,39,39,39,39,39,39,39,39,39,39,39,39,40,40,40,40,40,40,40,40,40,40,40,41,41,41,41,41,41,41,41,41,41,41,41,42,42,42,42,42,42,42,42,42,42,42,42,43,43,43,43,43,43,43,43,43,43,43,43,44,44,44,44,44,44,44,44,44,44,44,45,45,45,45,45,45,45,45,45,45,45,45,46,46,46,46,46,46,46,46,46,47,47,47,47,47,47,47,47,47,47,47,47,48,48,48,48,48,48,48,48,49,49,49,49,49,49,49,49,49,49,49,49,50,50,50,50,50,50,50,50,50,50,50,50)
#### cut
# group <- c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,11,12,12,12,12,13,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,17,18,18,18,18,18,18,18,18,18,19,19,19,19,19,19,19,19,19,20,20,20,20,20,20,20,21,21,21,21,21,21,21,21,21,21,21,22,22,22,22,22,22,22,22,22,23,23,23,23,23,23,23,23,23,23,23,24,24,24,24,24,24,24,24,24,25,25,25,25,25,25,25,25,26,26,26,26,26,26,26,26,26,26,26,26,27,27,27,27,27,27,27,27,27,27,28,28,28,28,28,28,28,28,28,28,28,29,29,29,29,29,29,29,29,29,29,29,30,30,30,30,30,30,30,31,31,31,31,31,31,31,31,31,32,32,32,32,32,32,32,32,32,32,32,33,33,33,33,33,33,33,33,33,33,34,34,34,34,34,34,34,34,34,34,34,35,35,35,35,35,35,35,35,35,35,35,36,36,36,36,36,36,36,36,36,36,36,36,37,37,37,37,37,37,37,37,37,38,38,38,38,38,38,38,38,39,39,39,39,39,39,39,39,39,40,40,40,40,40,40,41,41,41,41,41,41,41,41,41,41,41,41,42,42,42,42,42,42,42,43,43,43,43,43,43,44,44,44,44,44,44,45,45,45,45,45,45,45,45,45,46,46,46,46,46,46,46,47,47,47,47,47,47,47,47,47,48,48,48,48,48,48,48,49,49,49,49,49,49,49,49,49,50,50,50,50,50,50,50,50,50,50,50,50)

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

x_p = aggregate_function(x, group_p, 33)
x_hw = aggregate_function(x, group_hw, 12)

# init_centers <- kmeanspp(x, k, start = "normal", nstart = 1)
# write.csv(init_centers[["inicial.centers"]], file="../Result/new/inicial.csv")

# mymfgkm_ini <- mfgkm(x, e, init_centers$inicial.centers, group_p, lambda = lambda, eta = eta, maxiter = 500, delta = 1e-15, seed = 22)
# mymfgkm <- mfgkm(x, e, k, group_p, lambda = lambda, eta = eta, maxiter = 500, delta = 1e-15, seed = 22)
mymtwkm <- mtwkm(x, e, k, group_p, lambda = lambda, eta = eta, maxiter = 500, delta = 1e-15, seed = 22)
# myfgkm <- fgkm(x, k, group_p, lambda = lambda, eta = eta, maxiter = 500, delta = 1e-15, seed = 22)
# mytwkm <- twkm(x, k, group_p, lambda = lambda, eta = eta, maxiter = 500, delta = 1e-15, seed = 22)
# mykm <- kmeans(x, k)
# 
# output_function(mymfgkm_ini, "mfgkm_ini(p)")
# output_function(mymfgkm, "mfgkm(p)")
# output_function(mymtwkm, "mtwkm(p)")
# output_function(myfgkm, "fgkm(p)")
# output_function(mytwkm, "twkm(p)")
# output_function(mykm, "km(p)")

# evaluation_function(mymtwkm, x, x_p, x_hw, e, 50, nr, scale_center, scale_std)

# cost_mfgkm = cost_function(mymfgkm, x, e, group_p, nr, nc, ng, k, lambda, eta)
# cost_fgkm = cost_function(myfgkm, x, dum, group_p, nr, nc, ng, k, lambda, eta)

# lambda_range = c(1e-2, 1e-1, 1, 2, 4, 8, 12, 16, 24, 32, 48, 64, 80, 100, 500, 1000, 10000)
# eta_range = c(1e-3, 1e-2, 1e-1, 1, 2, 4, 8, 12, 16, 24, 32, 48, 64, 80, 100, 500, 1000)
# lambda_range = c(100, 500, 1000, 5000, 10000, 50000, 100000)
# eta_range = c(10, 50, 100, 500, 1000, 5000, 10000)
lambda_range = c(100000)
eta_range = c(100, 500, 1000, 5000, 10000)
k_range = c(8,10,12,14,16,18,20,22,24,26,28,30)
# tuning_function(mymtwkm, x, x_p, x_hw, e, group_p, nr, nc, ng, k_range, lambda_range, eta_range, 50, scale_center, scale_std, "mtwkm(p)")






