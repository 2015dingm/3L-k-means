inv_scale_function <- function(x){
  inv_scale_x <- matrix(nrow = dim(x)[1], ncol = dim(x)[2])
  for (c in seq(1, dim(x)[2])){
    for (r in seq(1, dim(x)[1])){
      inv_scale_x[r, c] = x[r, c] * scale_std[c] + scale_center[c]
    }
  }
  return(inv_scale_x)
}

aggregate_function <- function(x, group){
  # inv_scale_x = t(inv_scale_function(x))
  group = as.matrix(group)
  temp = cbind(t(x), group)
  x_g = t(aggregate(temp[,seq(1, dim(temp)[2]-1)],list(temp[,dim(temp)[2]]),sum))
  x_g <- x_g[-1,]
  
  return(x_g)
}

output_function <- function(result, file_name){
  write.csv(result$featureWeight, file=paste("../Result/new/", file_name, "_fw.csv", sep = ''))
  write.csv(result$groupWeight, file=paste("../Result/new/", file_name, "_gw.csv", sep = ''))
  write.csv(result$cluster, file=paste("../Result/new/", file_name, "_cluster.csv", sep = ''))

  centers_scaled = result$centers
  centers_noscaled <- inv_scale_function(centers_scaled)
  write.csv(centers_noscaled, file=paste("../Result/new/", file_name, "_center.csv", sep = ''))
}

cost_function <- function(result, x_scale, e_scale=dum, group, nr, nc, ng, k, lambda, eta){
  sum1 = 0
  sum2 = 0
  sum3 = 0
  for (i in seq(1,nr)){
    for (j in seq(1,nc)){
      sum1 = sum1 + e_scale[i] * result$groupWeight[group[j] * k + result$cluster[i]] * result$featureWeight[(j-1) * k + result$cluster[i]] * (x_scale[(j-1) * nr + i] - result$centers[(j-1) * k + result$cluster[i]])^2
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

cost_function_twkm <- function(result, x_scale, e_scale=dum, group, nr, nc, ng, k, lambda, eta){
  sum1 = 0
  sum2 = 0
  sum3 = 0
  for (i in seq(1,nr)){
    for (j in seq(1,nc)){
      sum1 = sum1 + e_scale[i] * result$groupWeight[group[j]+1] * result$featureWeight[j] * (x_scale[(j-1) * nr + i] - result$centers[(j-1) * k + result$cluster[i]])^2
    }
  }
  
  for (t in seq(1,ng)){
    sum2 = sum2 + result$groupWeight[t] * log(result$groupWeight[t])
  }
  
  for (j in seq(1,nc)){
    sum3 = sum3 + result$featureWeight[j] * log(result$featureWeight[j])
  }

  sum2 = sum2 * lambda
  sum3 = sum3 * eta
  
  cost = c(sum1, sum2, sum3)
  
  return(cost)
}

entropy_function <- function(size){
  res <- 0
  for(i in 1:length(size))
  {
    if(size[i]!=0)
      res <- res + size[i]*log(size[i])
  }
  return (-res)
}

distortion_function <- function(result, center, x, e, nr){
  sum = 0
  tot = 0
  for (i in seq(1, nr)){
    distortion = dist(rbind(x[i,], center[result$cluster[i],]), p=2)
    distortion = distortion * e[i,1]
    sum = sum + distortion
    tot = tot + e[i,1]
  }
  md = sum/tot
  
  return(md)
}

evaluation_function <- function(result, x, x_p, x_hw, e, ceiling, nr){
  k_real = length(result$size)
  
  # entropy
  size = result$size/nr
  entropy = entropy_function(size)
  
  # balance
  max = 0
  for (i in seq(1, k_real)){
    if (result$size[i] > max) {
      max <- result$size[i]
    }
  }
  balance = (nr-max)/nr
  
  # complexity
  complexity = (ceiling - k_real)/ceiling
  
  # distortion
  center = inv_scale_function(result$centers)
  center_p = aggregate_function(center, group_p)
  center_hw = aggregate_function(center, group_hw)
  
  mdph = distortion_function(result, center, x, e, nr)
  mdh = distortion_function(result, center_hw, x_hw, e, nr)
  mdp = distortion_function(result, center_p, x_p, e, nr)
  
  evaluation = c(entropy, balance, complexity, mdph, mdh, mdp)
  
  return(evaluation)
}

tuning_function <- function(x, x_scale, x_p, x_hw, e, e_scale, group, nr, nc, ng, k_range, lambda_range, eta_range, ceiling, file_name){
  df <- data.frame(lambda = numeric(0), eta = numeric(0), k = numeric(0), k_real = numeric(0), size = character(0), 
                   var_fw = numeric(0), var_gw = numeric(0), sum1 = numeric(0), sum2 = numeric(0), sum3 = numeric(0), 
                   entropy = numeric(0), balance = numeric(0), complexity = numeric(0), mdph = numeric(0), mdh = numeric(0), mdp = numeric(0), 
                   stringsAsFactors=FALSE)
  index = 1

  for (lam in lambda_range){
    for (et in eta_range){
      for (k in k_range){
        # init_centers <- kmeanspp(x, k, start = "normal")
        # mymfgkm <- mfgkm(x, e, init_centers$inicial.centers, group_p, lambda = lam, eta = et, maxiter = 500, delta = 1e-10, seed = 22)
        result <- mtwkm(x_scale, e_scale, k, group, lambda = lam, eta = et, maxiter = 500, delta = 1e-15, seed = 22)
        save(result, file = paste("../Result/new/model/", lam, '_', et, '_', k,  ".RData", sep = ''))
        fw = as.vector(result$featureWeight)
        gw = as.vector(result$groupWeight)
        k_real = length(result$size)
        size = ''
        for (i in seq(k_real)){
          size = paste(size, as.character(result$size[i]), sep = ',')
        }
        cost = cost_function_twkm(result, x_scale, e_scale, group, nr, nc, ng, k, lam, et)
        evaluation = evaluation_function(result, x, x_p, x_hw, e, ceiling, nr)
        
        df[index,] = c(lam, et, k, k_real, size, var(fw), var(gw), 
                       cost[1], cost[2], cost[3], 
                       evaluation[1], evaluation[2], evaluation[3], evaluation[4], evaluation[5], evaluation[5])
        print(c('lambda=', lam, 'eta=', et, 'k=', k, 'size=', size, 
                'var(fg)=', var(fw), 'var(gw)=', var(gw), 'sum1=', cost[1], 'sum2=', cost[2], 'sum3=', cost[3], 
                'entropy=', evaluation[1], 'balance=', evaluation[2], 'complexity=', evaluation[3], 
                'mdph=', evaluation[4], 'mdh=', evaluation[5], 'mdp=', evaluation[6]))
        
        index = index + 1
      }
    }
  }
  write.csv(df, file=paste("../Result/new/tune/tuning4_", file_name, ".csv", sep = ''))
}
