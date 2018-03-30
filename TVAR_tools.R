library(Rcpp)

sourceCpp("TVAR_tools.cpp")

predict_TVAR_by_M <- function(sample,t,order,max_log_bandwidth,coeffs){
  n_obs <- ncol(sample)
  res <- matrix(nrow = max_log_bandwidth+1,ncol=n_obs)
  for(i in (1:n_obs)){
    for(j in (1:(max_log_bandwidth+1))){
      res[j,i] <- predict_TVAR(sample[,i],t,coeffs[[i]][[j]])
    }
  }
  return(res)
}

predict_TVAR_romberg <- function(n_obs,predictions_by_M,log_bandwidth_min,log_bandwidth_max){
  k <- log_bandwidth_max-log_bandwidth_min
  weigths <- romberg_weights(k)
  res <- vector(mode="numeric",n_obs)
  for(i in (0:k)){
    res <- res + predictions_by_M[log_bandwidth_min+i+1,]*weigths[i+1]
  }
  return(res)
}

estim_coeffs_romberg <- function(n_obs,coeffs_by_M,log_bandwidth_min,log_bandwidth_max,order){
  k = log_bandwidth_max-log_bandwidth_min
  n_obs <- length(coeffs_by_M)
  weigths = romberg_weights(k)
  res = matrix(0,ncol=n_obs,nrow=order)
  for(j in (1:n_obs)){
    for(i in (0:k)){
      res[,j] <- res[,j] + coeffs_by_M[[j]][[log_bandwidth_min+i+1]]*weigths[i+1]
    }
  }
  return(res)
  
}