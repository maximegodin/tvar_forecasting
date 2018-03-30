#-------------------------------------------------------------------------
#  This file is part of the TVAR forecasting application https://github.com/maximegodin/tvar_forecasting
#
#  This application is governed by the MIT license. 
#  You can  use, modify and/ or redistribute this code under the terms
#  of the MIT license:  https://opensource.org/licenses/MIT
#
#  Maxime Godin, Amina Bouchafaa, Ecole Polytechnique
#  March 29th, 2018
#-------------------------------------------------------------------------

library(Rcpp)
library(reshape2)
library(ggplot2)
library(dplyr)
theme_set(theme_bw())

sourceCpp("TVAR_tools.cpp")
sourceCpp("local_yule_walker.cpp")

plot_coefficients <- function(reg_coeffs,order){
  step = .01
  time <- step*c(0:round(1/step))
  df <- as.data.frame(t(rbind(time,sapply(time,function(u) reg_coeffs(u)))))
  colnames(df)[2:(order+1)] <- c(1:order)
  df <- melt(df,id.vars = c("time"),
             variable.name="order", value.name="real")
  ggplot(df)+
    geom_line(aes(x=time, y=real))+
    facet_wrap(~ order, scales = "free_y",ncol=2)+
    xlab("Local time") + ylab("Coefficients") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          strip.text = element_text(size=13))
}

plot_taper <- function(taper){
  step = .001
  x <- 2*step*c(0:round(1/step)) - .5
  df = data.frame(x = x, y = sapply(x,taper))
  ggplot(df)+aes(x=x,y=y)+ geom_line() + xlim(-.5,1.5) +
    xlab("")+ylab("")+theme(plot.background = 
                              element_rect(fill = "transparent", colour = NA),
                            axis.title.y=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks.y=element_blank(),
                            axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank())+
    coord_fixed()
}

plot_sample <- function(sample,id=1){
  df <- data.frame(x = c(0:(nrow(sample)-1)),
                   y = sample[,id])
  ggplot(df)+
    aes(x=x,y=y)+
    geom_line()+
    xlab("Time")+ylab("TVAR sample")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
}

plot_pred_mse_yw <- function(predictions_by_M,best_predictions,current_log_M=0, romberg_MSE = NULL){
  n_obs <- length(best_predictions)
  log_Mmax <- nrow(predictions_by_M)-1
  delta <- matrix(nrow=log_Mmax+1,ncol=n_obs)
  for(i in c(1:n_obs)){
    delta[,i] <- predictions_by_M[,i] - best_predictions[i]
  }
  err <- delta**2
  mse <- rowMeans(err)
  MSE_by_M <- as.data.frame(cbind(c(0:log_Mmax),mse))
  colnames(MSE_by_M) <- c("log_bandwidth","MSE")
  p <- ggplot(MSE_by_M)+
    geom_line(aes(x=log_bandwidth,y=MSE))+
    geom_point(aes(x=log_bandwidth,y=MSE),color="black")+
    geom_point(aes(x=log_bandwidth,y=MSE),data = MSE_by_M[current_log_M+1,],size=3)+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))+
    ylab("Predictor estimation MSE")+xlab("log bandwidth")+scale_y_log10()
  
  if (!is.null(romberg_MSE)){
    p+geom_hline(yintercept=romberg_MSE)
  }
  else{
    p
  }
}


plot_pred_diag_yw <- function(pred_time,sample,predictions,best_predictions){
  
  df <- data.frame(pred = predictions,
                   best_pred = best_predictions,
                   real = sample[pred_time+1,])
  
  ggplot(df)+
    geom_point(aes(x=pred,y=best_pred))+
    ylab("Best predictions") + xlab("Predicted values")+
    geom_abline(slope=1,intercept=0,linetype="dashed")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))

}

plot_romb_pred_diag_yw <- function(pred_time,sample,predictions, romberg_predictions, best_predictions){
  
  df <- data.frame(normal = predictions,
                   romberg = romberg_predictions,
                   best_pred = best_predictions,
                   real = sample[pred_time+1,])
  
  df <- melt(df,id.vars=c("real","best_pred"),value.name="pred",variable.name="method")
  
  ggplot(df)+
    geom_point(aes(x=pred,y=best_pred,color=method))+
    ylab("Best predictions") + xlab("Predicted values")+
    geom_abline(slope=1,intercept=0,linetype="dashed")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.title=element_text(size=12),
          legend.text=element_text(size=12),
          legend.position="bottom")
}

plot_estim_mse_yw <- function(coeffs_by_M, reg_coeffs, current_log_M = 0, MSE_romberg = NULL){
  
  err <- sapply(coeffs_by_M,function(x) sapply(x, function(x) sqrt(sum((x-reg_coeffs)**2))))
  mse <- rowMeans(sapply(coeffs_by_M,function(x) sapply(x, function(x) sum((x-reg_coeffs)**2))))
  log_Mmax <- length(mse)-1
  MSE_by_M <- data.frame(cbind(c(0:log_Mmax),mse))
  colnames(MSE_by_M) <- c("log_bandwidth","MSE")
  p <- ggplot(MSE_by_M)+
    geom_line(aes(x=log_bandwidth,y=MSE))+
    geom_point(aes(x=log_bandwidth,y=MSE))+
    geom_point(aes(x=log_bandwidth,y=MSE),size=3,data = MSE_by_M[current_log_M+1,])+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))+
    ylab("Coefficients estimation MSE")+xlab("log bandwidth")+scale_y_log10()
  
  if (is.null(MSE_romberg)){
    p
  }else{
    p + geom_hline(yintercept = MSE_romberg)
  }
}