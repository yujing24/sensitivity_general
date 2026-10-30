########## generate the bootstrap t confidence interval ##########
# load the bootstrap dataset
rm(list=ls())
library(dplyr)

source("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm/1. function_inference.R")
setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm")

Alpha_seq = seq(-1.2, 0.8, by = 0.02)
A <- length(Alpha_seq)
spline_seq = seq(60, 460, by = 1)
range_value = c(1.2, 3)


## log, weight1 
{

  # load the parametric bootstrap estimation results
  data_trt <- readRDS("result/contour_plot/paraboot_trt_log_weight1")
  data_crl <- readRDS("result/contour_plot/paraboot_crl_log_weight1")
  
  # load the estimation based on the original synthetic data
  SIM_treat <- readRDS("result/contour_plot/result_trt_log_weight1.RData")
  SIM_control <- readRDS("result/contour_plot/result_crl_log_weight1.RData")
  SIM_treat <- SIM_treat[[1]]
  SIM_control <- SIM_control[[1]]
  
  
  ##############   contour plot for SIM - Month 6 - data processing ############## 
  {
    
    
    Estimate_vec <- rep(NA, A^2)
    Lower_ci <- rep(NA, A^2)
    Upper_ci <- rep(NA, A^2)
    Alpha_trt <- rep(NA, A^2)
    Alpha_crl <- rep(NA, A^2)
    
    for(a in 1:A){
      for(b in 1:A){
        
        res <- Wald_CI_effect(alpha_trt =  Alpha_seq[a],
                               alpha_crl =  Alpha_seq[b],
                               month = 6)
        
        h <- (a-1)*A + b
        
        Estimate_vec[h] <- res[1]
        
        Lower_ci[h] <- res[2]
        Upper_ci[h] <- res[3]
        
        Alpha_trt[h] <- Alpha_seq[a]
        Alpha_crl[h] <- Alpha_seq[b]
        
      }
      print(a)
    }
    
    plot_df <- data.frame(Alpha_trt,
                          Alpha_crl,
                          Estimate_vec, 
                          Lower_ci, 
                          Upper_ci)
    head(plot_df)
    
    saveRDS(plot_df, "result/contour_plot/Contour_plot_month6_SIM.rds")
    
    
  }
  
  
  ##############   contour plot for SIM - Month 12 - data processing ############## 
  {
    
    
    Estimate_vec <- rep(NA, A^2)
    Lower_ci <- rep(NA, A^2)
    Upper_ci <- rep(NA, A^2)
    Alpha_trt <- rep(NA, A^2)
    Alpha_crl <- rep(NA, A^2)
    
    for(a in 1:A){
      for(b in 1:A){
        
        res <- Wald_CI_effect(alpha_trt =  Alpha_seq[a],
                                alpha_crl =  Alpha_seq[b],
                                month = 12)
        
        h <- (a-1)*A + b
        
        Estimate_vec[h] <- res[1]
        
        Lower_ci[h] <- res[2]
        Upper_ci[h] <- res[3]
        
        Alpha_trt[h] <- Alpha_seq[a]
        Alpha_crl[h] <- Alpha_seq[b]
        
      }
      print(a)
    }
    
    plot_df <- data.frame(Alpha_trt,
                          Alpha_crl,
                          Estimate_vec, 
                          Lower_ci, 
                          Upper_ci)
    head(plot_df)
    
    saveRDS(plot_df, "result/contour_plot/Contour_plot_month12_SIM.rds")
    
    
  }
  
  
}


## log, weight2
{
  
  # load the parametric bootstrap estimation results
  data_trt <- readRDS("result/contour_plot/paraboot_trt_log_weight2")
  data_crl <- readRDS("result/contour_plot/paraboot_crl_log_weight2")
  
  # load the estimation based on the original synthetic data
  SIM_treat <- readRDS("result/contour_plot/result_trt_log_weight2.RData")
  SIM_control <- readRDS("result/contour_plot/result_crl_log_weight2.RData")
  SIM_treat <- SIM_treat[[1]]
  SIM_control <- SIM_control[[1]]
  
  
  ##############   contour plot for SIM - Month 6 - data processing ############## 
  {
    
    
    Estimate_vec <- rep(NA, A^2)
    Lower_ci <- rep(NA, A^2)
    Upper_ci <- rep(NA, A^2)
    Alpha_trt <- rep(NA, A^2)
    Alpha_crl <- rep(NA, A^2)
    
    for(a in 1:A){
      for(b in 1:A){
        
        res <- Wald_CI_effect(alpha_trt =  Alpha_seq[a],
                                alpha_crl =  Alpha_seq[b],
                                month = 6)
        
        h <- (a-1)*A + b
        
        Estimate_vec[h] <- res[1]
        
        Lower_ci[h] <- res[2]
        Upper_ci[h] <- res[3]
        
        Alpha_trt[h] <- Alpha_seq[a]
        Alpha_crl[h] <- Alpha_seq[b]
        
      }
      print(a)
    }
    
    plot_df <- data.frame(Alpha_trt,
                          Alpha_crl,
                          Estimate_vec, 
                          Lower_ci, 
                          Upper_ci)
    head(plot_df)
    
    saveRDS(plot_df, "result/contour_plot/Contour_plot_month6_SIM_weight2.rds")
    
    
  }
  
  
  ##############   contour plot for SIM - Month 12 - data processing ############## 
  {
    
    
    Estimate_vec <- rep(NA, A^2)
    Lower_ci <- rep(NA, A^2)
    Upper_ci <- rep(NA, A^2)
    Alpha_trt <- rep(NA, A^2)
    Alpha_crl <- rep(NA, A^2)
    
    for(a in 1:A){
      for(b in 1:A){
        
        res <- Wald_CI_effect(alpha_trt =  Alpha_seq[a],
                                alpha_crl =  Alpha_seq[b],
                                month = 12)
        
        h <- (a-1)*A + b
        
        Estimate_vec[h] <- res[1]
        
        Lower_ci[h] <- res[2]
        Upper_ci[h] <- res[3]
        
        Alpha_trt[h] <- Alpha_seq[a]
        Alpha_crl[h] <- Alpha_seq[b]
        
      }
      print(a)
    }
    
    plot_df <- data.frame(Alpha_trt,
                          Alpha_crl,
                          Estimate_vec, 
                          Lower_ci, 
                          Upper_ci)
    head(plot_df)
    
    saveRDS(plot_df, "result/contour_plot/Contour_plot_month12_SIM_weight2.rds")
    
    
  }
  
  
}

