rm(list=ls())

library(dplyr)
library(survival)
library(MAVE)
library(dfoptim)
library(orthogonalsplinebasis) 
library(BB)
library(statmod)  # For Gauss-Legendre quadrature
library(parallel)

select <- dplyr::select


source("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm/1. function_inference.R")


setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm")


#### CI_plot
{
  
  ### Log link ###
  {
    
    # treatment
    {
     
      # log, weight1 
      data_ARC_trt_continuous <- readRDS("ARC_trt_continuous")
      
      set.seed(123) # for treatment group
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_count_paraboot,
                                        SIM_model = "norm1coef",
                                        data = data_ARC_trt_continuous,
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_trt <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "log", # "log", "logit", "identity"
                                       IF_type = "weight1",
                                       SIM_model = "norm1coef", 
                                       Alpha_seq = c(-1.2, -0.8, -0.4, 0, 0.4, 0.8), # for log
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       mc.cores = 96)
      
      saveRDS(result_trt, "result/CI_plot/paraboot_trt_log_weight1")
      
      
      # log, weight2
      data_ARC_trt_continuous <- readRDS("ARC_trt_continuous")
      
      set.seed(234) # for treatment group
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_count_paraboot,
                                        data = data_ARC_trt_continuous,
                                        SIM_model = "norm1coef",
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_trt <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "log", # "log", "logit", "identity"
                                       IF_type = "weight2",
                                       SIM_model = "norm1coef", 
                                       Alpha_seq = c(-1.2, -0.8, -0.4, 0, 0.4, 0.8), # for log
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       mc.cores = 96)
      
      saveRDS(result_trt, "result/CI_plot/paraboot_trt_log_weight2")
      
    }
    
    
    # control
    {
      
      # log, weight1 
      data_ARC_crl_continuous <- readRDS("ARC_crl_continuous")
      
      set.seed(345) 
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_count_paraboot,
                                        data = data_ARC_crl_continuous,
                                        SIM_model = "norm1coef",
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_crl <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "log", # "log", "logit", "identity"
                                       IF_type = "weight1",
                                       SIM_model = "norm1coef", 
                                       Alpha_seq = c(-1.2, -0.8, -0.4, 0, 0.4, 0.8), # for log
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       mc.cores = 96)
      
      saveRDS(result_crl, "result/CI_plot/paraboot_crl_log_weight1")
      
      
      # log, weight2
      data_ARC_crl_continuous <- readRDS("ARC_crl_continuous")
      
      set.seed(456) 
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_count_paraboot,
                                        data = data_ARC_crl_continuous,
                                        SIM_model = "norm1coef",
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_crl <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "log", # "log", "logit", "identity"
                                       IF_type = "weight2",
                                       SIM_model = "norm1coef", 
                                       Alpha_seq = c(-1.2, -0.8, -0.4, 0, 0.4, 0.8), # for log
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       mc.cores = 96)
      
      saveRDS(result_crl, "result/CI_plot/paraboot_crl_log_weight2")
      
    }
    
    
  }
  
  
  ### Identity link ###
  {
    
    # treatment
    {
      
      # identity, weight1 
      data_ARC_trt_continuous <- readRDS("ARC_trt_continuous")
      
      set.seed(123) # for treatment group
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_count_paraboot,
                                        data = data_ARC_trt_continuous,
                                        SIM_model = "norm1coef",
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_trt <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "identity", # "log", "logit", "identity"
                                       IF_type = "weight1",
                                       SIM_model = "norm1coef",
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       mc.cores = 96)
      
      saveRDS(result_trt, "result/CI_plot/paraboot_trt_identity")
      
      
    }
    
    
    # control
    {
      
      # identity, weight1 
      data_ARC_crl_continuous <- readRDS("ARC_crl_continuous")
      
      set.seed(345) 
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_count_paraboot,
                                        data = data_ARC_crl_continuous,
                                        SIM_model = "norm1coef",
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_crl <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "identity", # "log", "logit", "identity"
                                       SIM_model = "norm1coef",
                                       IF_type = "weight1",
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       mc.cores = 96)
      
      saveRDS(result_crl, "result/CI_plot/paraboot_crl_identity")
      
      
    }
    
    
  }
  
  
  ### Logit link ###
  {
    
    # treatment
    {
      
      # logit, weight1 
      data_ARC_trt_binary <- readRDS("ARC_trt_binary")
      
      set.seed(123) # for treatment group
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_binary_paraboot,
                                        data = data_ARC_trt_binary,
                                        SIM_model = "norm1coef", 
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_trt <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "logit", # "log", "logit", "identity"
                                       SIM_model = "norm1coef", 
                                       Alpha_seq = c(-0.8, -0.4, 0, 0.4, 0.8, 1.2), # for binary
                                       IF_type = "weight1",
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       mc.cores = 96)
      
      saveRDS(result_trt, "result/CI_plot/paraboot_trt_logit_weight1")
      
      
      # logit, weight2
      data_ARC_trt_binary <- readRDS("ARC_trt_binary")
      
      set.seed(234) # for treatment group
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_binary_paraboot,
                                        data = data_ARC_trt_binary,
                                        SIM_model = "norm1coef", 
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_trt <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "logit", # "log", "logit", "identity"
                                       SIM_model = "norm1coef", 
                                       Alpha_seq = c(-0.8, -0.4, 0, 0.4, 0.8, 1.2), # for binary
                                       IF_type = "weight2",
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       mc.cores = 96)
      
      saveRDS(result_trt, "result/CI_plot/paraboot_trt_logit_weight2")
      
    }
    
    
    # control
    {
      
      # logit, weight1 
      data_ARC_crl_binary <- readRDS("ARC_crl_binary")
      
      set.seed(345) 
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_binary_paraboot,
                                        data = data_ARC_crl_binary,
                                        SIM_model = "norm1coef", 
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_crl <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "logit", # "log", "logit", "identity"
                                       SIM_model = "norm1coef", 
                                       Alpha_seq = c(-0.8, -0.4, 0, 0.4, 0.8, 1.2), # for binary
                                       IF_type = "weight1",
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       mc.cores = 96)
      
      saveRDS(result_crl, "result/CI_plot/paraboot_crl_logit_weight1")
      
      
      # logit, weight2
      data_ARC_crl_binary <- readRDS("ARC_crl_binary")
      
      set.seed(456) 
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_binary_paraboot,
                                        data = data_ARC_crl_binary,
                                        SIM_model = "norm1coef", 
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_crl <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "logit", # "log", "logit", "identity"
                                       SIM_model = "norm1coef", 
                                       Alpha_seq = c(-0.8, -0.4, 0, 0.4, 0.8, 1.2), # for binary
                                       IF_type = "weight2",
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       mc.cores = 96)
      
      saveRDS(result_crl, "result/CI_plot/paraboot_crl_logit_weight2")
      
    }
    
    
  }
  
}


#### Contour_plot
{
  
  ### Log link ###
  {
    
    # treatment
    {
      
      # log, weight1 
      data_ARC_trt_continuous <- readRDS("ARC_trt_continuous")
      
      set.seed(123) # for treatment group
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_count_paraboot,
                                        data = data_ARC_trt_continuous,
                                        SIM_model = "norm1coef", 
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_trt <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "log", # "log", "logit", "identity"
                                       IF_type = "weight1",
                                       SIM_model = "norm1coef", 
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       Alpha_seq = seq(-1.2, 0.8, by = 0.02),# for log
                                       mc.cores = 96)
      
      saveRDS(result_trt, "result/contour_plot/paraboot_trt_log_weight1")
      
      
      # log, weight2
      data_ARC_trt_continuous <- readRDS("ARC_trt_continuous")
      
      set.seed(234) # for treatment group
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_count_paraboot,
                                        data = data_ARC_trt_continuous,
                                        SIM_model = "norm1coef", 
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_trt <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "log", # "log", "logit", "identity"
                                       IF_type = "weight2",
                                       SIM_model = "norm1coef", 
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       Alpha_seq = seq(-1.2, 0.8, by = 0.02), # for log
                                       mc.cores = 96)
      
      saveRDS(result_trt, "result/contour_plot/paraboot_trt_log_weight2")
      
    }
    
    
    # control
    {
      
      # log, weight1 
      data_ARC_crl_continuous <- readRDS("ARC_crl_continuous")
      
      set.seed(345) 
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_count_paraboot,
                                        data = data_ARC_crl_continuous,
                                        SIM_model = "norm1coef", 
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_crl <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "log", # "log", "logit", "identity"
                                       IF_type = "weight1",
                                       SIM_model = "norm1coef", 
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       Alpha_seq = seq(-1.2, 0.8, by = 0.02), # for log
                                       mc.cores = 96)
      
      saveRDS(result_crl, "result/contour_plot/paraboot_crl_log_weight1")
      
      
      # log, weight2
      data_ARC_crl_continuous <- readRDS("ARC_crl_continuous")
      
      set.seed(456) 
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_count_paraboot,
                                        data = data_ARC_crl_continuous,
                                        SIM_model = "norm1coef", 
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_crl <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "log", # "log", "logit", "identity"
                                       IF_type = "weight2",
                                       SIM_model = "norm1coef", 
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       Alpha_seq = seq(-1.2, 0.8, by = 0.02), # for log
                                       mc.cores = 96)
      
      saveRDS(result_crl, "result/contour_plot/paraboot_crl_log_weight2")
      
    }
    
    
  }
  
  
  ### Identity link ###
  {
    
    # treatment
    {
      
      # identity, weight1 
      data_ARC_trt_continuous <- readRDS("ARC_trt_continuous")
      
      set.seed(123) # for treatment group
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_count_paraboot,
                                        data = data_ARC_trt_continuous,
                                        SIM_model = "norm1coef", 
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_trt <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "identity", # "log", "logit", "identity"
                                       IF_type = "weight1",
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       Alpha_seq = seq(-1.2, 0.8, by = 0.02),
                                       mc.cores = 96)
      
      saveRDS(result_trt, "result/contour_plot/paraboot_trt_identity")
      
      
    }
    
    
    # control
    {
      
      # identity, weight1 
      data_ARC_crl_continuous <- readRDS("ARC_crl_continuous")
      
      set.seed(345) 
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_count_paraboot,
                                        data = data_ARC_crl_continuous,
                                        SIM_model = "norm1coef", 
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_crl <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "identity", # "log", "logit", "identity"
                                       IF_type = "weight1",
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       Alpha_seq = seq(-1.2, 0.8, by = 0.02),
                                       mc.cores = 96)
      
      saveRDS(result_crl, "result/contour_plot/paraboot_crl_identity")
      
      
    }
    
    
  }
  
  
  ### Logit link ###
  {
    
    # treatment
    {
      
      # log, weight1 
      data_ARC_trt_binary <- readRDS("ARC_trt_binary")
      
      set.seed(123) # for treatment group
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_binary_paraboot,
                                        data = data_ARC_trt_binary,
                                        SIM_model = "norm1coef", 
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_trt <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "logit", # "log", "logit", "identity"
                                       IF_type = "weight1",
                                       SIM_model = "norm1coef", 
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       Alpha_seq = seq(-0.8, 1.2, by = 0.02),  # for binary
                                       mc.cores = 96)
      
      saveRDS(result_trt, "result/contour_plot/paraboot_trt_logit_weight1")
      
      
      # log, weight2
      data_ARC_trt_binary <- readRDS("ARC_trt_binary")
      
      set.seed(234) # for treatment group
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_binary_paraboot,
                                        data = data_ARC_trt_binary,
                                        SIM_model = "norm1coef", 
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_trt <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "logit", # "log", "logit", "identity"
                                       IF_type = "weight2",
                                       SIM_model = "norm1coef", 
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       Alpha_seq = seq(-0.8, 1.2, by = 0.02),  # for binary
                                       mc.cores = 96)
      
      saveRDS(result_trt, "result/contour_plot/paraboot_trt_logit_weight2")
      
    }
    
    
    # control
    {
      
      # log, weight1 
      data_ARC_crl_binary <- readRDS("ARC_crl_binary")
      
      set.seed(345) 
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_binary_paraboot,
                                        data = data_ARC_crl_binary,
                                        SIM_model = "norm1coef", 
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_crl <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "logit", # "log", "logit", "identity"
                                       IF_type = "weight1",
                                       SIM_model = "norm1coef", 
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       Alpha_seq = seq(-0.8, 1.2, by = 0.02),  # for binary
                                       mc.cores = 96)
      
      saveRDS(result_crl, "result/contour_plot/paraboot_crl_logit_weight1")
      
      
      # log, weight2
      data_ARC_crl_binary <- readRDS("ARC_crl_binary")
      
      set.seed(456) 
      B <- 200
      Index_sim <- sample(1:50000000, B, replace = F)
      
      result_simu <- parallel::mclapply(Index_sim,
                                        inner_ARC_SIR_sim_data_binary_paraboot,
                                        data = data_ARC_crl_binary,
                                        SIM_model = "norm1coef", 
                                        mc.cores = 96)
      
      call_function_with_data <- function(index, fun, ...){
        
        data <- result_simu[[index]][[1]] # for the original data
        fun(data, ...)
      }
      
      # Point estimation for each parametric bootstrap dataset
      Index_est <- seq(1, B, 1)
      
      result_crl <- parallel::mclapply(Index_est, 
                                       call_function_with_data,
                                       f_data_analysis_line,
                                       link_type = "logit", # "log", "logit", "identity"
                                       IF_type = "weight2",
                                       SIM_model = "norm1coef", 
                                       spline_seq_orig = c(90, 180, 270, 360),
                                       Alpha_seq = seq(-0.8, 1.2, by = 0.02),  # for binary
                                       mc.cores = 96)
      
      saveRDS(result_crl, "result/contour_plot/paraboot_crl_logit_weight2")
      
    }
    
    
  }
  
}




