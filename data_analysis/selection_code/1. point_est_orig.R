rm(list=ls())

library(dplyr)
library(survival)
library(MAVE)
library(dfoptim)
library(parallel)
library(orthogonalsplinebasis) 
library(BB)
library(statmod)  # For Gauss-Legendre quadrature

source("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm/1. function_selection.R")


## onestep - multiplier bootstrap (grid search method)
{
  
  setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm")
  
  
  # log link
  {
   
    # trt, weight1  
    data_ARC_trt_continuous <- readRDS("ARC_trt_continuous")
    result_trt <-  fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_continuous,
                                               Alpha_seq = seq(-3, 3, by = 0.01),
                                               link_type = "log", # "log", "logit", "identity"
                                               SIM_model = "norm1coef",
                                               IF_type = "weight1")
    saveRDS(result_trt, "select_alpha/onestep_adaptive/orig_trt_log_weight1.RData")
    
    # trt, weight2
    data_ARC_trt_continuous <- readRDS("ARC_trt_continuous")
    result_trt <-  fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_continuous,
                                               Alpha_seq = seq(-3, 3, by = 0.01),
                                               link_type = "log", # "log", "logit", "identity"
                                               SIM_model = "norm1coef",
                                               IF_type = "weight2")
    saveRDS(result_trt, "select_alpha/onestep_adaptive/orig_trt_log_weight2.RData")
    
    
    # crl, weight1 
    data_ARC_trt_continuous <- readRDS("ARC_crl_continuous")
    result_trt <-  fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_continuous,
                                               Alpha_seq = seq(-3, 3, by = 0.01),
                                               link_type = "log", # "log", "logit", "identity"
                                               SIM_model = "norm1coef",
                                               IF_type = "weight1")
    saveRDS(result_trt, "select_alpha/onestep_adaptive/orig_crl_log_weight1.RData")
    
    
    # crl, weight2
    data_ARC_trt_continuous <- readRDS("ARC_crl_continuous")
    result_trt <-  fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_continuous,
                                               Alpha_seq = seq(-3, 3, by = 0.01),
                                               link_type = "log", # "log", "logit", "identity"
                                               SIM_model = "norm1coef",
                                               IF_type = "weight2")
    saveRDS(result_trt, "select_alpha/onestep_adaptive/orig_crl_log_weight2.RData")
    
    
  }
  
  
  # identity link
  {
    
    data_ARC_trt_continuous <- readRDS("ARC_trt_continuous")
    result_trt <-  fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_continuous,
                                               Alpha_seq = seq(-3, 3, by = 0.01),
                                               link_type = "identity", # "log", "logit", "identity"
                                               SIM_model = "norm1coef")
    saveRDS(result_trt, "select_alpha/onestep_adaptive/orig_trt_identity.RData")
  
    
    data_ARC_trt_continuous <- readRDS("ARC_crl_continuous")
    result_trt <-  fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_continuous,
                                               Alpha_seq = seq(-3, 3, by = 0.01),
                                               link_type = "identity", # "log", "logit", "identity"
                                               SIM_model = "norm1coef")
    saveRDS(result_trt, "select_alpha/onestep_adaptive/orig_crl_identity.RData")
    
    
  }
  
  
  # logit link 
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm")
    data_ARC_trt_binary <- readRDS("ARC_trt_binary")
    
    result_trt <-  fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_binary,
                                               Alpha_seq = seq(-3, 3, by = 0.01),
                                               link_type = "logit", # "log", "logit", "identity"
                                               SIM_model = "norm1coef",
                                               IF_type = "weight1",
                                               range = c(0.2, 0.5))
    saveRDS(result_trt, "select_alpha/onestep_adaptive/orig_trt_logit_weight1.RData")
    
    
    result_trt <-  fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_binary,
                                               Alpha_seq = seq(-3, 3, by = 0.01),
                                               link_type = "logit", # "log", "logit", "identity"
                                               SIM_model = "norm1coef",
                                               IF_type = "weight2",
                                               range = c(0.2, 0.5))
    saveRDS(result_trt, "select_alpha/onestep_adaptive/orig_trt_logit_weight2.RData")
    
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm")
    data_ARC_trt_binary <- readRDS("ARC_crl_binary")
    
    result_trt <-  fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_binary,
                                               Alpha_seq = seq(-3, 3, by = 0.01),
                                               link_type = "logit", # "log", "logit", "identity"
                                               SIM_model = "norm1coef",
                                               IF_type = "weight1",
                                               range = c(0.2, 0.5))
    saveRDS(result_trt, "select_alpha/onestep_adaptive/orig_crl_logit_weight1.RData")
    
    
    result_trt <-  fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_binary,
                                               Alpha_seq = seq(-3, 3, by = 0.01),
                                               link_type = "logit", # "log", "logit", "identity"
                                               SIM_model = "norm1coef",
                                               IF_type = "weight2",
                                               range = c(0.2, 0.5))
    saveRDS(result_trt, "select_alpha/onestep_adaptive/orig_crl_logit_weight2.RData")
    
    
  }
  

}








