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
Boot_num <- 200
setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm/select_alpha/onestep_adaptive")

# log link
{
  
  search_method = "grid"
  # search_method = "adaptive"
  
  # trt, weight1
  orig_trt_log_weight1 <- readRDS("orig_trt_log_weight1.RData")
  
  set.seed(123)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_trt_log_weight1,
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "log", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(1.2, 3),
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_trt_log_weight1_", search_method, ".RData"))
  
  
  # trt, weight2
  orig_trt_log_weight2 <- readRDS("orig_trt_log_weight2.RData")
  
  set.seed(234)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_trt_log_weight2,
                                    # Alpha_seq = seq(-1.5, 1.5, by = 0.01),
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "log", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(1.2, 3),
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_trt_log_weight2_", search_method, ".RData"))
  
  
  # crl, weight1
  orig_crl_log_weight1 <- readRDS("orig_crl_log_weight1.RData")
  
  set.seed(345)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_crl_log_weight1,
                                    # Alpha_seq = seq(-1.5, 1.5, by = 0.01),
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "log", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(1.2, 3),
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_crl_log_weight1_", search_method, ".RData"))
  
  
  # crl, weight2
  orig_crl_log_weight2 <- readRDS("orig_crl_log_weight2.RData")
  
  set.seed(456)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_crl_log_weight2,
                                    # Alpha_seq = seq(-1.5, 1.5, by = 0.01),
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "log", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(1.2, 3),
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_crl_log_weight2_", search_method, ".RData"))
  
  
}

# logit link 
{
  
  # search_method = "grid"
  search_method = "adaptive"
  
  # trt, weight1
  orig_trt_log_weight1 <- readRDS("orig_trt_logit_weight1.RData")
  
  set.seed(123)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_trt_log_weight1,
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "logit", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(0.2, 0.5),
                                    mc.cores = 64)
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_trt_logit_weight1_", search_method, ".RData"))
  
  
  # trt, weight2
  orig_trt_log_weight2 <- readRDS("orig_trt_logit_weight2.RData")
  
  set.seed(234)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_trt_log_weight2,
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "logit", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(0.2, 0.5),
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_trt_logit_weight2_", search_method, ".RData"))
  
  
  # crl, weight1
  orig_crl_log_weight1 <- readRDS("orig_crl_logit_weight1.RData")
  
  set.seed(345)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_crl_log_weight1,
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "logit", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(0.2, 0.5),
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_crl_logit_weight1_", search_method, ".RData"))
  
  
  # crl, weight2
  orig_crl_log_weight2 <- readRDS("orig_crl_logit_weight2.RData")
  
  set.seed(456)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_crl_log_weight2,
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "logit", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(0.2, 0.5),
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_crl_logit_weight2_", search_method, ".RData"))
  
}

# identity link 
{
  
  # search_method = "grid"
  search_method = "adaptive"
  
  
  # trt
  orig_trt_identity <- readRDS("orig_trt_identity.RData")
  
  set.seed(567)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_trt_identity,
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "identity", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(1.2, 3),
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_trt_identity_", search_method, ".RData"))
  
  
  
  # crl
  orig_crl_identity <- readRDS("orig_crl_identity.RData")
  
  set.seed(678)
  Boot_seed_all <- sample(500000000, size = Boot_num, replace = T)
  
  result_bias <- parallel::mclapply(Boot_seed_all, 
                                    f_onestep_bias_adaptive,
                                    orig_input = orig_crl_identity,
                                    Alpha_seq = seq(-3, 3, by = 0.01),
                                    weight_type = "rademacher", # rademacher, normal, exp
                                    link_type = "identity", # "log", "logit", "identity"
                                    search_method = search_method, # "binary", "grid" or "adaptive"
                                    range = c(1.2, 3),
                                    mc.cores = 64)
  
  saveRDS(result_bias, paste0("bias_correct_", Boot_num, "/corr_crl_identity_", search_method, ".RData"))
  
  
  
}





