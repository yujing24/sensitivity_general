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


Boot_outer <- 400
Boot_inner <- 50


setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm")


# log link, weight 1
{

  search_method <- "adaptive"
  # search_method <- "grid"
  
  # trt
  data_ARC_trt_continuous <- readRDS("ARC_trt_continuous")
  orig_input <- fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_continuous,
                                            Alpha_seq = seq(-3, 3, by = 0.01),
                                            range = c(1.2, 3), # c(1.2, 3) or c(0.2, 0.5)
                                            SIM_model = "norm1coef", # "norm1coef", "fixed_coef", "fixed_band"
                                            link_type = "log", # "log", "logit", "identity"
                                            IF_type = "weight1")
  
  set.seed(123)
  Boot_seed_all <- sample(500000000, size = Boot_outer, replace = T)
  
  start <- Sys.time()
  result_double <- parallel::mclapply(Boot_seed_all, 
                                      f_double_boot_onestep_adaptive,
                                      orig_input = orig_input,
                                      Alpha_seq = seq(-3, 3, by = 0.01),
                                      weight_type = "rademacher", # rademacher, normal, exp
                                      link_type = "log",
                                      IF_type = "weight1",
                                      search_method = search_method, # "adaptive" or "grid" 
                                      L = Boot_inner, # inner bootstrap size
                                      mc.cores = 96)
  end <- Sys.time()
  end - start
  
  saveRDS(result_double, paste0("select_alpha/onestep_adaptive/double_", Boot_outer, "_", Boot_inner,"/double_trt_log_weight1_", search_method, ".RData"))
  
  
  
  # crl
  data_ARC_trt_continuous <- readRDS("ARC_crl_continuous")
  orig_input <- fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_continuous,
                                            Alpha_seq = seq(-3, 3, by = 0.01),
                                            range = c(1.2, 3), # c(1.2, 3) or c(0.2, 0.5)
                                            SIM_model = "norm1coef", # "norm1coef", "fixed_coef", "fixed_band"
                                            link_type = "log", # "log", "logit", "identity"
                                            IF_type = "weight1")
  
  set.seed(234)
  Boot_seed_all <- sample(500000000, size = Boot_outer, replace = T)
  
  start <- Sys.time()
  result_double <- parallel::mclapply(Boot_seed_all, 
                                      f_double_boot_onestep_adaptive,
                                      orig_input = orig_input,
                                      Alpha_seq = seq(-3, 3, by = 0.01),
                                      weight_type = "rademacher", # rademacher, normal, exp
                                      link_type = "log",
                                      IF_type = "weight1",
                                      search_method = search_method, # "adaptive" or "grid" 
                                      L = Boot_inner, # inner bootstrap size
                                      mc.cores = 96)
  end <- Sys.time()
  end - start
  
  saveRDS(result_double, paste0("select_alpha/onestep_adaptive/double_", Boot_outer, "_", Boot_inner,"/double_crl_log_weight1_", search_method, ".RData"))
  
  
}


# log link, weight 2
{
  
  search_method <- "adaptive"
  # search_method <- "grid"
  
  # trt
  data_ARC_trt_continuous <- readRDS("ARC_trt_continuous")
  orig_input <- fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_continuous,
                                            Alpha_seq = seq(-3, 3, by = 0.01),
                                            range = c(1.2, 3), # c(1.2, 3) or c(0.2, 0.5)
                                            SIM_model = "norm1coef", # "norm1coef", "fixed_coef", "fixed_band"
                                            link_type = "log", # "log", "logit", "identity"
                                            IF_type = "weight2")
  
  set.seed(345)
  Boot_seed_all <- sample(500000000, size = Boot_outer, replace = T)
  
  start <- Sys.time()
  result_double <- parallel::mclapply(Boot_seed_all, 
                                      f_double_boot_onestep_adaptive,
                                      orig_input = orig_input,
                                      Alpha_seq = seq(-3, 3, by = 0.01),
                                      weight_type = "rademacher", # rademacher, normal, exp
                                      link_type = "log",
                                      IF_type = "weight2",
                                      search_method = search_method, # "adaptive" or "grid" 
                                      L = Boot_inner, # inner bootstrap size
                                      mc.cores = 96)
  end <- Sys.time()
  end - start
  
  saveRDS(result_double, paste0("select_alpha/onestep_adaptive/double_", Boot_outer, "_", Boot_inner,"/double_trt_log_weight2_", search_method, ".RData"))
  
  
  
  # crl
  data_ARC_trt_continuous <- readRDS("ARC_crl_continuous")
  orig_input <- fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_continuous,
                                            Alpha_seq = seq(-3, 3, by = 0.01),
                                            range = c(1.2, 3), # c(1.2, 3) or c(0.2, 0.5)
                                            SIM_model = "norm1coef", # "norm1coef", "fixed_coef", "fixed_band"
                                            link_type = "log", # "log", "logit", "identity"
                                            IF_type = "weight2")
  
  set.seed(456)
  Boot_seed_all <- sample(500000000, size = Boot_outer, replace = T)
  
  start <- Sys.time()
  result_double <- parallel::mclapply(Boot_seed_all, 
                                      f_double_boot_onestep_adaptive,
                                      orig_input = orig_input,
                                      Alpha_seq = seq(-3, 3, by = 0.01),
                                      weight_type = "rademacher", # rademacher, normal, exp
                                      link_type = "log",
                                      IF_type = "weight2",
                                      search_method = search_method, # "adaptive" or "grid" 
                                      L = Boot_inner, # inner bootstrap size
                                      mc.cores = 96)
  end <- Sys.time()
  end - start
  
  saveRDS(result_double, paste0("select_alpha/onestep_adaptive/double_", Boot_outer, "_", Boot_inner,"/double_crl_log_weight2_", search_method, ".RData"))
  
  
}


# logit link, weight 1
{
  
  search_method <- "adaptive"
  # search_method <- "grid"
  
  
  # trt
  data_ARC_trt_continuous <- readRDS("ARC_trt_binary")
  orig_input <- fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_continuous,
                                            link_type = "logit", # "log", "logit", "identity"
                                            SIM_model = "norm1coef", # "norm1coef", "fixed_coef", "fixed_band"
                                            Alpha_seq = seq(-3, 3, by = 0.01),
                                            range = c(0.2, 0.5), # c(1.2, 3) or c(0.2, 0.5)
                                            IF_type = "weight1")
  
  set.seed(123)
  Boot_seed_all <- sample(500000000, size = Boot_outer, replace = T)
  
  start <- Sys.time()
  result_double <- parallel::mclapply(Boot_seed_all, 
                                      f_double_boot_onestep_adaptive,
                                      orig_input = orig_input,
                                      Alpha_seq = seq(-3, 3, by = 0.01),
                                      range = c(0.2, 0.5), # c(1.2, 3) or c(0.2, 0.5)
                                      weight_type = "rademacher", # rademacher, normal, exp
                                      link_type = "logit",
                                      IF_type = "weight1",
                                      search_method = search_method, # "adaptive" or "grid" 
                                      L = Boot_inner, # inner bootstrap size
                                      mc.cores = 96)
  end <- Sys.time()
  end - start
  
  saveRDS(result_double, paste0("select_alpha/onestep_adaptive/double_", Boot_outer, "_", Boot_inner,"/double_trt_logit_weight1_", search_method, ".RData"))
  
  
  
  # crl
  data_ARC_trt_continuous <- readRDS("ARC_crl_binary")
  orig_input <- fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_continuous,
                                            Alpha_seq = seq(-3, 3, by = 0.01),
                                            SIM_model = "norm1coef", # "norm1coef", "fixed_coef", "fixed_band"
                                            link_type = "logit", # "log", "logit", "identity"
                                            range = c(0.2, 0.5), # c(1.2, 3) or c(0.2, 0.5)
                                            IF_type = "weight1")
  
  set.seed(234)
  Boot_seed_all <- sample(500000000, size = Boot_outer, replace = T)
  
  start <- Sys.time()
  result_double <- parallel::mclapply(Boot_seed_all, 
                                      f_double_boot_onestep_adaptive,
                                      orig_input = orig_input,
                                      Alpha_seq = seq(-3, 3, by = 0.01),
                                      range = c(0.2, 0.5), # c(1.2, 3) or c(0.2, 0.5)
                                      weight_type = "rademacher", # rademacher, normal, exp
                                      link_type = "logit",
                                      IF_type = "weight1",
                                      search_method = search_method, # "adaptive" or "grid" 
                                      L = Boot_inner, # inner bootstrap size
                                      mc.cores = 96)
  end <- Sys.time()
  end - start
  
  saveRDS(result_double, paste0("select_alpha/onestep_adaptive/double_", Boot_outer, "_", Boot_inner,"/double_crl_logit_weight1_", search_method, ".RData"))
  
  
}


# logit link, weight 2
{
  
  search_method <- "adaptive"
  # search_method <- "grid"
  
  # trt
  data_ARC_trt_continuous <- readRDS("ARC_trt_binary")
  orig_input <- fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_continuous,
                                            link_type = "logit", # "log", "logit", "identity"
                                            SIM_model = "norm1coef", # "norm1coef", "fixed_coef", "fixed_band"
                                            Alpha_seq = seq(-3, 3, by = 0.01),
                                            range = c(0.2, 0.5), # c(1.2, 3) or c(0.2, 0.5)
                                            IF_type = "weight2")
  
  set.seed(123)
  Boot_seed_all <- sample(500000000, size = Boot_outer, replace = T)
  
  start <- Sys.time()
  result_double <- parallel::mclapply(Boot_seed_all, 
                                      f_double_boot_onestep_adaptive,
                                      orig_input = orig_input,
                                      Alpha_seq = seq(-3, 3, by = 0.01),
                                      range = c(0.2, 0.5), # c(1.2, 3) or c(0.2, 0.5)
                                      weight_type = "rademacher", # rademacher, normal, exp
                                      link_type = "logit",
                                      IF_type = "weight2",
                                      search_method = search_method, # "adaptive" or "grid" 
                                      L = Boot_inner, # inner bootstrap size
                                      mc.cores = 96)
  end <- Sys.time()
  end - start
  
  saveRDS(result_double, paste0("select_alpha/onestep_adaptive/double_", Boot_outer, "_", Boot_inner,"/double_trt_logit_weight2_", search_method, ".RData"))
  
  
  
  # crl
  data_ARC_trt_continuous <- readRDS("ARC_crl_binary")
  orig_input <- fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_continuous,
                                            link_type = "logit", # "log", "logit", "identity"
                                            SIM_model = "norm1coef", # "norm1coef", "fixed_coef", "fixed_band"
                                            Alpha_seq = seq(-3, 3, by = 0.01),
                                            range = c(0.2, 0.5), # c(1.2, 3) or c(0.2, 0.5)
                                            IF_type = "weight2")
  
  set.seed(234)
  Boot_seed_all <- sample(500000000, size = Boot_outer, replace = T)
  
  start <- Sys.time()
  result_double <- parallel::mclapply(Boot_seed_all, 
                                      f_double_boot_onestep_adaptive,
                                      orig_input = orig_input,
                                      Alpha_seq = seq(-3, 3, by = 0.01),
                                      range = c(0.2, 0.5), # c(1.2, 3) or c(0.2, 0.5)
                                      weight_type = "rademacher", # rademacher, normal, exp
                                      link_type = "logit",
                                      IF_type = "weight2",
                                      search_method = search_method, # "adaptive" or "grid" 
                                      L = Boot_inner, # inner bootstrap size
                                      mc.cores = 96)
  end <- Sys.time()
  end - start
  
  saveRDS(result_double, paste0("select_alpha/onestep_adaptive/double_", Boot_outer, "_", Boot_inner,"/double_crl_logit_weight2_", search_method, ".RData"))
  
  
}



# identity 
{
  
  search_method <- "adaptive"
  # search_method <- "grid"
  
  # trt
  data_ARC_trt_continuous <- readRDS("ARC_trt_continuous")
  orig_input <- fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_continuous,
                                            Alpha_seq = seq(-3, 3, by = 0.01),
                                            range = c(1.2, 3), # c(1.2, 3) or c(0.2, 0.5)
                                            SIM_model = "norm1coef", # "norm1coef", "fixed_coef", "fixed_band"
                                            link_type = "identity" # "log", "logit", "identity"
                                            )
  
  set.seed(567)
  Boot_seed_all <- sample(500000000, size = Boot_outer, replace = T)
  
  start <- Sys.time()
  result_double <- parallel::mclapply(Boot_seed_all, 
                                      f_double_boot_onestep_adaptive,
                                      orig_input = orig_input,
                                      Alpha_seq = seq(-3, 3, by = 0.01),
                                      weight_type = "rademacher", # rademacher, normal, exp
                                      link_type = "identity",
                                      search_method = search_method, # "adaptive" or "grid" 
                                      L = Boot_inner, # inner bootstrap size
                                      mc.cores = 96)
  end <- Sys.time()
  end - start
  
  saveRDS(result_double, paste0("select_alpha/onestep_adaptive/double_", Boot_outer, "_", Boot_inner,"/double_trt_identity_", search_method, ".RData"))
  
  
  
  # crl
  data_ARC_trt_continuous <- readRDS("ARC_crl_continuous")
  orig_input <- fun_est_IF_onestep_analysis(df_sim = data_ARC_trt_continuous,
                                            Alpha_seq = seq(-3, 3, by = 0.01),
                                            range = c(1.2, 3), # c(1.2, 3) or c(0.2, 0.5)
                                            SIM_model = "norm1coef", # "norm1coef", "fixed_coef", "fixed_band"
                                            link_type = "identity" # "log", "logit", "identity"
                                            )
  
  set.seed(678)
  Boot_seed_all <- sample(500000000, size = Boot_outer, replace = T)
  
  start <- Sys.time()
  result_double <- parallel::mclapply(Boot_seed_all, 
                                      f_double_boot_onestep_adaptive,
                                      orig_input = orig_input,
                                      Alpha_seq = seq(-3, 3, by = 0.01),
                                      weight_type = "rademacher", # rademacher, normal, exp
                                      link_type = "identity",
                                      search_method = search_method, # "adaptive" or "grid" 
                                      L = Boot_inner, # inner bootstrap size
                                      mc.cores = 96)
  
  end <- Sys.time()
  end - start
  
  saveRDS(result_double, paste0("select_alpha/onestep_adaptive/double_", Boot_outer, "_", Boot_inner,"/double_crl_identity_", search_method, ".RData"))
  
  
}
