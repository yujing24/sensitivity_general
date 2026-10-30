rm(list=ls())

library(dplyr)
library(survival)
library(MAVE)
library(dfoptim)
library(orthogonalsplinebasis) 
library(BB)
library(statmod)  # For Gauss-Legendre quadrature

select <- dplyr::select



source("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm/1. function_inference.R")


# log and identity
{
  
  # Original point estimation - include IF variance, less alpha
  {
    
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm")
    
    # trt_group
    {
      
      data_ARC_trt_continuous <- readRDS("ARC_trt_continuous")
      
      
      result_trt_log_weight1 <- f_data_analysis_point(data = data_ARC_trt_continuous,
                                                      link_type = "log", # "log", "logit", "identity"
                                                      Alpha_seq = c(-1.2, -0.8, -0.4, 0, 0.4, 0.8), # for log
                                                      IF_type = "weight1",
                                                      SIM_model = "norm1coef"
                                                      )
      saveRDS(result_trt_log_weight1, "result/CI_plot/result_trt_log_weight1.RData")
      
      result_trt_log_weight2 <- f_data_analysis_point(data = data_ARC_trt_continuous,
                                                      link_type = "log", # "log", "logit", "identity"
                                                      Alpha_seq = c(-1.2, -0.8, -0.4, 0, 0.4, 0.8), # for log
                                                      IF_type = "weight2",
                                                      SIM_model = "norm1coef")
      saveRDS(result_trt_log_weight2, "result/CI_plot/result_trt_log_weight2.RData")
      
      result_trt_identity <- f_data_analysis_point(data = data_ARC_trt_continuous,
                                                   link_type = "identity", # "log", "logit", "identity"
                                                   Alpha_seq = c(-1.2, -0.8, -0.4, 0, 0.4, 0.8), # for log
                                                   IF_type = "weight2",
                                                   SIM_model = "norm1coef")
      saveRDS(result_trt_identity, "result/CI_plot/result_trt_identity.RData")
      
    }
    
    
    # crl_group
    {
      
      
      data_ARC_crl_continuous <- readRDS("ARC_crl_continuous")
      
      
      result_crl_log_weight1 <- f_data_analysis_point(data = data_ARC_crl_continuous,
                                                      link_type = "log", # "log", "logit", "identity"
                                                      Alpha_seq = c(-1.2, -0.8, -0.4, 0, 0.4, 0.8), # for log
                                                      IF_type = "weight1",
                                                      SIM_model = "norm1coef")
      saveRDS(result_crl_log_weight1, "result/CI_plot/result_crl_log_weight1.RData")
      
      result_crl_log_weight2 <- f_data_analysis_point(data = data_ARC_crl_continuous,
                                                      link_type = "log", # "log", "logit", "identity"
                                                      Alpha_seq = c(-1.2, -0.8, -0.4, 0, 0.4, 0.8), # for log
                                                      IF_type = "weight2",
                                                      SIM_model = "norm1coef")
      saveRDS(result_crl_log_weight2, "result/CI_plot/result_crl_log_weight2.RData")
      
      result_crl_identity <- f_data_analysis_point(data = data_ARC_crl_continuous,
                                                   link_type = "identity", # "log", "logit", "identity"
                                                   Alpha_seq = c(-1.2, -0.8, -0.4, 0, 0.4, 0.8), # for log
                                                   IF_type = "weight2",
                                                   SIM_model = "norm1coef")
      saveRDS(result_crl_identity, "result/CI_plot/result_crl_identity.RData")
      
    }
    
  }
  
  
  # Original line estimation - exclude IF variance, less alpha
  {
    
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm")
    
    # trt_group
    {
      
      data_ARC_trt_continuous <- readRDS("ARC_trt_continuous")
      
      result_trt_log_weight1 <- f_data_analysis_line(data = data_ARC_trt_continuous,
                                                     link_type = "log", # "log", "logit", "identity"
                                                     Alpha_seq = c(-1.2, -0.8, -0.4, 0, 0.4, 0.8), # for log
                                                     IF_type = "weight1",
                                                     SIM_model = "norm1coef")
      saveRDS(result_trt_log_weight1, "result/line_plot/result_trt_log_weight1.RData")
      
      
      result_trt_log_weight2 <- f_data_analysis_line(data = data_ARC_trt_continuous,
                                                     link_type = "log", # "log", "logit", "identity"
                                                     Alpha_seq = c(-1.2, -0.8, -0.4, 0, 0.4, 0.8), # for log
                                                     IF_type = "weight2",
                                                     SIM_model = "norm1coef")
      saveRDS(result_trt_log_weight2, "result/line_plot/result_trt_log_weight2.RData")
      
      
      result_trt_identity <- f_data_analysis_line(data = data_ARC_trt_continuous,
                                                  link_type = "identity", # "log", "logit", "identity"
                                                  Alpha_seq = c(-1.2, -0.8, -0.4, 0, 0.4, 0.8), # for log
                                                  IF_type = "weight2",
                                                  SIM_model = "norm1coef")
      saveRDS(result_trt_identity, "result/line_plot/result_trt_identity.RData")
      
    }
    
    # crl_group
    {
      
      data_ARC_crl_continuous <- readRDS("ARC_crl_continuous")
      
      result_crl_log_weight1 <- f_data_analysis_line(data = data_ARC_crl_continuous,
                                                     link_type = "log", # "log", "logit", "identity"
                                                     Alpha_seq = c(-1.2, -0.8, -0.4, 0, 0.4, 0.8), # for log
                                                     IF_type = "weight1",
                                                     SIM_model = "norm1coef")
      saveRDS(result_crl_log_weight1, "result/line_plot/result_crl_log_weight1.RData")
      
      result_crl_log_weight2 <- f_data_analysis_line(data = data_ARC_crl_continuous,
                                                     link_type = "log", # "log", "logit", "identity"
                                                     Alpha_seq = c(-1.2, -0.8, -0.4, 0, 0.4, 0.8), # for log
                                                     IF_type = "weight2",
                                                     SIM_model = "norm1coef")
      saveRDS(result_crl_log_weight2, "result/line_plot/result_crl_log_weight2.RData")
      
      result_crl_identity <- f_data_analysis_line(data = data_ARC_crl_continuous,
                                                  link_type = "identity", # "log", "logit", "identity"
                                                  Alpha_seq = c(-1.2, -0.8, -0.4, 0, 0.4, 0.8), # for log
                                                  IF_type = "weight2",
                                                  SIM_model = "norm1coef")
      saveRDS(result_crl_identity, "result/line_plot/result_crl_identity.RData")
      
    }
    
    
  }
  
  
  # Original contour estimation - exclude IF variance, more alpha
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm")
    
    # trt_group
    {
      
      data_ARC_trt_continuous <- readRDS("ARC_trt_continuous")
      
      result_trt_log_weight1 <- f_data_analysis_line(data = data_ARC_trt_continuous,
                                                     Alpha_seq = seq(-1.2, 0.8, by = 0.02),
                                                     spline_seq_orig = c(90, 180, 270, 360),
                                                     link_type = "log", # "log", "logit", "identity"
                                                     IF_type = "weight1",
                                                     SIM_model = "norm1coef")
      saveRDS(result_trt_log_weight1, "result/contour_plot/result_trt_log_weight1.RData")
      
      
      result_trt_log_weight2 <- f_data_analysis_line(data = data_ARC_trt_continuous,
                                                     Alpha_seq = seq(-1.2, 0.8, by = 0.02),
                                                     spline_seq_orig = c(90, 180, 270, 360),
                                                     link_type = "log", # "log", "logit", "identity"
                                                     IF_type = "weight2",
                                                     SIM_model = "norm1coef")
      saveRDS(result_trt_log_weight2, "result/contour_plot/result_trt_log_weight2.RData")
      
      
      result_trt_identity <- f_data_analysis_line(data = data_ARC_trt_continuous,
                                                  Alpha_seq = seq(-1.2, 0.8, by = 0.02),
                                                  spline_seq_orig = c(90, 180, 270, 360),
                                                  link_type = "identity", # "log", "logit", "identity"
                                                  IF_type = "weight2",
                                                  SIM_model = "norm1coef")
      saveRDS(result_trt_identity, "result/contour_plot/result_trt_identity.RData")
      
    }
    
    # crl_group
    {
      
      data_ARC_crl_continuous <- readRDS("ARC_crl_continuous")
      
      result_crl_log_weight1 <- f_data_analysis_line(data = data_ARC_crl_continuous,
                                                     Alpha_seq = seq(-1.2, 0.8, by = 0.02),
                                                     spline_seq_orig = c(90, 180, 270, 360),
                                                     link_type = "log", # "log", "logit", "identity"
                                                     IF_type = "weight1",
                                                     SIM_model = "norm1coef")
      saveRDS(result_crl_log_weight1, "result/contour_plot/result_crl_log_weight1.RData")
      
      result_crl_log_weight2 <- f_data_analysis_line(data = data_ARC_crl_continuous,
                                                     Alpha_seq = seq(-1.2, 0.8, by = 0.02),
                                                     spline_seq_orig = c(90, 180, 270, 360),
                                                     link_type = "log", # "log", "logit", "identity"
                                                     IF_type = "weight2",
                                                     SIM_model = "norm1coef")
      saveRDS(result_crl_log_weight2, "result/contour_plot/result_crl_log_weight2.RData")
      
      result_crl_identity <- f_data_analysis_line(data = data_ARC_crl_continuous,
                                                  Alpha_seq = seq(-1.2, 0.8, by = 0.02),
                                                  spline_seq_orig = c(90, 180, 270, 360),
                                                  link_type = "identity", # "log", "logit", "identity"
                                                  IF_type = "weight2",
                                                  SIM_model = "norm1coef")
      saveRDS(result_crl_identity, "result/contour_plot/result_crl_identity.RData")
      
    }
    
    
  } 
  
  
}



# logit
{
  
  # Original point estimation
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm")
    
    # trt_group
    {
      
      
      data_ARC_trt_binary <- readRDS("ARC_trt_binary")
      
      
      result_trt_logit_weight1 <- f_data_analysis_point(data = data_ARC_trt_binary,
                                                        Alpha_seq = c(-0.8, -0.4, 0, 0.4, 0.8, 1.2), # for binary
                                                        link_type = "logit", # "log", "logit", "identity"
                                                        IF_type = "weight1",
                                                        SIM_model = "norm1coef")
      saveRDS(result_trt_logit_weight1, "result/CI_plot/result_trt_logit_weight1.RData")
      
      result_trt_logit_weight2 <- f_data_analysis_point(data = data_ARC_trt_binary,
                                                        Alpha_seq = c(-0.8, -0.4, 0, 0.4, 0.8, 1.2), # for binary
                                                        link_type = "logit", # "log", "logit", "identity"
                                                        IF_type = "weight2",
                                                        SIM_model = "norm1coef")
      saveRDS(result_trt_logit_weight2, "result/CI_plot/result_trt_logit_weight2.RData")
      
      
    }
    
    
    # crl_group
    {
      
      
      data_ARC_crl_binary <- readRDS("ARC_crl_binary")
      
      
      result_crl_logit_weight1 <- f_data_analysis_point(data = data_ARC_crl_binary,
                                                        Alpha_seq = c(-0.8, -0.4, 0, 0.4, 0.8, 1.2), # for binary
                                                        link_type = "logit", # "log", "logit", "identity"
                                                        IF_type = "weight1",
                                                        SIM_model = "norm1coef")
      saveRDS(result_crl_logit_weight1, "result/CI_plot/result_crl_logit_weight1.RData")
      
      result_crl_logit_weight2 <- f_data_analysis_point(data = data_ARC_crl_binary,
                                                        Alpha_seq = c(-0.8, -0.4, 0, 0.4, 0.8, 1.2), # for binary
                                                        link_type = "logit", # "log", "logit", "identity"
                                                        IF_type = "weight2",
                                                        SIM_model = "norm1coef")
      saveRDS(result_crl_logit_weight2, "result/CI_plot/result_crl_logit_weight2.RData")
      
      
    }
    
    
  }
  
  
  # Original line estimation
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm")
    
    # trt_group
    {
      
      data_ARC_trt_binary <- readRDS("ARC_trt_binary")
      
      result_trt_logit_weight1 <- f_data_analysis_line(data = data_ARC_trt_binary,
                                                       Alpha_seq = c(-0.8, -0.4, 0, 0.4, 0.8, 1.2), # for binary
                                                       link_type = "logit", # "log", "logit", "identity"
                                                       IF_type = "weight1",
                                                       SIM_model = "norm1coef")
      saveRDS(result_trt_logit_weight1, "result/line_plot/result_trt_logit_weight1.RData")
      
      result_trt_logit_weight2 <- f_data_analysis_line(data = data_ARC_trt_binary,
                                                       Alpha_seq = c(-0.8, -0.4, 0, 0.4, 0.8, 1.2), # for binary
                                                       link_type = "logit", # "log", "logit", "identity"
                                                       IF_type = "weight2",
                                                       SIM_model = "norm1coef")
      saveRDS(result_trt_logit_weight2, "result/line_plot/result_trt_logit_weight2.RData")
      
    }
    
    # crl_group
    {
      
      data_ARC_crl_binary <- readRDS("ARC_crl_binary")
      
      result_crl_log_weight1 <- f_data_analysis_line(data = data_ARC_crl_binary,
                                                     link_type = "logit", # "log", "logit", "identity"
                                                     Alpha_seq = c(-0.8, -0.4, 0, 0.4, 0.8, 1.2), # for binary
                                                     IF_type = "weight1",
                                                     SIM_model = "norm1coef")
      saveRDS(result_crl_log_weight1, "result/line_plot/result_crl_logit_weight1.RData")
      
      result_crl_log_weight2 <- f_data_analysis_line(data = data_ARC_crl_binary,
                                                     link_type = "logit", # "log", "logit", "identity"
                                                     Alpha_seq = c(-0.8, -0.4, 0, 0.4, 0.8, 1.2), # for binary
                                                     IF_type = "weight2",
                                                     SIM_model = "norm1coef")
      saveRDS(result_crl_log_weight2, "result/line_plot/result_crl_logit_weight2.RData")
      
      
    }
    
    
  }
  
  
  # Original contour estimation
  {
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm")
    
    # trt_group
    {
      
      data_ARC_trt_binary <- readRDS("ARC_trt_binary")
      
      result_trt_log_weight1 <- f_data_analysis_line(data = data_ARC_trt_binary,
                                                     Alpha_seq = seq(-0.8, 1.2, by = 0.02), # for binary
                                                     spline_seq_orig = c(90, 180, 270, 360),
                                                     link_type = "logit", # "log", "logit", "identity"
                                                     IF_type = "weight1",
                                                     SIM_model = "norm1coef")
      saveRDS(result_trt_log_weight1, "result/contour_plot/result_trt_logit_weight1.RData")
      
      
      result_trt_log_weight2 <- f_data_analysis_line(data = data_ARC_trt_binary,
                                                     Alpha_seq = seq(-0.8, 1.2, by = 0.02), # for binary
                                                     spline_seq_orig = c(90, 180, 270, 360),
                                                     link_type = "logit", # "log", "logit", "identity"
                                                     IF_type = "weight2",
                                                     SIM_model = "norm1coef")
      saveRDS(result_trt_log_weight2, "result/contour_plot/result_trt_logit_weight2.RData")
      
    }
    
    # crl_group
    {
      
      data_ARC_crl_binary <- readRDS("ARC_crl_binary")
      
      result_crl_log_weight1 <- f_data_analysis_line(data = data_ARC_crl_binary,
                                                     Alpha_seq = seq(-0.8, 1.2, by = 0.02), # for binary
                                                     spline_seq_orig = c(90, 180, 270, 360),
                                                     link_type = "logit", # "log", "logit", "identity"
                                                     IF_type = "weight1",
                                                     SIM_model = "norm1coef")
      saveRDS(result_crl_log_weight1, "result/contour_plot/result_crl_logit_weight1.RData")
      
      result_crl_log_weight2 <- f_data_analysis_line(data = data_ARC_crl_binary,
                                                     Alpha_seq = seq(-0.8, 1.2, by = 0.02), # for binary
                                                     spline_seq_orig = c(90, 180, 270, 360),
                                                     link_type = "logit", # "log", "logit", "identity"
                                                     IF_type = "weight2",
                                                     SIM_model = "norm1coef")
      saveRDS(result_crl_log_weight2, "result/contour_plot/result_crl_logit_weight2.RData")
      
    }
    
    
  } 
}




