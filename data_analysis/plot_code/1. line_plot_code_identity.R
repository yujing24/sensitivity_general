rm(list=ls())
library(ggplot2) 
library(dplyr)

setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis_norm")

Alpha_seq <- c(-1.2, -0.8, -0.4, 0, 0.4, 0.8)
spline_seq <- seq(60, 460, by = 1)
range_value <- c(1.2, 3)


## identity link
{
  
  df_trt <- readRDS("result/line_plot/result_trt_identity.RData")
  df_crl <- readRDS("result/line_plot/result_crl_identity.RData")
  
  ### control group
  plot_df_crl <- data.frame(Time = rep(spline_seq, 6),
                            Mean=c(df_crl[[1]][1, -1], df_crl[[1]][2, -1], df_crl[[1]][3, -1],  
                                   df_crl[[1]][4, -1], df_crl[[1]][5, -1], df_crl[[1]][6, -1]),
                            Alpha=c(rep("-1.2", length(spline_seq)),
                                    rep("-0.8", length(spline_seq)),
                                    rep("-0.4", length(spline_seq)),
                                    rep("0", length(spline_seq)),
                                    rep("0.4", length(spline_seq)),
                                    rep("0.8", length(spline_seq))
                            ))
  
  orig_crl <- data.frame(Time = orig_simu_crl$time,
                         Mean = orig_simu_crl$Asthma_control,
                         Alpha = rep("original data", nrow(orig_simu_crl)))
  
  line_plot_crl <- ggplot(plot_df_crl, aes(x = Time, y = Mean, colour = Alpha))+
    geom_line(linewidth = 1.5) + 
    # geom_point(data = orig_crl, size = 0.5) +
    theme_bw()+
    expand_limits(x=c(30, 480),y=c(0, 5))+
    scale_x_continuous(breaks=seq(90, 450, by=90))+
    scale_colour_manual(values = c("#253494", "#41b6c4", 
                                   "#a1dab4", "#ffffb2",
                                   "#fd8d3c", "#bd0026"),
                        # "#f03b20",# "#fecc5c",# "#2c7fb8", # "grey50"),
                        breaks=c("0.8", "0.4", "0", "-0.4", "-0.8", "-1.2")) + 
    # "original data"))+
    labs(x="Time (Days)", y="ACQ Score", color = expression(alpha)) +
    labs(color=expression(alpha))+
    theme(axis.title=element_text(size=25),
          legend.title=element_text(size=20),
          axis.text=element_text(size=18),
          legend.text=element_text(size=18))+
    # theme(legend.position = "none")+
    annotate("text", x = 400, y = 4.5,label="PT Arm",size=6) +
    geom_hline(yintercept=range_value[1],lty="dotted",size=2)+
    geom_hline(yintercept=range_value[2],lty="dotted",size=2)
  
  line_plot_crl
  
  ggsave(plot = line_plot_crl, 
         width = 8, 
         height = 6, 
         dpi = 500, 
         filename = "result/plots/line_plot/line_crl_identity.pdf")
  
  
  ### treatment group
  
  plot_df_trt <- data.frame(Time = rep(spline_seq, 6),
                            Mean=c(df_trt[[1]][1, -1], df_trt[[1]][2, -1], df_trt[[1]][3, -1],  
                                   df_trt[[1]][4, -1], df_trt[[1]][5, -1], df_trt[[1]][6, -1]),
                            Alpha=c(rep("-1.2", length(spline_seq)),
                                    rep("-0.8", length(spline_seq)),
                                    rep("-0.4", length(spline_seq)),
                                    rep("0", length(spline_seq)),
                                    rep("0.4", length(spline_seq)),
                                    rep("0.8", length(spline_seq))
                            ))
  
  orig_trt <- data.frame(Time = orig_simu_trt$time,
                         Mean = orig_simu_trt$Asthma_control,
                         Alpha = rep("original data", nrow(orig_simu_trt)))
  
  line_plot_trt <- ggplot(plot_df_trt, aes(x = Time, y = Mean, colour = Alpha))+
    geom_line(linewidth = 1.5) + 
    # geom_point(data = orig_trt, size = 0.5)+
    theme_bw()+
    expand_limits(x=c(30, 480),y=c(0, 5))+
    scale_x_continuous(breaks=seq(90, 450, by=90))+
    scale_colour_manual(values = c("#253494", "#41b6c4", 
                                   "#a1dab4", "#ffffb2",
                                   "#fd8d3c", "#bd0026"),
                        # "#f03b20",# "#fecc5c",# "#2c7fb8", # "grey50"),
                        breaks=c("0.8", "0.4", "0", "-0.4", "-0.8", "-1.2")) + 
    # "original data"))+
    labs(x="Time (Days)", y="ACQ Score", color = expression(alpha)) +
    labs(color=expression(alpha))+
    theme(axis.title=element_text(size=25),
          legend.title=element_text(size=20),
          axis.text=element_text(size=18),
          legend.text=element_text(size=18))+
    theme(legend.position = "none")+
    annotate("text", x = 400, y = 4.5,label="PT+HV Arm",size=6)+
    geom_hline(yintercept=range_value[1],lty="dotted",size=2)+
    geom_hline(yintercept=range_value[2],lty="dotted",size=2)
  
  line_plot_trt
  
  
  ggsave(plot = line_plot_trt, 
         width = 7, 
         height = 6, 
         dpi = 500, 
         filename = "result/plots/line_plot/line_trt_identity.pdf")
  
  
}
