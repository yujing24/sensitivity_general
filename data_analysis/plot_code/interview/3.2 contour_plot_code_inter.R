rm(list=ls())
library(ggplot2)
library(metR)
library(isoband)
library(dplyr)

setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/data_analysis")

# treat: -0.6, 0.54
# control: -0.6, 0.26

##############   Month 6  ############## 
{
  
  df_m6 <- readRDS("result/contour_plot_inter/Contour_plot_month6_SIM.rds")

  # check the estimated value and confidence interval under AAR assumption 
  df_m6_restrict <- filter(df_m6, 
                           Alpha_trt <= 0.54,
                           Alpha_trt >= -0.6,
                           Alpha_crl <= 0.26,
                           Alpha_crl >= -0.6)
  max(df_m6_restrict$Estimate_vec)
  min(df_m6_restrict$Estimate_vec)
  
  
  SIM_p1 <- ggplot(df_m6_restrict, aes(x = Alpha_crl, y = Alpha_trt, z = Estimate_vec))+
    geom_contour_fill(binwidth=0.4, show.legend = F) +
    scale_fill_divergent(mid="white",low="forestgreen",high="orange") +
    geom_contour2(aes(z=Estimate_vec,label=..level..),
                  colour="black", size=0.5, binwidth=0.2, skip=0.2,
                  margin=grid::unit(c(5,5,5,5),"pt"),
                  label.placer = label_placer_fraction(frac=0.38,
                                                       rot_adjuster=angle_fixed()),
                  label_size=7) + # in package we will remove this geom_contour2() function
    # geom_contour2(data = df2_m6_restrict, aes(x=Alpha_crl, y=Alpha_trt, z=CI_info, label=..level..)) +
    theme_bw() +
    scale_colour_divergent(mid="yellow",low="forestgreen",high="orange")+
    labs(x="Alpha for PT Arm",y="Alpha for PT+HV Arm")+
    theme(axis.title = element_text(size=25),
          axis.text = element_text(size=18))+
    geom_vline(xintercept=0,lty="dotted",size=1)+
    geom_hline(yintercept=0,lty="dotted",size=1)+
    geom_point(x=0, y=0, size=3)+
    scale_x_continuous(breaks=c(-0.6, -0.3, 0, 0.3, 0.6))+
    scale_y_continuous(breaks=c(-0.6, -0.3, 0, 0.3, 0.6))
  
  SIM_p1
  
  ggsave(plot = SIM_p1, 
         width = 7.5, 
         height = 6, 
         dpi = 500, 
         filename = "result/plots_inter/contour_plot/contour_month6_mean.pdf")

  
  df2_m6 <- mutate(df_m6_restrict , CI_info = Lower_ci*(Lower_ci>0) + Upper_ci*(Upper_ci<0))
  
  min(df2_m6$CI_info)
  max(df2_m6$CI_info)

  SIM_p2 <- ggplot(df2_m6, aes(x = Alpha_crl, y=Alpha_trt, z=CI_info))+
    geom_contour_fill(breaks = c(seq(-1.1, 1.1, by = 0.2))) +
    scale_fill_divergent(mid="white",low="forestgreen",high="orange",
                         midpoint = 0) +
    theme_bw()+
    theme(legend.title = element_blank())+
    labs(x="Alpha for PT Arm",y="Alpha for PT+HV Arm")+
    theme(axis.title = element_text(size=25),
          axis.text = element_text(size=18),
          legend.text=element_text(size=14),
          legend.key.size = unit(1,"cm"))+
    geom_vline(xintercept=0,lty="dotted",size=1)+
    geom_hline(yintercept=0,lty="dotted",size=1)+
    geom_point(x=0,y=0,size=3)+
    scale_x_continuous(breaks=c(-0.6, -0.3, 0, 0.3, 0.6))+
    scale_y_continuous(breaks=c(-0.6, -0.3, 0, 0.3, 0.6))
  
  SIM_p2
  
  ggsave(plot = SIM_p2, 
         width = 8, 
         height = 6, 
         dpi = 500, 
         filename = "result/plots_inter/contour_plot/contour_month6_diff.pdf")

  
}


##############   Month 12  ############## 
{
  
  
  #################################
  df_m12 <- readRDS("result/contour_plot/Contour_plot_month12_SIM.rds")
  
  df_m12_restrict <- filter(df_m12, 
                           Alpha_trt <= 0.54,
                           Alpha_trt >= -0.6,
                           Alpha_crl <= 0.26,
                           Alpha_crl >= -0.6)
  max(df_m12_restrict$Estimate_vec)
  min(df_m12_restrict$Estimate_vec)
  
  SIM_p3 <- ggplot(df_m12_restrict, aes(x=Alpha_crl, y=Alpha_trt, z=Estimate_vec))+
    geom_contour_fill(binwidth=0.4, show.legend = F) +
    scale_fill_divergent(mid="white",low="forestgreen",high="orange") +
    geom_contour2(aes(z=Estimate_vec,label=..level..),
                  colour="black", size=0.5, binwidth=0.2, skip=0.2,
                  margin=grid::unit(c(5,5,5,5),"pt"),
                  label.placer = label_placer_fraction(frac=0.38,
                                                       rot_adjuster=angle_fixed()),
                  label_size=7) + 
    theme_bw()+
    scale_colour_divergent(mid="yellow",low="forestgreen",high="orange")+
    labs(x="Alpha for PT Arm",y="Alpha for PT+HV Arm")+
    theme(axis.title = element_text(size=25),
          axis.text = element_text(size=18))+
    geom_vline(xintercept=0,lty="dotted",size=1)+
    geom_hline(yintercept=0,lty="dotted",size=1)+
    geom_point(x=0, y=0, size=3)+
    scale_x_continuous(breaks=c(-0.6, -0.3, 0, 0.3, 0.6))+
    scale_y_continuous(breaks=c(-0.6, -0.3, 0, 0.3, 0.6))
  
  SIM_p3
  
  ggsave(plot = SIM_p3, 
         width = 7, 
         height = 6, 
         dpi = 500, 
         filename = "result/plots_inter/contour_plot/contour_month12_mean.pdf")
  
  
  df2_m12 <- mutate(df_m12_restrict, CI_info = Lower_ci*(Lower_ci>0) + Upper_ci*(Upper_ci<0))
  
  min(df2_m12$CI_info)
  max(df2_m12$CI_info)
  
  SIM_p4 <- ggplot(df2_m12,aes(x=Alpha_crl, y=Alpha_trt, z=CI_info))+
    geom_contour_fill(breaks = c(seq(-0.9, 0.9, by = 0.2))) +
    scale_fill_divergent(mid="white",low="forestgreen",high="orange",
                         midpoint = 0) +
    theme_bw()+
    theme(legend.title = element_blank())+
    labs(x="Alpha for PT Arm",y="Alpha for PT+HV Arm")+
    theme(axis.title = element_text(size=25),
          axis.text = element_text(size=18),
          legend.text=element_text(size=14),
          legend.key.size = unit(1,"cm"))+
    geom_vline(xintercept=0,lty="dotted",size=1)+
    geom_hline(yintercept=0,lty="dotted",size=1)+
    geom_point(x=0,y=0,size=3)+
    scale_x_continuous(breaks=c(-0.6, -0.3, 0, 0.3, 0.6))+
    scale_y_continuous(breaks=c(-0.6, -0.3, 0, 0.3, 0.6))
  
  SIM_p4
  
  ggsave(plot = SIM_p4, 
         width = 8, 
         height = 6, 
         dpi = 500, 
         filename = "result/plots_inter/contour_plot/contour_month12_diff.pdf")
  
  
}

