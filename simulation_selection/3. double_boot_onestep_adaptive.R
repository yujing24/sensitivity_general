rm(list=ls())

source("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/selection_code/double_boot_functions.R")

Boot_outer <- 400
Boot_inner <- 50



# log link
{
  
  
  sample_size = 200
  group = "control" # "control"
  
  # Count outcome - negative binomial model
  setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/parametric_sim_data")
  dataall <- readRDS(paste0("SIM_sim_data_", group, "_orig_",sample_size, "_neg"))
  
  link_type = "log"
  IF_type = "weight2"
  
  if(group == "treat"){
    grp = "trt"
  }else{
    grp = "crl"
  }
  if(link_type == "log"){
    range = c(1.2, 3)
  }else{
    range = c(0.2, 0.5) 
  }
  
  
  for(i in c()){
    
    print(i)
    
    # Generate original information
    # start <- Sys.time()
    orig_input <- fun_est_IF_onestep_double(index = i, 
                                            link_type = link_type, # "log", "logit", "identity"
                                            IF_type =  IF_type,
                                            range = range)
    # end <- Sys.time()
    # end - start
    # 1.859556 mins
    
    # Do outer bootstrap
    if(group == "treat"){
      set.seed(i) # treatment group
    }else{
      set.seed(i*123) # control group
    }
    
    Boot_seed_all <- sample(500000000, size = Boot_outer, replace = T)
    
    # start <- Sys.time()
    result_double <- parallel::mclapply(Boot_seed_all, 
                                        f_double_boot_onestep_adaptive,
                                        orig_input = orig_input,
                                        link_type = link_type,
                                        IF_type = IF_type,
                                        range = range,
                                        search_method = "adaptive", # "adaptive" or "grid"
                                        L = Boot_inner, # inner bootstrap size
                                        mc.cores = 96)
    
    # end <- Sys.time()
    # end - start
    # 5.250098 mins
    
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new/double_boot_400_50")
    saveRDS(result_double, file = paste0(link_type, "_", sample_size,
                                       "/", 
                                       group, "_", IF_type, 
                                       "/", 
                                       "adaptive_", grp, "_", sample_size, "_", link_type, "_", IF_type,"_", i, "_double"))
    
  }
  
  
}



# logit link - weight1
{
  
  
  sample_size = 200
  group = "control" # "control"
  
  # Count outcome - negative binomial model
  setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/parametric_sim_data")
  dataall <- readRDS(paste0("SIM_sim_data_", group, "_binary_", sample_size, "_logistic"))
  
  link_type = "logit"
  IF_type = "weight1"
  
  if(group == "treat"){
    grp = "trt"
  }else{
    grp = "crl"
  }
  if(link_type == "log"){
    range = c(1.2, 3)
  }else{
    range = c(0.2, 0.5) 
  }
  
  
  for(i in c(328:330)){
    
    print(i)
    
    # Generate original information
    start <- Sys.time()
    orig_input <- fun_est_IF_onestep_double(index = i, 
                                            link_type = link_type, # "log", "logit", "identity"
                                            IF_type =  IF_type,
                                            range = range)
    end <- Sys.time()
    end - start
    # 2.256671 mins
    
    # Do outer bootstrap
    if(group == "treat"){
      set.seed(i) # treatment group
    }else{
      set.seed(i*123) # control group
    }
    
    Boot_seed_all <- sample(500000000, size = Boot_outer, replace = T)
    
    start <- Sys.time()
    result_double <- parallel::mclapply(Boot_seed_all, 
                                        f_double_boot_onestep_adaptive,
                                        orig_input = orig_input,
                                        link_type = link_type,
                                        IF_type = IF_type,
                                        range = range,
                                        search_method = "adaptive", # "adaptive" or "grid"
                                        L = Boot_inner, # inner bootstrap size
                                        mc.cores = 96)
    
    end <- Sys.time()
    end - start
    # 
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new/double_boot_400_50")
    saveRDS(result_double, file = paste0(link_type, "_", sample_size,
                                         "/", 
                                         group, "_", IF_type, 
                                         "/", 
                                         "double_", grp, "_", sample_size, "_", link_type, "_", IF_type,"_", i, "_corr"))
    
  }
  
  
}



# logit link - weight2
{
  
  
  sample_size = 200
  group = "treat" # "control"
  
  # Count outcome - negative binomial model
  setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/parametric_sim_data")
  dataall <- readRDS(paste0("SIM_sim_data_", group, "_binary_", sample_size, "_logistic"))
  
  link_type = "logit"
  IF_type = "weight2"
  
  if(group == "treat"){
    grp = "trt"
  }else{
    grp = "crl"
  }
  if(link_type == "log"){
    range = c(1.2, 3)
  }else{
    range = c(0.2, 0.5) 
  }
  
  
  for(i in c(310)){
    
    print(i)
    
    # Generate original information
    start <- Sys.time()
    orig_input <- fun_est_IF_onestep_double(index = i, 
                                            link_type = link_type, # "log", "logit", "identity"
                                            IF_type =  IF_type,
                                            range = range)
    end <- Sys.time()
    end - start
    # 2.256671 mins
    
    # Do outer bootstrap
    if(group == "treat"){
      set.seed(i) # treatment group
    }else{
      set.seed(i*123) # control group
    }
    
    Boot_seed_all <- sample(500000000, size = Boot_outer, replace = T)
    
    start <- Sys.time()
    result_double <- parallel::mclapply(Boot_seed_all, 
                                        f_double_boot_onestep_adaptive,
                                        orig_input = orig_input,
                                        link_type = link_type,
                                        IF_type = IF_type,
                                        range = range,
                                        search_method = "adaptive", # "adaptive" or "grid"
                                        L = Boot_inner, # inner bootstrap size
                                        mc.cores = 96)
    
    end <- Sys.time()
    end - start
    # 
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new/double_boot_400_50")
    saveRDS(result_double, file = paste0(link_type, "_", sample_size,
                                         "/", 
                                         group, "_", IF_type, 
                                         "/", 
                                         "double_", grp, "_", sample_size, "_", link_type, "_", IF_type,"_", i, "_corr"))
    
  }
  
  
}



# identity link
{
  
  
  sample_size = 200
  group = "control" # "control"
  
  # Count outcome - negative binomial model
  setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/parametric_sim_data")
  dataall <- readRDS(paste0("SIM_sim_data_", group, "_orig_",sample_size, "_neg"))
  
  link_type = "identity"
  
  if(group == "treat"){
    grp = "trt"
  }else{
    grp = "crl"
  }
  if(link_type == "identity"){
    range = c(1.2, 3)
  }else{
    range = c(0.2, 0.5) 
  }
  
  
  for(i in c(1:500)){
    
    print(i)
    
    # Generate original information
    # start <- Sys.time()
    orig_input <- fun_est_IF_onestep_double(index = i, 
                                            link_type = link_type, # "log", "logit", "identity"
                                            IF_type =  IF_type,
                                            range = range)
    # end <- Sys.time()
    # end - start
    # 1.859556 mins
    
    # Do outer bootstrap
    if(group == "treat"){
      set.seed(i) # treatment group
    }else{
      set.seed(i*123) # control group
    }
    
    Boot_seed_all <- sample(500000000, size = Boot_outer, replace = T)
    
    start <- Sys.time()
    result_double <- parallel::mclapply(Boot_seed_all, 
                                        f_double_boot_identity_adaptive,
                                        orig_input = orig_input,
                                        link_type = link_type,
                                        range = range,
                                        search_method = "adaptive", # "adaptive" or "grid"
                                        L = Boot_inner, # inner bootstrap size
                                        mc.cores = 96)
    
    end <- Sys.time()
    end - start
    # 5.250098 mins
    
    
    setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result/select_alpha_new/double_boot_400_50")
    saveRDS(result_double, file = paste0(link_type, "_", sample_size,
                                         "/", 
                                         group,
                                         "/", 
                                         "adaptive_", grp, "_", sample_size, "_", link_type, "_", i, "_double"))
    
    
  }
  
  
}

