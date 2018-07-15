 
 three_level_mixed_wrapper <- function(formula, 
                                       data, 
                                       c_id, 
                                       s_id, 
                                       model, 
                                       wd, 
                                       ...){
   
   # formula = a two sided formula for model formulation
   #           ordinal response is expected to be numeric
   # data = data frame 
   # c_id = numeric id column for clusters
   # s_id = numeric id column for individuals
   # model = bridge, normal, t, normal_t, t_normal, normal_t_no_sigma_v, two_bridge, fixed
   # wd = directory for the folder where the stan codes are located
   # ... = to be passed into the stan function from rstan
   
   ### load rstan package
   library(rstan)
   rstan_options(auto_write = TRUE)
   
   ### direct R to the folder where the codes are located in
   setwd(wd)
   
   ### data to be passed to stan   
   y <- as.numeric(model.frame(formula, data = data)[, 1])
   x <- model.matrix(formula, data)[, -1, drop = FALSE]
   
   if(model %in% c("bridge", "normal", "t", "normal_t", "t_normal", "normal_t_no_sigma_v")){
     
     cluster_id <- data[, c_id]
     subj_id <- data[, s_id]
     ntot <- nrow(data)
     p <- ncol(x)
     ncluster <- length(unique(cluster_id))
     nsubj <- length(unique(subj_id))
     k <- nlevels(factor(y))
     
     dat <- list(y = y, 
                 x = x,
                 cluster_id = cluster_id,
                 subj_id = subj_id,
                 ntot = ntot,
                 p = p,
                 ncluster = ncluster,
                 nsubj = nsubj,
                 k = k) 
     
   }else if(model %in% "two_bridge"){
     
     subj_id <- data[, s_id]
     ntot <- nrow(data)
     p <- ncol(x)
     nsubj <- length(unique(subj_id))
     k <- nlevels(factor(y))
     
     dat <- list(y = y, 
                 x = x,
                 subj_id = subj_id,
                 ntot = ntot,
                 p = p,
                 nsubj = nsubj,
                 k = k) 
     
   }else if(model %in% "fixed"){
     
     ntot <- nrow(data)
     p <- ncol(x)
     k <- nlevels(factor(y))
     
     dat <- list(y = y, 
                 x = x,
                 ntot = ntot,
                 p = p,
                 k = k) 
     
   }
   

   ### modified bridge distribution for u and bridge for v
   if(model == "bridge"){
 
       source("bridge_ordinal_mixed_threelev_u_v_prior_on_var.R")
       res <- stan(model_code = bridge_ordinal_mixed_threelev_u_v_prior_on_var, 
                   data = dat,
                   ...)

   }
   
   ### normal distribution for both u and v
   if(model == "normal"){
       
         source("normal_ordinal_mixed_threelev_reparam.R")
         res <- stan(model_code = normal_ordinal_mixed_threelev_reparam, 
                     data = dat,
                     ...)
  
   }
   
   ### t distribution for both u and v
   if(model == "t"){
       
       source("t_ordinal_mixed_threelev_reparam_prior_on_var_inv_nu.R")
       res <- stan(model_code = t_ordinal_mixed_threelev_reparam_prior_on_var_inv_nu, 
                   data = dat,
                   ...)

   }
   
   ### normal for U, t for V
   if(model == "normal_t"){
     source("normal_t_ordinal_mixed_threelev_reparam_prior_on_var_inv_nu.R")
     res <- stan(model_code = normal_t_ordinal_mixed_threelev_reparam_prior_on_var_inv_nu, 
                 data = dat,
                 ...)
   }
   
   ### t for U, normal for V
   if(model == "t_normal"){
     source("t_normal_ordinal_mixed_threelev_reparam_prior_on_var_inv_nu.R")
     res <- stan(model_code = t_normal_ordinal_mixed_threelev_reparam_prior_on_var_inv_nu, 
                 data = dat,
                 ...)
   }
   
   ### normal for U, t for V - no sigma for V
   if(model == "normal_t_no_sigma_v"){
     source("normal_t_ordinal_mixed_threelev_reparam_prior_on_var_inv_nu_no_sigma_v.R")
     res <- stan(model_code = normal_t_ordinal_mixed_threelev_reparam_prior_on_var_inv_nu_no_sigma_v, 
                 data = dat,
                 ...)
   }
   
   ## two level model: no u, bridge distributed v
   
   if(model == "two_bridge"){
     source("bridge_ordinal_mixed_twolev_v_prior_on_var.R")
     res <- stan(model_code = bridge_ordinal_mixed_twolev_v_prior_on_var, 
                 data = dat,
                 ...)
   }
   
   if(model == "fixed"){
     source("ordinal_fixed.R")
     res <- stan(model_code = ordinal_fixed, 
                 data = dat,
                 ...)
   }
   
   return(res)
 }

 