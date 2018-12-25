
 three_level_mixed_wrapper <- function(formula,
                                       data,
                                       c_id,
                                       s_id,
                                       model,
                                       ...){

   # formula = a two sided formula for model formulation
   #           ordinal response is expected to be numeric
   # data = data frame
   # c_id = numeric id column for clusters
   # s_id = numeric id column for individuals
   # model = bridge, normal, t, normal_t, t_normal, normal_t_no_sigma_v, two_bridge, fixed
   # ... = to be passed into the stan function from rstan

   ### data to be passed to stan
   y <- as.numeric(model.frame(formula, data = data)[, 1])
   x <- model.matrix(formula, data)[, -1, drop = FALSE]

   if(model %in% c("bridge", "normal", "t", "normal_t", "t_normal", "normal_t_no_sigma_v")){

     nrepeat_c <- data[, c_id] %>% table %>% as.numeric
     nrepeat_s <- data[, s_id] %>% table %>% as.numeric

     ncluster <- data[, c_id] %>% unique %>% length
     nsubj    <- data[, s_id] %>% unique %>% length

     cluster_id <- rep(1:ncluster, nrepeat_c)#data[, c_id]
     subj_id <- rep(1:nsubj, nrepeat_s)#data[, s_id]

     ntot <- nrow(data)
     p <- ncol(x)
     #ncluster <- length(unique(cluster_id))
     #nsubj <- length(unique(subj_id))
     k <- nlevels(factor(y))

     cumsum_nrepeat_c <- cumsum(nrepeat_c)
     cumsum_nrepeat_s <- cumsum(nrepeat_s)

     ind_c <- cbind(c(1, (cumsum_nrepeat_c[-ncluster] + 1)), cumsum_nrepeat_c)
     ind_s <- cbind(c(1, (cumsum_nrepeat_s[-nsubj] + 1)), cumsum_nrepeat_s)

     dat <- list(y = y,
                 x = x,
                 cluster_id = cluster_id,
                 subj_id = subj_id,
                 ntot = ntot,
                 p = p,
                 ncluster = ncluster,
                 nsubj = nsubj,
                 k = k,
                 ind_c = ind_c,
                 ind_s = ind_s,
                 nrepeat_c = nrepeat_c,
                 nrepeat_s = nrepeat_s)

   }else if(model %in% "two_bridge"){

     nrepeat_s <- data[, s_id] %>% table %>% as.numeric
     nsubj     <- data[, s_id] %>% unique %>% length

     subj_id <- rep(1:nsubj, nrepeat_s)#data[, s_id]

     ntot <- nrow(data)
     p <- ncol(x)

     k <- nlevels(factor(y))

     cumsum_nrepeat_s <- cumsum(nrepeat_s)
     ind_s <- cbind(c(1, (cumsum_nrepeat_s[-nsubj] + 1)), cumsum_nrepeat_s)

     #subj_id <- data[, s_id]
     #ntot <- nrow(data)
     #p <- ncol(x)
     #nsubj <- length(unique(subj_id))
     #k <- nlevels(factor(y))

     dat <- list(y = y,
                 x = x,
                 subj_id = subj_id,
                 ntot = ntot,
                 p = p,
                 nsubj = nsubj,
                 k = k,
                 ind_s = ind_s,
                 nrepeat_s = nrepeat_s)

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
       res <- stan(model_code = bridge_ordinal_mixed_threelev_u_v_prior_on_var,
                   data = dat,
                   ...)
   }

   ### normal distribution for both u and v
   if(model == "normal"){
         res <- stan(model_code = normal_ordinal_mixed_threelev_reparam,
                     data = dat,
                     ...)
   }

   ### t distribution for both u and v
   if(model == "t"){
       res <- stan(model_code = t_ordinal_mixed_threelev_reparam_prior_on_var_inv_nu,
                   data = dat,
                   ...)
   }

   ### normal for U, t for V
   if(model == "normal_t"){
     res <- stan(model_code = normal_t_ordinal_mixed_threelev_reparam_prior_on_var_inv_nu,
                 data = dat,
                 ...)
   }

   ### t for U, normal for V
   if(model == "t_normal"){
     res <- stan(model_code = t_normal_ordinal_mixed_threelev_reparam_prior_on_var_inv_nu,
                 data = dat,
                 ...)
   }

   ### normal for U, t for V - no sigma for V
   if(model == "normal_t_no_sigma_v"){
     res <- stan(model_code = normal_t_ordinal_mixed_threelev_reparam_prior_on_var_inv_nu_no_sigma_v,
                 data = dat,
                 ...)
   }

   ## two level model: no u, bridge distributed v

   if(model == "two_bridge"){
     res <- stan(model_code = bridge_ordinal_mixed_twolev_v_prior_on_var,
                 data = dat,
                 ...)
   }

   if(model == "fixed"){
     res <- stan(model_code = ordinal_fixed,
                 data = dat,
                 ...)
   }
   return(res)
 }

