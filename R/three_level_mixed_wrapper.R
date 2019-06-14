
 three_level_mixed_wrapper <- function(formula,
                                       data,
                                       c_id,
                                       s_id,
                                       model,
                                       timeVar = NULL,
                                       kappa = NULL,
                                       nu_v = NULL,
                                       ...){

   # formula = a two sided formula for model formulation
   #           ordinal response is expected to be numeric
   # data = data frame
   # c_id = numeric id column for clusters
   # s_id = numeric id column for individuals
   # model = bridge, normal, t, normal_t, t_normal, normal_t_no_sigma_v, normal_t_fixed_nu, two_bridge, fixed
   # ... = to be passed into the stan function from rstan

   ## checks
   if(model == "normal_t_fixed_nu" & is.null(nu_v)){
     stop("Provide nu_v")
   }

   if(model %in% c("bridge_gauss_copula_3lev",
                   "bridge_gauss_copula_2lev",
                   "bridge_t_copula_3lev")){
      if(is.null(timeVar)) stop("Provide timeVar")
      if(is.null(kappa)) stop("Provide kappa (between 0 and 2)")
   }

   ### data to be passed to stan
   y <- as.numeric(model.frame(formula, data = data)[, 1])
   x <- model.matrix(formula, data)[, -1, drop = FALSE]

   if(model %in% c("bridge", "bridge_gauss_copula_3lev", "bridge_t_copula_3lev",
                   "normal", "t", "normal_t", "t_normal",
                   "normal_t_no_sigma_v", "normal_t_fixed_nu")){

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

     if(model %in% c("bridge_gauss_copula_3lev", "bridge_t_copula_3lev")){
     dat$time  <- as.array(data[, timeVar])
     dat$kappa <- kappa
     }

  if(model == "normal_t_fixed_nu"){
    dat$nu_v <- nu_v
  }

   }else if(model %in% c("two_bridge", "bridge_gauss_copula_2lev")){

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

     if(model == "bridge_gauss_copula_2lev"){
        dat$time <- data[, timeVar]
        dat$kappa <- kappa
     }

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
     mod <- rstan::stan_model(model_code = bridge_ordinal_mixed_threelev_u_v_prior_on_var,
                              auto_write = TRUE)
     res <- rstan::sampling(mod, data = dat, ...)
   }

   if(model == "bridge_gauss_copula_3lev"){
      mod <- rstan::stan_model(model_code = bridge_ordinal_mixed_threelev_gauss_copula,
                               auto_write = TRUE)
      res <- rstan::sampling(mod, data = dat, ...)
   }

   if(model == "bridge_gauss_copula_2lev"){
      mod <- rstan::stan_model(model_code = bridge_ordinal_mixed_twolev_gauss_copula,
                               auto_write = TRUE)
      res <- rstan::sampling(mod, data = dat, ...)
   }

   if(model == "bridge_t_copula_3lev"){
      mod <- rstan::stan_model(model_code = bridge_ordinal_mixed_threelev_t_copula,
                               auto_write = TRUE)
      res <- rstan::sampling(mod, data = dat, ...)
   }

   ### normal distribution for both u and v
   if(model == "normal"){
     mod <- rstan::stan_model(model_code = normal_ordinal_mixed_threelev_reparam,
                              auto_write = TRUE)
     res <- rstan::sampling(mod, data = dat, ...)
   }

   ### t distribution for both u and v
   if(model == "t"){
     mod <- rstan::stan_model(model_code = t_ordinal_mixed_threelev_reparam_prior_on_var_inv_nu,
                              auto_write = TRUE)
     res <- rstan::sampling(mod, data = dat, ...)
   }

   ### normal for U, t for V
   if(model == "normal_t"){
     mod <- rstan::stan_model(model_code = normal_t_ordinal_mixed_threelev_reparam_prior_on_var_inv_nu,
                              auto_write = TRUE)
     res <- rstan::sampling(mod, data = dat, ...)
   }

   ### t for U, normal for V
   if(model == "t_normal"){
     mod <- rstan::stan_model(model_code = t_normal_ordinal_mixed_threelev_reparam_prior_on_var_inv_nu,
                              auto_write = TRUE)
     res <- rstan::sampling(mod, data = dat, ...)
   }

   ### normal for U, t for V - no sigma for V
   if(model == "normal_t_no_sigma_v"){
     mod <- rstan::stan_model(model_code = normal_t_ordinal_mixed_threelev_reparam_prior_on_var_inv_nu_no_sigma_v,
                              auto_write = TRUE)
     res <- rstan::sampling(mod, data = dat, ...)
   }

   ### normal for U, t for V - fixed nu
   if(model == "normal_t_fixed_nu"){
     mod <- rstan::stan_model(model_code = normal_t_ordinal_mixed_threelev_reparam_prior_on_var_fixed_nu,
                              auto_write = TRUE)
     res <- rstan::sampling(mod, data = dat, ...)
   }

   ## two level model: no u, bridge distributed v

   if(model == "two_bridge"){
     mod <- rstan::stan_model(model_code = bridge_ordinal_mixed_twolev_v_prior_on_var,
                              auto_write = TRUE)
     res <- rstan::sampling(mod, data = dat, ...)
   }

   ## fixed effects model
   if(model == "fixed"){
     mod <- rstan::stan_model(model_code = ordinal_fixed,
                              auto_write = TRUE)
     res <- rstan::sampling(mod, data = dat, ...)
   }

   return(res)

 }

