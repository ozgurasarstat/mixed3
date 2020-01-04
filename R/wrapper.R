
 wrapper <- function(formula,
                     data,
                     c_id,
                     s_id,
                     model,
                     timeVar,
                     ...
                     ){

  # formula = a two sided formula for model formulation ordinal response is expected to be numeric
  # data = a data frame
  # c_id = numeric id column for clusters
  # s_id = numeric id column for individuals
  # model = bridge_threelev, bridge_twolev, normal_threelev, fixed
  # ... = to be passed into the stan function from rstan

  ## sort the data wrt c_id, s_id, timeVar
  data <- data[order(data[, c_id], data[, s_id], data[, timeVar]), ]

  ### data to be passed to stan
  y <- as.numeric(model.frame(formula, data = data)[, 1])
  x <- model.matrix(formula, data)[, -1, drop = FALSE]
  ntot <- nrow(data)
  p <- ncol(x)
  k <- nlevels(factor(y))

  if(model != "fixed"){

    cluster_id_orig <- data[, c_id]
    subj_id_orig    <- data[, s_id]

    cluster_id_unique <- unique(cluster_id_orig)
    subj_id_unique    <- unique(subj_id_orig)

    ncluster <- length(cluster_id_unique)
    nsubj    <- length(subj_id_unique)

    cluster_id_seq <- 1:ncluster
    subj_id_seq    <- 1:nsubj

    if(model %in% c("bridge_threelev", "normal_threelev")){
      cluster_id <- subj_id <- rep(0, ntot)
      for(i in 1:ntot){
        cluster_id[i] <- cluster_id_seq[which(cluster_id_orig[i] == cluster_id_unique)]
        subj_id[i]    <- subj_id_seq[which(subj_id_orig[i] == subj_id_unique)]
      }

      dat <- list(ntot = ntot,
                  cluster_id = cluster_id,
                  subj_id = subj_id,
                  y = y,
                  p = p,
                  x = x,
                  ncluster = ncluster,
                  nsubj = nsubj,
                  k = k)

    }else if(model == "bridge_twolev"){
      subj_id <- rep(0, ntot)
      for(i in 1:ntot){
        subj_id[i]    <- subj_id_seq[which(subj_id_orig[i] == subj_id_unique)]
      }

      dat <- list(ntot = ntot,
                  subj_id = subj_id,
                  y = y,
                  p = p,
                  x = x,
                  nsubj = nsubj,
                  k = k)

    }

  }else if(model == "fixed"){

    dat <- list(y = y,
                x = x,
                ntot = ntot,
                p = p,
                k = k)

  }

  ## fits

  if(model == "bridge_threelev"){
    mod <- rstan::stan_model(model_code = bridge_threelev, auto_write = TRUE)
    res <- rstan::sampling(mod, data = dat, ...)
  }

  if(model == "normal_threelev"){
    mod <- rstan::stan_model(model_code = normal_threelev, auto_write = TRUE)
    res <- rstan::sampling(mod, data = dat, ...)
  }

  if(model == "bridge_twolev"){
    mod <- rstan::stan_model(model_code = bridge_twolev, auto_write = TRUE)
    res <- rstan::sampling(mod, data = dat, ...)
  }

  if(model == "fixed"){
    mod <- rstan::stan_model(model_code = ordinal_fixed, auto_write = TRUE)
    res <- rstan::sampling(mod, data = dat, ...)
  }

  return(res)

}
