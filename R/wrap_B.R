
wrap_B <- function(formula,
                   data,
                   ec_id,
                   c_id,
                   s_id,
                   timeVar,
                   kappa = list(kappa1 = 1, kappa2 = 1),
                   model = "bridge",
                   reff = list(U = "intercept", V = "serial"),
                   ...){

  if(model != "fixed"){
    if(reff$U == "none" & reff$V == "none"){
      model <- "fixed"
    }
  }

  if(model != "fixed"){
    if(reff$U != "intercept" | reff$V != "intercept"){
      nrepeat_ec <- as.numeric(table(data[, ec_id]))
      cumsum_nrepeat_ec <- cumsum(nrepeat_ec)
      ind_ec <- cbind(c(1, (cumsum_nrepeat_ec[-length(unique(data[, ec_id]))] + 1)), cumsum_nrepeat_ec)
    }
  }

  y <- as.numeric(model.frame(formula, data = data)[, 1])
  x <- model.matrix(formula, data)[, -1, drop = FALSE]

  ntot <- nrow(data)

  cluster_id_orig <- data[, c_id]
  subj_id_orig <- data[, s_id]

  cluster_id_unique <- unique(cluster_id_orig)
  subj_id_unique <- unique(subj_id_orig)

  ncluster <- length(cluster_id_unique)
  nsubj <- length(subj_id_unique)

  cluster_id_seq <- 1:ncluster
  subj_id_seq <- 1:nsubj

  cluster_id <- subj_id <- rep(0, ntot)

  for(i in 1:ntot){
   cluster_id[i] <- cluster_id_seq[which(cluster_id_orig[i] == cluster_id_unique)]
   subj_id[i] <- subj_id_seq[which(subj_id_orig[i] == subj_id_unique)]
  }

  dat <- list(ntot = ntot,
              subj_id = subj_id,
              y = as.array(y),
              p = ncol(x),
              x = x,
              k = nlevels(factor(y)))

  if(model != "fixed"){

    if(reff$U != "none"){
      dat$cluster_id = cluster_id
    }

    if(reff$U != "intercept" | reff$V != "intercept"){
      dat$n_ec <- length(unique(data[, ec_id]))
      dat$ind_ec <- ind_ec
      dat$nrepeat_ec <- nrepeat_ec
      dat$time <- as.array(data[, timeVar])
    }

    if(reff$U == "serial"){
      dat$kappa1 <- kappa[["kappa1"]]
    }

    if(reff$V == "serial"){
      dat$kappa2 <- kappa[["kappa2"]]
    }

  }

  if(model == "normal"){
    if(reff$U == "intercept" & reff$V == "serial"){
      mod <- stan_model(model_code = mixed_normal_ordinal_intercept_serial,
                        auto_write = TRUE)
    }
  }else if(model == "bridge"){
    if(reff$U == "serial" & reff$V == "serial"){
      mod <- stan_model(model_code = mixed_bridge_ordinal_serial_serial_B,
                        auto_write = TRUE)
    }else if(reff$U == "intercept" & reff$V == "serial"){
      mod <- stan_model(model_code = mixed_bridge_ordinal_intercept_serial,
                        auto_write = TRUE)
    }else if(reff$U == "serial" & reff$V == "intercept"){
      mod <- stan_model(model_code = mixed_bridge_ordinal_serial_intercept,
                        auto_write = TRUE)
    }
    else if(reff$U == "intercept" & reff$V == "intercept"){
      mod <- stan_model(model_code = mixed_bridge_ordinal_intercept_intercept,
                        auto_write = TRUE)
    }else if(reff$U == "none" & reff$V == "serial"){
      mod <- stan_model(model_code = mixed_bridge_ordinal_none_serial,
                        auto_write = TRUE)
    }else if(reff$U == "none" & reff$V == "intercept"){
      mod <- stan_model(model_code = mixed_bridge_ordinal_none_intercept,
                        auto_write = TRUE)
    }

  }else if(model == "fixed"){
    mod <- stan_model(model_code = ordinal_fixed,
                      auto_write = TRUE)
  }

  res <- sampling(mod, data = dat, ...)
  return(res)

}
