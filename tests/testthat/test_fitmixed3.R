
##
## Ozgur Asar (ozgurasarstat@gmail.com)
## Bayesian analysis of Turkish Income and Living Conditions data, using clustered
## longitudinal ordinal modelling with Bridge distributed random-effects
## 19.03.2019
##

## Install the package from github
library(mixed3)

## load the data-set
data(sim_data)

# below: try model = "normal", model = "t", model = "normal_t", model = "t_normal",
#            model = "normal_t_no_sigma_v", model = "two_bridge", model = "fixed"

fit <-  fit_mixed3(formula = y ~ t + cov,
                   data = sim_data,
                   ec_id = "c_id",
                   c_id = "c_id",
                   s_id = "s_id",
                   timeVar = "t",
                   model = "normal",
                   ref = list(U = "intercept", V = "intercept"),
                   #pars = c("alpha", "alpha_c",  "beta"),
                   iter = 2000,   # increase?
                   chains = 1,    # increase?
                   warmup = 1000, # increase?
                   cores = 4)     # increase?

# tells the parameter names - ignore v_vec, u_vec, linpred, lp__
print(names(extract(fit)))

# summarise results
print(fit, pars = c("alpha", "beta", "sigmasq2", "delta2"), digits = 3)

# traceplots
traceplot(res_br, pars = c("alpha"), inc_warmup = TRUE)
traceplot(res_br, pars = c("alphamarg"), inc_warmup = TRUE)
traceplot(res_br, pars = c("beta"), inc_warmup = TRUE)
traceplot(res_br, pars = c("betamarg"), inc_warmup = TRUE)
traceplot(res_br, pars = c("phi_ustar", "phi_v"), inc_warmup = TRUE)

# checking chains
check_hmc_diagnostics(res_br)
