
##
## Ozgur Asar (ozgurasarstat@gmail.com)
## Bayesian analysis of Turkish Income and Living Conditions data, using clustered
## longitudinal ordinal modelling with Bridge distributed random-effects
## 19.03.2019
##

## Install the package from github
install.packages("devtools")
library(devtools)

install_github("ozgurasarstat/mixed3")
library(mixed3)

## load the data-set
data(sim_data)

# below: try model = "normal", model = "t", model = "normal_t", model = "t_normal",
#            model = "normal_t_no_sigma_v", model = "two_bridge", model = "fixed"
res_br <-  three_level_mixed_wrapper(formula = y ~ t + cov,
                                     data = sim_data,
                                     c_id = "c_id",
                                     s_id = "s_id",
                                     model = "normal_serial_cor_on_V",
                                     timeVar = "t",
                                     kappa = 1,
                                     pars = c("alpha", "beta",  "delta"),
                                     iter = 2000,   # increase?
                                     chains = 4,    # increase?
                                     warmup = 1000, # increase?
                                     cores = 4)     # increase?

# tells the parameter names - ignore v_vec, u_vec, linpred, lp__
print(names(extract(res_br)))

# summarise results
print(res_br, pars = c("alpha", "beta", "alphamarg", "betamarg", "phi_ustar", "phi_v"), digits = 3)

# traceplots
traceplot(res_br, pars = c("alpha"), inc_warmup = TRUE)
traceplot(res_br, pars = c("alphamarg"), inc_warmup = TRUE)
traceplot(res_br, pars = c("beta"), inc_warmup = TRUE)
traceplot(res_br, pars = c("betamarg"), inc_warmup = TRUE)
traceplot(res_br, pars = c("phi_ustar", "phi_v"), inc_warmup = TRUE)

# checking chains
check_hmc_diagnostics(res_br)
