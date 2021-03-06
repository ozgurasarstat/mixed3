
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

# below: try model = "fixed", model = "bridge_twolev", model = "bridge_threelev",
#            model = "normal_threelev"

res_fixed <-  wrapper(formula = y ~ t + cov,
                   data = sim_data,
                   c_id = "c_id",
                   s_id = "s_id",
                   model = "fixed",
                   timeVar = "t",
                   pars = c("alpha", "alpha_c",  "beta"),
                   iter = 2000,
                   chains = 4,
                   warmup = 1000,
                   cores = 4)

print(res_fixed, pars = c("alpha", "alpha_c", "beta"), digits = 3)

res_br2 <-  wrapper(formula = y ~ t + cov,
                    data = sim_data,
                    c_id = "c_id",
                    s_id = "s_id",
                    model = "bridge_twolev",
                    timeVar = "t",
                    pars = c("alpha", "alpha_c", "alphamarg", "beta", "betamarg", "v_vec", "phi_v"),
                    iter = 2000,
                    chains = 4,
                    warmup = 1000,
                    cores = 4)

print(res_br2, pars = c("alpha", "alpha_c", "alphamarg", "beta", "betamarg", "phi_v"), digits = 3)

res_br3 <-  wrapper(formula = y ~ t + cov,
                    data = sim_data,
                    c_id = "c_id",
                    s_id = "s_id",
                    model = "bridge_threelev",
                    timeVar = "t",
                    pars = c("alpha", "alpha_c", "alphamarg", "betamarg", "beta", "u_vec", "v_vec", "phi_ustar", "phi_v"),
                    iter = 2000,
                    chains = 4,
                    warmup = 1000,
                    cores = 4)

print(res_br3, pars = c("alpha", "alpha_c", "alphamarg", "beta", "betamarg", "phi_ustar", "phi_v"), digits = 3)


res_nor <-  wrapper(formula = y ~ t + cov,
                    data = sim_data,
                    c_id = "c_id",
                    s_id = "s_id",
                    model = "normal_threelev",
                    timeVar = "t",
                    pars = c("alpha", "alpha_c", "beta", "u_vec", "v_vec"),
                    iter = 2000,
                    chains = 4,
                    warmup = 1000,
                    cores = 4)

print(res_nor, pars = c("alpha", "alpha_c", "beta"), digits = 3)


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
