
## set the working directory
setwd("C:\\Users\\Ozgur Asar\\Desktop\\Desktop\\dropbox\\on-going works\\Gelir ve Yasam Kosullari\\Submit - Statistical Modelling\\SILC - Statistical Modelling\\1st revision\\SILC analysis")

## import the data-set
data <- read.csv("comb_10_13_processed.csv", header = TRUE, stringsAsFactors = FALSE)

data$y         <- with(data, ifelse(FS010_3lev_f == "poor", 3, ifelse(FS010_3lev_f == "fair", 2, 1)))
data$year      <- factor(data$FB010_f, levels = c("2010", "2011", "2012", "2013"))
data$log.mhdi  <- log(data$HG110_adj_h + 1)
data$gender    <- factor(data$FK090_fk_ch, levels = c("male", "female"))
data$marital   <- factor(data$FB100_new_f, levels = c("married", "never_married", "others"))
data$age       <- factor(data$FK070_new_fk, levels = c(">34 & <=64", "<=34", ">=65"))
data$education <- factor(data$FE030_new_f, levels = c("higher_edu", "primary_less", "secondary_high"))
data$work      <- factor(data$FI010_new_f, levels = c("full_part", "housekeeper", "retired", "student", "unemployed", "other"))
data$age2      <- scale(data$FK070_fk)
data$work2     <- ifelse(data$FI010_new_f == "unemployed", 0, 1)

data <- data[order(data$HKIMLIK_f, data$FKIMLIK_f), ]

data <- data.frame(GHKIMLIK_f = substr(data$HKIMLIK_f, 1, 5), data, stringsAsFactors = FALSE)

data$time <- data$FB010_f - 2010

data_sub <- data[1:2000, ]

library(mixed3)

fit <- wrap_B(formula = y ~ as.factor(gender) + log.mhdi,
              data = data_sub,
              ec_id = "GHKIMLIK_f",
              c_id = "HKIMLIK_f",
              s_id = "FKIMLIK_f",
              timeVar = "time",
              kappa = list(kappa1 = 1, kappa2 = 1),
              model = "bridge",
              reff = list(U = "intercept", V = "serial"),
              iter = 1000,
              warmup = 500,
              chains = 4,
              cores = 4)

print(fit, pars = c("alpha", "alphamarg", "beta", "betamarg",
                    "phi", "sigma1", "delta1", "sigma2", "delta2"))
traceplot(fit, pars = c("alpha", "alphamarg", "beta", "betamarg",
                    "phi", "sigma1", "delta1", "sigma2", "delta2"))
