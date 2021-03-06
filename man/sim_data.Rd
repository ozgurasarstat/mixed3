% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_data.R
\docType{data}
\name{sim_data}
\alias{sim_data}
\title{Simulated data-set.}
\format{A data frame with 4000 rows and 5 variables:
\describe{
  \item{c_id}{cluster id}
  \item{s_id}{subject id}
  \item{t}{time}
  \item{cov}{covariate}
  \item{y}{outcome}
}}
\usage{
sim_data
}
\description{
A dataset containing three-level ordinal outcomes
for 250 clusters, 1000 subjects, 4 time points
}
\details{
{
 c_id: cluster id (1, ..., 250)
s_id: subject id (1, ..., 1000)
t: time (0, 1, 2, 3)
cov: a covariate simulated from normal
y: ordinal outcome taking 1, 2, or 3
The model for simulation is
logit[ P(Y_{ijk} <= a | x_{ijk}, U_i, V_{ij}) ] = alpha_a - t * beta_1 - cov * beta_2 - U_i - V_{ij},
where
- alpha_1 = -1, alpha_2 = 1
- beta_1 = 1, beta_2 = 0.5
- U_i = U_i^*/phi_V, U_i^* ~ Bridge(phi_{U^*}), phi_{U^*} = 0.8, phi_V = 0.6
- V_{ij} ~ Bridge(phi_{V})
}
}
\keyword{datasets}
