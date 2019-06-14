bridge_ordinal_mixed_threelev_gauss_copula = "

functions{

//pdf of bridge distribution
real bridge_lpdf(vector x, real theta){
int N;
real lpdf;
N = rows(x);
lpdf = 0;
for (n in 1:N)
  lpdf += log(0.5/pi()) + log(sin(theta * pi())) - log(cosh(theta * x[n]) + cos(theta * pi()));
return lpdf;
}

//pdf of modified bridge
real modified_bridge_lpdf(vector x, real theta1, real theta2){ // theta1: phi_{U*}, theat2: phi_{V}
int N;
real lpdf;
N = rows(x);
lpdf = 0;
for(n in 1:N)
  lpdf += log(0.5/pi()) + log(sin(theta1 * pi())) - log(cosh(theta1 * theta2 * x[n]) + cos(theta1 * pi())) + log(theta2);
return lpdf;
}

//inverse of the cdf of bridge
real inv_cdf_bridge(real x, real theta){
real out;
out = log( sin(theta * pi() * x) / sin(theta * pi() * (1 - x)) ) / theta;
return out;
}

//create cholesk decomposition of covariance matrix for z terms
matrix chol_covmat_exp(vector t, int mi, real delta, real kappa){

matrix[mi, mi] out;
matrix[mi, mi] L;

for (i in 1:(mi-1)){
out[i, i] = 1;
for (j in (i+1):mi){
out[i, j] = exp(-(fabs(t[i] - t[j])/delta)^kappa);
out[j, i] = exp(-(fabs(t[j] - t[i])/delta)^kappa);
}
}
out[mi, mi] = 1;

L = cholesky_decompose(out);
return L;
}

}

data{
int<lower = 1> ntot;//total number of observations
int cluster_id[ntot]; //id column - first
int subj_id[ntot]; //id column - second
int<lower = 0> y[ntot]; // responses
int<lower = 1> p; // number of columns of the x matrix
matrix[ntot, p] x; // design matrix for fixed effects
int<lower = 1> ncluster; // number of clusters
int<lower = 1> nsubj; // number of individuals
int<lower = 3> k; //number of categories for the ordinal variable
int ind_c[ncluster, 2];
int ind_s[nsubj, 2];
int nrepeat_c[ncluster];
int nrepeat_s[nsubj];
vector[ntot] time;
real<lower = 0, upper = 2> kappa;
}

parameters{
ordered[k - 1] alpha;
vector[p] beta;
vector[ncluster] u;
vector[ntot] zstar;
real<lower = 0> sd_u;
real<lower = 0> sd_v;
real<lower = 0> delta;
}

transformed parameters{
real<lower = 0, upper = 1> phi_ustar;
real<lower = 0, upper = 1> phi_v;
vector[ntot] v_vec;
vector[ntot] z_vec;
vector[ntot] u_vec;
vector[ntot] linpred;

phi_v     = 1/sqrt( 3 * sd_v^2/(pi()^2) + 1);
phi_ustar = 1/sqrt( 3 * sd_u^2 * phi_v^2/(pi()^2) + 1);

//for(i in 1:ntot){
//u_vec[i] = u[cluster_id[i]];
//v_vec[i] = v[subj_id[i]];
//}

for(i in 1:ncluster){
u_vec[ind_c[i, 1]:ind_c[i, 2]] = rep_vector(u[i], nrepeat_c[i]);
}

for(i in 1:nsubj){
z_vec[ind_s[i, 1]:ind_s[i, 2]] =
  chol_covmat_exp(time[ind_s[i, 1]:ind_s[i, 2]], nrepeat_s[i], delta, kappa) *
  zstar[ind_s[i, 1]:ind_s[i, 2]];
}

for(i in 1:ntot){
v_vec[i] = inv_cdf_bridge(Phi(z_vec[i]), phi_v);
}

linpred = x * beta + u_vec + v_vec;

}

model{

alpha ~ cauchy(0, 5);
beta ~ cauchy(0, 5);

sd_u ~ cauchy(0, 5);
sd_v ~ cauchy(0, 5);

zstar ~ std_normal();
u ~ modified_bridge(phi_ustar, phi_v);

delta ~ cauchy(0, 5);

y ~ ordered_logistic(linpred, alpha);

}

generated quantities{

vector[p] betamarg;
vector[k - 1] alphamarg;
vector[ntot] linpred2;
vector[ntot] log_lik;

alphamarg = alpha * phi_ustar * phi_v;
betamarg = beta * phi_ustar * phi_v;

linpred2 = x * beta + u_vec + v_vec;
for(i in 1:ntot) log_lik[i] = ordered_logistic_lpmf(y[i] | linpred2[i], alpha);

}
"
