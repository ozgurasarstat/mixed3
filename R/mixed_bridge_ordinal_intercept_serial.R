mixed_bridge_ordinal_intercept_serial = "

functions{

//create cholesky decomposition of varcov matrix
matrix chol_cov(vector s_id,
                vector t,
                int mi,
                real sigma2sq,
                real delta2,
                real kappa2){

matrix[mi, mi] mat2;
matrix[mi, mi] L;

for(i in 1:mi){
for(j in 1:mi){

if(s_id[i] == s_id[j]){
if(i == j){
mat2[i, j] = sigma2sq;
}else{
mat2[i, j] = sigma2sq * exp(-(fabs(t[i] - t[j])/delta2)^kappa2);
}
}else{
mat2[i, j] = 0;
}

}
}

L = cholesky_decompose(mat2);
return L;
}

//inverse of the cdf of bridge
real inv_cdf_bridge(real x, real theta){
real out;
out = log( sin(theta * pi() * x) / sin(theta * pi() * (1 - x)) ) / theta;
return out;
}

real modified_bridge_lpdf(vector x, real theta1, real theta2){ // theta1: phi_{V*}, theat2: phi_{U}
int N;
real lpdf;
N = rows(x);
lpdf = 0;
for(n in 1:N)
  lpdf += log(0.5/pi()) + log(sin(theta1 * pi())) - log(cosh(theta1 * theta2 * x[n]) + cos(theta1 * pi())) + log(theta2);
return lpdf;
}

}

data{
int<lower = 1> ntot;     // total number of observations
int<lower = 1> ncluster;
int cluster_id[ntot]; // id column - first
vector[ntot] subj_id;    // id column - second
vector[ntot] time;       // time
int<lower = 0, upper = 5> y[ntot]; // responses
int<lower = 1> p;        // number of columns of the x matrix
matrix[ntot, p] x;       // desgin matrix for fixed effects
int<lower = 1> n_ec;     // number of extended families
int<lower = 3> k;        // number of categories for the ordinal variable
int ind_ec[n_ec, 2];     // matrix of indices for which observations belong to extended families
int nrepeat_ec[n_ec];    // number of observations that belong to extended families
real<lower = 0, upper = 2> kappa2;
}

transformed data{
matrix[ntot, p] xc;
vector[p] xmeans;

for(i in 1:p){
xmeans[i] = mean(x[, i]);
xc[, i] = x[, i] - xmeans[i];
}

}

parameters{
ordered[k - 1] alpha_c;
vector[p] beta;
vector[ncluster] u;
vector[ntot] zstar;
real<lower = 0> sigma2;
real<lower = 0> delta2;
real<lower = 0> sd_u;
real<lower = 0> sd_v;
}

transformed parameters{

vector[ntot] u_vec;
vector[ntot] z;
vector[ntot] v_vec;
real<lower = 0, upper = 1> phi_ustar;
real<lower = 0, upper = 1> phi_v;
real<lower = 0> sigma2sq;

phi_v     = 1/sqrt( 3 * sd_v^2/(pi()^2) + 1);
phi_ustar = 1/sqrt( 3 * sd_u^2 * phi_v^2/(pi()^2) + 1);

sigma2sq = sigma2^2;

for(i in 1:n_ec){
z[ind_ec[i, 1]:ind_ec[i, 2]] =
  chol_cov(subj_id[ind_ec[i, 1]:ind_ec[i, 2]],
           time[ind_ec[i, 1]:ind_ec[i, 2]],
           nrepeat_ec[i],
           sigma2sq,
           delta2,
           kappa2) *
  zstar[ind_ec[i, 1]:ind_ec[i, 2]];
}

for(i in 1:ntot){
v_vec[i] = inv_cdf_bridge(Phi(z[i]), phi_v);
u_vec[i] = u[cluster_id[i]];
}

}

model{

vector[ntot] linpred = xc * beta + u_vec + v_vec;

alpha_c ~ cauchy(0, 5);
beta ~ cauchy(0, 5);

sigma2 ~ cauchy(0, 5);
delta2 ~ cauchy(0, 5);

sd_u ~ cauchy(0, 5);
sd_v ~ cauchy(0, 5);

zstar ~ std_normal();
u ~ modified_bridge(phi_ustar, phi_v);

for(i in 1:ntot){
target += ordered_logistic_lpmf(y[i] | linpred[i], alpha_c);
}

}

generated quantities{

vector[k - 1] alpha = alpha_c + dot_product(xmeans, beta);
vector[k - 1] alphamarg = alpha * phi_ustar * phi_v;
vector[p] betamarg = beta * phi_ustar * phi_v;

}

"
