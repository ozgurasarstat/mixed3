mixed_bridge_ordinal_serial_intercept = "

functions{

//create cholesky decomposition of varcov matrix
matrix chol_cov(vector c_id,
                vector t,
                int mi,
                real sigma1,
                real delta1,
                real kappa1){

matrix[mi, mi] mat1;
matrix[mi, mi] L;

for(i in 1:mi){
for(j in 1:mi){

if(c_id[i] == c_id[j]){
mat1[i, j] = (sigma1^2) * exp(-(fabs(t[i] - t[j])/delta1)^kappa1);
}else{
mat1[i, j] = 0;
}

}
}

L = cholesky_decompose(mat1);

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
int<lower = 1> nsubj;
vector[ntot] cluster_id; // id column - first
int subj_id[ntot];
vector[ntot] time;       // time
int<lower = 0, upper = 5> y[ntot]; // responses
int<lower = 1> p;        // number of columns of the x matrix
matrix[ntot, p] x;       // desgin matrix for fixed effects
int<lower = 1> n_ec;     // number of extended families
int<lower = 3> k;        // number of categories for the ordinal variable
int ind_ec[n_ec, 2];     // matrix of indices for which observations belong to extended families
int nrepeat_ec[n_ec];    // number of observations that belong to extended families
real<lower = 0, upper = 2> kappa1;
}

parameters{
ordered[k - 1] alpha;
vector[p] beta;
vector[ntot] zstar;
vector[nsubj] v;
real<lower = 0> sigma1;
real<lower = 0> delta1;
real<lower = 0> sd_u;
real<lower = 0> sd_v;
}

transformed parameters{

vector[ntot] z;
vector[ntot] u_vec;
vector[ntot] v_vec;
real<lower = 0, upper = 1> phi_u;
real<lower = 0, upper = 1> phi_vstar;

phi_u = 1/sqrt(3*sd_u/(pi()^2) + 1);
phi_vstar = 1/sqrt( 3 * sd_v^2 * phi_u^2/(pi()^2) + 1);

for(i in 1:n_ec){
z[ind_ec[i, 1]:ind_ec[i, 2]] =
  chol_cov(cluster_id[ind_ec[i, 1]:ind_ec[i, 2]],
           time[ind_ec[i, 1]:ind_ec[i, 2]],
           nrepeat_ec[i],
           sigma1,
           delta1,
           kappa1) *
  zstar[ind_ec[i, 1]:ind_ec[i, 2]];
}

for(i in 1:ntot){
u_vec[i] = inv_cdf_bridge(Phi(z[i]), phi_u);
v_vec[i] = v[subj_id[i]];
}

}

model{

alpha ~ cauchy(0, 5);
beta ~ cauchy(0, 5);

sigma1 ~ cauchy(0, 5);
delta1 ~ cauchy(0, 5);

sd_u ~ cauchy(0, 5);
sd_v ~ cauchy(0, 5);

zstar ~ std_normal();
v ~ modified_bridge(phi_vstar, phi_u);

y ~ ordered_logistic(x * beta + u_vec + v_vec, alpha);

}

generated quantities{

real<lower = 0> sigmasq1;
vector[p] betamarg;
vector[k - 1] alphamarg;

sigmasq1 = sigma1^2;

alphamarg = alpha * phi_u * phi_vstar;
betamarg = beta * phi_u * phi_vstar;

}

"
