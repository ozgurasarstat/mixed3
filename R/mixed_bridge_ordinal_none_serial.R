mixed_bridge_ordinal_none_serial = "

functions{

//create cholesky decomposition of varcov matrix
matrix chol_cov(vector s_id,
                vector t,
                int mi,
                real sigma2,
                real delta2,
                real kappa2){

matrix[mi, mi] mat2;
matrix[mi, mi] L;

for(i in 1:mi){
for(j in 1:mi){

if(s_id[i] == s_id[j]){
mat2[i, j] = (sigma2^2) * exp(-(fabs(t[i] - t[j])/delta2)^kappa2);
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

}

data{
int<lower = 1> ntot;     // total number of observations
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
vector[ntot] zstar;
real<lower = 0> sigma2;
real<lower = 0> delta2;
real<lower = 0> sd_v;
}

transformed parameters{

vector[ntot] z;
vector[ntot] v_vec;
real<lower = 0, upper = 1> phi_v;

phi_v = 1/sqrt(3*sd_v/(pi()^2) + 1);

for(i in 1:n_ec){
z[ind_ec[i, 1]:ind_ec[i, 2]] =
  chol_cov(subj_id[ind_ec[i, 1]:ind_ec[i, 2]],
           time[ind_ec[i, 1]:ind_ec[i, 2]],
           nrepeat_ec[i],
           sigma2,
           delta2,
           kappa2) *
  zstar[ind_ec[i, 1]:ind_ec[i, 2]];
}

for(i in 1:ntot){
v_vec[i] = inv_cdf_bridge(Phi(z[i]), phi_v);
}

}

model{

vector[ntot] linpred = xc * beta + v_vec;

alpha_c ~ cauchy(0, 5);
beta ~ cauchy(0, 5);

sigma2 ~ cauchy(0, 5);

delta2 ~ cauchy(0, 5);

sd_v ~ cauchy(0, 5);

zstar ~ std_normal();

for(i in 1:ntot){
target += ordered_logistic_lpmf(y[i] | linpred[i], alpha_c);
}

}

generated quantities{

vector[k - 1] alpha = alpha_c + dot_product(xmeans, beta);

real<lower = 0> sigmasq2 = sigma2^2;
vector[p] betamarg = alpha * phi_v;
vector[k - 1] alphamarg = beta * phi_v;

}

"
