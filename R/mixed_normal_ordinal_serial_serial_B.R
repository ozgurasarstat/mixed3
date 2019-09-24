mixed_normal_ordinal_serial_serial_B = "

functions{

//create cholesky decomposition of varcov matrix
matrix chol_cov(vector c_id,
                vector s_id,
                vector t,
                int mi,
                real sigma1,
                real sigma2,
                real delta1,
                real delta2,
                real kappa1,
                real kappa2){

matrix[mi, mi] mat1;
matrix[mi, mi] mat2;
matrix[mi, mi] L;

for(i in 1:mi){
for(j in 1:mi){

if(Sigma1_spec == 0){
if(c_id[i] == c_id[j]){
mat1[i, j] = (sigma1^2) * exp(-(fabs(t[i] - t[j])/delta1)^kappa1);
}else{
mat1[i, j] = 0;
}
}else if(Sigma1_spec == 1){
if(c_id[i] == c_id[j]){
mat1[i, j] = sigma1^2;
}else{
mat1[i, j] = 0;
}
}

if(Sigma2_spec == 0){
if(s_id[i] == s_id[j]){
mat2[i, j] = (sigma2^2) * exp(-(fabs(t[i] - t[j])/delta2)^kappa2);
}else{
mat2[i, j] = 0;
}
}else if(Sigma2_spec == 1){
if(s_id[i] == s_id[j]){
mat2[i, j] = sigma2^2;
}else{
mat2[i, j] = 0;
}
}

}
}

if(Sigma1_spec == 2){
L = cholesky_decompose(mat2);
}else{
L = cholesky_decompose(mat1 + mat2);
}

return L;
}

}

data{
int<lower = 1> ntot;     // total number of observations
vector[ntot] cluster_id; // id column - first
vector[ntot] subj_id;    // id column - second
vector[ntot] time;       // time
int<lower = 0, upper = 5> y[ntot]; // responses
int<lower = 1> p;        // number of columns of the x matrix
matrix[ntot, p] x;       // desgin matrix for fixed effects
int<lower = 1> n_ec;     // number of extended families
int<lower = 3> k;        // number of categories for the ordinal variable
int ind_ec[n_ec, 2];     // matrix of indices for which observations belong to extended families
int nrepeat_ec[n_ec];    // number of observations that belong to extended families
real<lower = 0, upper = 2> kappa1;
real<lower = 0, upper = 2> kappa2;
}

transformed data{
matrix[ntot, p] xc;
vector[p] xmeans;

for(i in 1:p){
xmeans[i] = mean(x[, i]);
xc[, i] = x[, i] - xmeans[i];
}

parameters{
ordered[k - 1] alpha_c;
vector[p] beta;
vector[ntot] zstar;
real<lower = 0> sigma1;
real<lower = 0> sigma2;
real<lower = 0> delta1;
real<lower = 0> delta2;
}

transformed parameters{

vector[ntot] b;

for(i in 1:n_ec){
b[ind_ec[i, 1]:ind_ec[i, 2]] =
  chol_cov(cluster_id[ind_ec[i, 1]:ind_ec[i, 2]],
           subj_id[ind_ec[i, 1]:ind_ec[i, 2]],
           time[ind_ec[i, 1]:ind_ec[i, 2]],
           nrepeat_ec[i],
           sigma1,
           sigma2,
           delta1,
           delta2,
           kappa1,
           kappa2) *
  zstar[ind_ec[i, 1]:ind_ec[i, 2]];
}

}

model{
vector[ntot] linpred = xc * beta + b;

alpha_c ~ cauchy(0, 5);
beta ~ cauchy(0, 5);

sigma1 ~ cauchy(0, 5);
sigma2 ~ cauchy(0, 5);

delta1 ~ cauchy(0, 5);
delta2 ~ cauchy(0, 5);

zstar ~ std_normal();

for(i in 1:ntot){
target += ordered_logistic_lpmf(y[i] | linpred[i], alpha_c);
}

}

generated quantities{

vector[k - 1] alpha = alpha_c + dot_product(xmeans, beta);

real<lower = 0> sigmasq1 = sigma1^2;
real<lower = 0> sigmasq2 = sigma2^2;

}

"
