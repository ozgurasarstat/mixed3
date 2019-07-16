mixed_normal_ordinal_intercept_serial = "

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

parameters{
ordered[k - 1] alpha;
vector[p] beta;
vector[ncluster] ustar;
vector[ntot] vstar;
real<lower = 0> sigma1;
real<lower = 0> sigma2;
real<lower = 0> delta2;
}

transformed parameters{

vector[ntot] u_vec;
vector[ntot] v_vec;

for(i in 1:n_ec){
v_vec[ind_ec[i, 1]:ind_ec[i, 2]] =
  chol_cov(subj_id[ind_ec[i, 1]:ind_ec[i, 2]],
           time[ind_ec[i, 1]:ind_ec[i, 2]],
           nrepeat_ec[i],
           sigma2,
           delta2,
           kappa2) *
  vstar[ind_ec[i, 1]:ind_ec[i, 2]];
}

for(i in 1:ntot){
u_vec[i] = sigma1 * ustar[cluster_id[i]];
}

}

model{

alpha ~ cauchy(0, 5);
beta ~ cauchy(0, 5);

sigma1 ~ cauchy(0, 5);
sigma2 ~ cauchy(0, 5);

delta2 ~ cauchy(0, 5);

vstar ~ std_normal();
ustar ~ std_normal();

y ~ ordered_logistic(x * beta + u_vec + v_vec, alpha);

}

generated quantities{

real<lower = 0> sigmasq1;
real<lower = 0> sigmasq2;

sigmasq1 = sigma1^2;
sigmasq2 = sigma2^2;

}

"
