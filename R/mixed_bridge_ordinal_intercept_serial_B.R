mixed_bridge_ordinal_intercept_serial_B = "

functions{

//create cholesky decomposition of varcov matrix
matrix chol_cov(vector c_id,
                vector s_id,
                vector t,
                int mi,
                real sigma1,
                real sigma2,
                real delta2,
                real kappa2){

matrix[mi, mi] mat1;
matrix[mi, mi] mat2;
matrix[mi, mi] L;

for(i in 1:mi){
for(j in 1:mi){

if(c_id[i] == c_id[j]){
mat1[i, j] = (sigma1^2);
}else{
mat1[i, j] = 0;
}

if(s_id[i] == s_id[j]){
mat2[i, j] = (sigma2^2) * exp(-(fabs(t[i] - t[j])/delta2)^kappa2);
}else{
mat2[i, j] = 0;
}

}
}

L = cholesky_decompose(mat1 + mat2);

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
real<lower = 0, upper = 2> kappa2;
}

parameters{
ordered[k - 1] alpha;
vector[p] beta;
vector[ntot] zstar;
real<lower = 0> sigma1;
real<lower = 0> sigma2;
real<lower = 0> delta2;
real<lower = 0> sd_b;
}

transformed parameters{

vector[ntot] z;
vector[ntot] b;
real<lower = 0, upper = 1> phi;

phi = 1/sqrt(3*sd_b/(pi()^2) + 1);

for(i in 1:n_ec){
z[ind_ec[i, 1]:ind_ec[i, 2]] =
  chol_cov(cluster_id[ind_ec[i, 1]:ind_ec[i, 2]],
           subj_id[ind_ec[i, 1]:ind_ec[i, 2]],
           time[ind_ec[i, 1]:ind_ec[i, 2]],
           nrepeat_ec[i],
           sigma1,
           sigma2,
           delta2,
           kappa2) *
  zstar[ind_ec[i, 1]:ind_ec[i, 2]];
}

for(i in 1:ntot){
b[i] = inv_cdf_bridge(Phi(z[i]), phi);
}

}

model{

alpha ~ cauchy(0, 5);
beta ~ cauchy(0, 5);

sigma1 ~ cauchy(0, 5);
sigma2 ~ cauchy(0, 5);

delta2 ~ cauchy(0, 5);

sd_b ~ cauchy(0, 5);

zstar ~ std_normal();

y ~ ordered_logistic(x * beta + b, alpha);

}

generated quantities{

real<lower = 0> sigmasq1;
real<lower = 0> sigmasq2;
vector[p] betamarg;
vector[k - 1] alphamarg;

sigmasq1 = sigma1^2;
sigmasq2 = sigma2^2;

alphamarg = alpha * phi;
betamarg = beta * phi;

}

"
