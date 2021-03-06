mixed_bridge_ordinal_none_intercept = "

functions{

real bridge_lpdf(vector x, real theta){
int N;
real lpdf;
N = rows(x);
lpdf = 0;
for (n in 1:N)
  lpdf += log(0.5/pi()) + log(sin(theta * pi())) - log(cosh(theta * x[n]) + cos(theta * pi()));
return lpdf;
}

}

data{
int<lower = 1> ntot;     // total number of observations
int<lower = 1> nsubj;
int subj_id[ntot];    // id column - second
int<lower = 0, upper = 5> y[ntot]; // responses
int<lower = 1> p;        // number of columns of the x matrix
matrix[ntot, p] x;       // desgin matrix for fixed effects
int<lower = 3> k;        // number of categories for the ordinal variable
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
vector[nsubj] v;
real<lower = 0> sd_v;
}

transformed parameters{

vector[ntot] v_vec;
real<lower = 0, upper = 1> phi_v;

phi_v     = 1/sqrt( 3 * sd_v^2/(pi()^2) + 1);

for(i in 1:ntot){
v_vec[i] = v[subj_id[i]];
}

}

model{

vector[ntot] linpred = xc * beta + v_vec;

alpha_c ~ cauchy(0, 5);
beta ~ cauchy(0, 5);

sd_v ~ cauchy(0, 5);

v ~ bridge(phi_v);

for(i in 1:ntot){
target += ordered_logistic_lpmf(y[i] | linpred[i], alpha_c);
}

}

generated quantities{

vector[k - 1] alpha = alpha_c + dot_product(xmeans, beta);
vector[k - 1] alphamarg = alpha * phi_v;
vector[p] betamarg = beta * phi_v;

}

"
