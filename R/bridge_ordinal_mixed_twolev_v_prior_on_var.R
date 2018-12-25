bridge_ordinal_mixed_twolev_v_prior_on_var = "

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
int<lower = 1> ntot;//total number of observations
int subj_id[ntot]; //id column - second
int<lower = 0> y[ntot]; // responses
int<lower = 1> p; // number of columns of the x matrix
matrix[ntot, p] x; // design matrix for fixed effects
int<lower = 1> nsubj; // number of individuals
int<lower = 3> k; //number of categories for the ordinal variable
int ind_s[nsubj, 2];
int nrepeat_s[nsubj];
}

parameters{
ordered[k - 1] alpha;
vector[p] beta;
vector[nsubj] v;
real<lower = 0> sd_v;
}

transformed parameters{
real<lower = 0, upper = 1> phi_v;
vector[ntot] v_vec;
vector[ntot] linpred;

phi_v = 1/sqrt( 3 * sd_v^2/(pi()^2) + 1);

//for(i in 1:ntot){
//v_vec[i] = v[subj_id[i]];
//}

for(i in 1:nsubj){
v_vec[ind_s[i, 1]:ind_s[i, 2]] = rep_vector(v[i], nrepeat_s[i]);
}

linpred = x * beta + v_vec;

}

model{

alpha ~ cauchy(0, 5);
beta ~ cauchy(0, 5);

sd_v ~ cauchy(0, 5);

v ~ bridge(phi_v);

y ~ ordered_logistic(linpred, alpha);

}

generated quantities{

vector[p] betamarg;
vector[k - 1] alphamarg;
vector[ntot] linpred2;
vector[ntot] log_lik;

alphamarg = alpha * phi_v;
betamarg = beta * phi_v;

linpred2 = x * beta + v_vec;
for(i in 1:ntot) log_lik[i] = ordered_logistic_lpmf(y[i] | linpred2[i], alpha);

}
"
