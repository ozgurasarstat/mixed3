ordinal_fixed = "

data{
int<lower = 1> ntot;//total number of observations
int<lower = 0> y[ntot]; // responses
int<lower = 1> p; // number of columns of the x matrix
matrix[ntot, p] x; // design matrix for fixed effects
int<lower = 3> k; //number of categories for the ordinal variable
}

parameters{
ordered[k - 1] alpha;
vector[p] beta;
}

transformed parameters{

vector[ntot] linpred;

linpred = x * beta;

}

model{

alpha ~ cauchy(0, 5);
beta ~ cauchy(0, 5);

y ~ ordered_logistic(linpred, alpha);

}

generated quantities{

vector[ntot] linpred2;
vector[ntot] log_lik;

linpred2 = x * beta;
for(i in 1:ntot) log_lik[i] = ordered_logistic_lpmf(y[i] | linpred2[i], alpha);

}

"
