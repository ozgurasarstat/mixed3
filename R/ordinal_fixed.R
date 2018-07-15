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

for (n in 1:ntot)
y[n] ~ ordered_logistic(linpred[n], alpha);

}

"