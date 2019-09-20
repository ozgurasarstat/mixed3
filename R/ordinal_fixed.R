ordinal_fixed = "

data{
int<lower = 1> ntot;//total number of observations
int<lower = 0> y[ntot]; // responses
int<lower = 1> p; // number of columns of the x matrix
matrix[ntot, p] x; // design matrix for fixed effects
int<lower = 3> k; //number of categories for the ordinal variable
}

//center the design matrix
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
}

model{
vector[ntot] linpred = xc * beta;

alpha_c ~ cauchy(0, 5);
beta ~ cauchy(0, 5);

for(i in 1:ntot){
target += ordered_logistic_lpmf(y[i] | linpred[i], alpha_c);
}

}

generated quantities{
vector[k - 1] alpha = alpha_c + dot_product(xmeans, beta);
}

"
