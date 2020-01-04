normal_threelev = "

data{
int<lower = 1> ntot;//total number of observations
int cluster_id[ntot]; //id column - first
int subj_id[ntot]; //id column - second
int<lower = 0, upper = 5> y[ntot]; // responses
int<lower = 1> p; // number of columns of the x matrix
matrix[ntot, p] x; // desgin matrix for fixed effects
int<lower = 1> ncluster; // number of clusters
int<lower = 1> nsubj; // number of individuals
int<lower = 3> k; //number of categories for the ordinal variable
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
vector[ncluster] z_u;
vector[nsubj] z_v;
real<lower = 0> sigma_u;
real<lower = 0> sigma_v;
}

transformed parameters{

vector[ncluster] u;
vector[nsubj] v;
vector[ntot] v_vec;
vector[ntot] u_vec;

u = z_u * sigma_u;
v = z_v * sigma_v;

for(i in 1:ntot){
u_vec[i] = u[cluster_id[i]];
v_vec[i] = v[subj_id[i]];
}

}

model{

vector[ntot] linpred = xc * beta + u_vec + v_vec;

alpha ~ cauchy(0, 5);
beta ~ cauchy(0, 5);

sigma_u ~ cauchy(0, 5);
sigma_v ~ cauchy(0, 5);

z_u ~ normal(0, 1);
z_v ~ normal(0, 1);

for(i in 1:ntot){
target += ordered_logistic_lpmf(y[i] | linpred[i], alpha_c);
}

}

generated quantities{

vector[k - 1] alpha = alpha_c + dot_product(xmeans, beta);

}
"
