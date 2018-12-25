t_normal_ordinal_mixed_threelev_reparam_prior_on_var_inv_nu = "

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

parameters{
ordered[k - 1] alpha;
vector[p] beta;
vector[ncluster] ustar;
vector[nsubj] z_v;
real<lower = 0> sd_u;
real<lower = 0.01, upper = 0.5> inv_nu_u;
real<lower = 0> sigma_v;
}

transformed parameters{

real<lower = 0> sigma_u;
real<lower = 2, upper = 100> nu_u;
vector[ncluster] u;
vector[nsubj] v;
vector[ntot] v_vec;
vector[ntot] u_vec;
vector[ntot] linpred;

nu_u = 1/inv_nu_u;
sigma_u = sd_u * sqrt((nu_u - 2)/nu_u);

u = ustar * sigma_u;
v = z_v * sigma_v;

for(i in 1:ntot){
u_vec[i] = u[cluster_id[i]];
v_vec[i] = v[subj_id[i]];
}

linpred = x * beta + u_vec + v_vec;

}

model{

alpha ~ cauchy(0, 5);
beta ~ cauchy(0, 5);

inv_nu_u ~ beta(1, 1);

sd_u ~ cauchy(0, 5);
sigma_v ~ cauchy(0, 5);

ustar ~ student_t(nu_u, 0, 1);
z_v ~ normal(0, 1);

y ~ ordered_logistic(linpred, alpha);

}

"
