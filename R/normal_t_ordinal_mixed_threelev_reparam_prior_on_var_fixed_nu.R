normal_t_ordinal_mixed_threelev_reparam_prior_on_var_fixed_nu = "

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
int ind_c[ncluster, 2];
int ind_s[nsubj, 2];
int nrepeat_c[ncluster];
int nrepeat_s[nsubj];
real<lower = 2> nu_v;
}

parameters{
ordered[k - 1] alpha;
vector[p] beta;
vector[ncluster] ustar;
vector[nsubj] vstar;
real<lower = 0> sigma_u;
//real<lower = 0.01, upper = 0.5> inv_nu_u;
real<lower = 0> sd_v;
//real<lower = 0.01, upper = 0.5> inv_nu_v;
}

transformed parameters{

//real<lower = 0> sigma_u;
real<lower = 0> sigma_v;
//real<lower = 2, upper = 100> nu_u;
//real<lower = 2, upper = 100> nu_v;
vector[ncluster] u;
vector[nsubj] v;
vector[ntot] v_vec;
vector[ntot] u_vec;
vector[ntot] linpred;

//nu_u = 1/inv_nu_u;
//nu_v = 1/inv_nu_v;

//sigma_u = sd_u * sqrt((nu_u-2)/nu_u);
sigma_v = sd_v * sqrt((nu_v-2)/nu_v);

u = ustar * sigma_u;
v = vstar * sigma_v;

//for(i in 1:ntot){
//u_vec[i] = u[cluster_id[i]];
//v_vec[i] = v[subj_id[i]];
//}

for(i in 1:ncluster){
u_vec[ind_c[i, 1]:ind_c[i, 2]] = rep_vector(u[i], nrepeat_c[i]);
}

for(i in 1:nsubj){
v_vec[ind_s[i, 1]:ind_s[i, 2]] = rep_vector(v[i], nrepeat_s[i]);
}

linpred = x * beta + u_vec + v_vec;

}

model{

alpha ~ cauchy(0, 5);
beta ~ cauchy(0, 5);

//inv_nu_u ~ beta(1, 1);
//inv_nu_v ~ beta(1, 1);

sigma_u ~ cauchy(0, 5);
sd_v ~ cauchy(0, 5);

ustar ~ normal(0, 1);
vstar ~ student_t(nu_v, 0, 1);

y ~ ordered_logistic(linpred, alpha);

}

generated quantities{

vector[ntot] linpred2;
vector[ntot] log_lik;

linpred2 = x * beta + u_vec + v_vec;
for(i in 1:ntot) log_lik[i] = ordered_logistic_lpmf(y[i] | linpred2[i], alpha);

}

"
