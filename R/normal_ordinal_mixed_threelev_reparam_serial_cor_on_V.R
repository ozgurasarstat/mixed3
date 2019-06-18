normal_ordinal_mixed_threelev_reparam_serial_cor_on_V = "

functions{

//create cholesk decomposition of covariance matrix for z terms
matrix chol_covmat_exp(vector t, int mi, real sigma, real delta, real kappa){

matrix[mi, mi] out;
matrix[mi, mi] L;

for (i in 1:(mi-1)){
out[i, i] = sigma^2;
for (j in (i+1):mi){
out[i, j] = (sigma^2) * exp(-(fabs(t[i] - t[j])/delta)^kappa);
out[j, i] = (sigma^2) * exp(-(fabs(t[j] - t[i])/delta)^kappa);
}
}
out[mi, mi] = sigma^2;

L = cholesky_decompose(out);
return L;
}

}

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
vector[ntot] time;
real<lower = 0, upper = 2> kappa;
}

parameters{
ordered[k - 1] alpha;
vector[p] beta;
vector[ncluster] z_u;
vector[ntot] zstar_v;
real<lower = 0> sigma_u;
real<lower = 0> sigma_v;
real<lower = 0> delta;
}

transformed parameters{

vector[ncluster] u;
//vector[nsubj] v;
vector[ntot] v_vec;
vector[ntot] u_vec;
vector[ntot] linpred;

u = z_u * sigma_u;
//v = z_v * sigma_v;

//for(i in 1:ntot){
//u_vec[i] = u[cluster_id[i]];
//v_vec[i] = v[subj_id[i]];
//}

for(i in 1:ncluster){
u_vec[ind_c[i, 1]:ind_c[i, 2]] = rep_vector(u[i], nrepeat_c[i]);
}

for(i in 1:nsubj){
v_vec[ind_s[i, 1]:ind_s[i, 2]] =
  chol_covmat_exp(time[ind_s[i, 1]:ind_s[i, 2]], nrepeat_s[i], sigma_v, delta, kappa) *
  zstar_v[ind_s[i, 1]:ind_s[i, 2]];
}

linpred = x * beta + u_vec + v_vec;

}

model{

alpha ~ cauchy(0, 5);
beta ~ cauchy(0, 5);

sigma_u ~ cauchy(0, 5);
sigma_v ~ cauchy(0, 5);

z_u ~ normal(0, 1);
zstar_v ~ normal(0, 1);

delta ~ cauchy(0, 5);

y ~ ordered_logistic(linpred, alpha);

}

generated quantities{

vector[ntot] linpred2;
vector[ntot] log_lik;

linpred2 = x * beta + u_vec + v_vec;
for(i in 1:ntot) log_lik[i] = ordered_logistic_lpmf(y[i] | linpred2[i], alpha);

}

"
