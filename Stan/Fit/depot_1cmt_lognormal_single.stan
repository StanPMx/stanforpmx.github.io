// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model
// Single subject
// lognormal error - DV = CP*exp(eps)
// Closed form solution using a self-written function

functions{
  
  real depot_1cmt(real dose, real cl, real v, real ka, 
                  real time_since_dose){
    
    real ke = cl/v;
    
    real cp = dose/v * ka/(ka - ke) * 
              (exp(-ke*time_since_dose) - exp(-ka*time_since_dose));
    
    return cp;
    
  }
  
}

data{
  
  int n_obs;
  real<lower = 0> dose;
  array[n_obs] real time;
  real time_of_first_dose;
  vector[n_obs] dv;
  
  real<lower = 0> scale_cl;     // Prior Scale parameter for CL
  real<lower = 0> scale_v;      // Prior Scale parameter for V
  real<lower = 0> scale_ka;     // Prior Scale parameter for KA
  
  real<lower = 0> scale_sigma;  // Prior Scale parameter for lognormal error
  
  int n_pred;                   // Number of new times at which to make a prediction
  array[n_pred] real time_pred; // New times at which to make a prediction
 
}
transformed data{ 
  
  vector[n_obs] time_since_dose = to_vector(time) - time_of_first_dose;
  vector[n_pred] time_since_dose_pred = to_vector(time_pred) - 
                                        time_of_first_dose;
  
}
parameters{  
  
  real<lower = 0> CL;
  real<lower = 0> V;
  real<lower = CL/V> KA;
  
  real<lower = 0> sigma;
  
}
transformed parameters{

  vector[n_obs] ipred;
  
  for(i in 1:n_obs){
    ipred[i] = depot_1cmt(dose, CL, V, KA, time_since_dose[i]);
  }
  
}

model{ 
  
  // Priors
  CL ~ cauchy(0, scale_cl);
  V ~ cauchy(0, scale_v);
  KA ~ normal(0, scale_ka) T[CL/V, ];
  
  sigma ~ normal(0, scale_sigma);
  
  // Likelihood
  dv ~ lognormal(log(ipred), sigma);

}
generated quantities{
  
  real<lower = 0> KE = CL/V;
  real<lower = 0> sigma_sq = square(sigma);

  vector[n_obs] dv_ppc;
  vector[n_obs] log_lik;
  vector[n_pred] cp;
  vector[n_pred] dv_pred;

  vector[n_obs] ires = log(dv) - log(ipred);
  vector[n_obs] iwres = ires/sigma;
  
  for(i in 1:n_obs){
    dv_ppc[i] = lognormal_rng(log(ipred[i]), sigma);
    log_lik[i] = lognormal_lpdf(dv[i] | log(ipred[i]), sigma);
  }

  for(j in 1:n_pred){
    if(time_since_dose_pred[j] <= 0){
      cp[j] = 0;
      dv_pred[j] = 0;
    }else{
      cp[j] = depot_1cmt(dose, CL, V, KA, time_since_dose_pred[j]);
      dv_pred[j] = lognormal_rng(log(cp[j]), sigma);
    }
  }
}

