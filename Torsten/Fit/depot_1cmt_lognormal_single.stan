// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model
// Single subject
// lognormal error - DV = CP*exp(eps)
// Closed form solution using a Torsten function

data{
  
  int n_total;
  int n_obs;
  array[n_obs] int i_obs;
  array[n_total] real amt;
  array[n_total] int cmt;
  array[n_total] int evid;
  array[n_total] real rate;
  array[n_total] real ii;
  array[n_total] int addl;
  array[n_total] int ss;
  array[n_total] real time;
  vector<lower = 0>[n_total] dv;
  
  real<lower = 0> scale_cl;     // Prior Scale parameter for CL
  real<lower = 0> scale_v;      // Prior Scale parameter for V
  real<lower = 0> scale_ka;     // Prior Scale parameter for KA
  
  real<lower = 0> scale_sigma;  // Prior Scale parameter for lognormal error
  
  // These are data variables needed to make predictions at unobserved 
  // timepoints
  int n_pred;               // Number of new times at which to make a prediction
  array[n_pred] real amt_pred;
  array[n_pred] int cmt_pred;
  array[n_pred] int evid_pred;
  array[n_pred] real rate_pred;
  array[n_pred] real ii_pred;
  array[n_pred] int addl_pred;
  array[n_pred] int ss_pred;
  array[n_pred] real time_pred;

}
transformed data{
  
  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  
}
parameters{  
  
  real<lower = 0> CL;
  real<lower = 0> V;
  real<lower = CL/V> KA;
  
  real<lower = 0> sigma;
  
}
transformed parameters{
  
  vector[n_obs] ipred;
   
  { 
    vector[n_total] dv_ipred;
    matrix[n_total, 2] x_ipred = pmx_solve_onecpt(time,
                                                  amt,
                                                  rate,
                                                  ii,
                                                  evid,
                                                  cmt,
                                                  addl,
                                                  ss,
                                                  {CL, V, KA})';
                                                  
    dv_ipred = x_ipred[, 2] ./ V;
    ipred = dv_ipred[i_obs];
  }
  
}
model{ 

  // Priors
  CL ~ cauchy(0, scale_cl);
  V ~ cauchy(0, scale_v);
  KA ~ normal(0, scale_ka) T[CL/V, ];
  
  sigma ~ normal(0, scale_sigma);
  
  // Likelihood
  dv_obs ~ lognormal(log(ipred), sigma);

}
generated quantities{

  real<lower = 0> KE = CL/V;
  real<lower = 0> sigma_sq = square(sigma);

  vector[n_obs] dv_ppc;
  vector[n_obs] log_lik;
  vector[n_pred] cp;
  vector[n_pred] dv_pred;

  vector[n_obs] ires = log(dv_obs) - log(ipred);
  vector[n_obs] iwres = ires/sigma;
  
  for(i in 1:n_obs){
    dv_ppc[i] = lognormal_rng(log(ipred[i]), sigma);
    log_lik[i] = lognormal_lpdf(dv[i] | log(ipred[i]), sigma);
  }

  {
    
    matrix[n_pred, 2] x_cp;
    array[3] real theta_params = {CL, V, KA}; 
    
    x_cp = pmx_solve_onecpt(time_pred,
                            amt_pred,
                            rate_pred,
                            ii_pred,
                            evid_pred,
                            cmt_pred,
                            addl_pred,
                            ss_pred,
                            theta_params)';
                               
    cp = x_cp[, 2] ./ V;
    
  }

  for(i in 1:n_pred){
    if(cp[i] == 0){
      dv_pred[i] = 0;
    }else{
      dv_pred[i] = lognormal_rng(log(cp[i]), sigma);  
    }
    
  }
}



     