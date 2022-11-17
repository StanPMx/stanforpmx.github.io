// IV infusion
// Two-compartment PK Model
// IIV on CL, VC, Q, VP (full covariance matrix)
// proportional plus additive error - DV = CP(1 + eps_p) + eps_a
// Closed form solution using Torsten

data{
  
  int n_subjects;
  int n_total;
  int n_obs;
  array[n_obs] int i_obs;
  array[n_total] int ID;
  array[n_total] real amt;
  array[n_total] int cmt;
  array[n_total] int evid;
  array[n_total] real rate;
  array[n_total] real ii;
  array[n_total] int addl;
  array[n_total] int ss;
  array[n_total] real time;
  vector[n_total] dv;
  array[n_subjects] int subj_start;
  array[n_subjects] int subj_end;

  real<lower = 0> scale_tvcl;      // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;      // Prior Scale parameter for VC
  real<lower = 0> scale_tvq;       // Prior Scale parameter for Q
  real<lower = 0> scale_tvvp;      // Prior Scale parameter for VP
  
  real<lower = 0> scale_omega_cl; // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_vc; // Prior scale parameter for omega_vc
  real<lower = 0> scale_omega_q;  // Prior scale parameter for omega_q
  real<lower = 0> scale_omega_vp; // Prior scale parameter for omega_vp
  
  real<lower = 0> lkj_df_omega;   // Prior degrees of freedom for omega cor mat
  
  real<lower = 0> scale_sigma_p;  // Prior Scale parameter for proportional error
  real<lower = 0> scale_sigma_a;  // Prior Scale parameter for additive error
  real<lower = 0> lkj_df_sigma;   // Prior degrees of freedom for sigma cor mat
 
}
transformed data{ 
  
  vector[n_obs] dv_obs = dv[i_obs]; 
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  int n_random = 4;     // Number of random effects
  int n_cmt = 3;        // Number of compartments in the model (includes depot)
  
  array[n_random] real scale_omega = {scale_omega_cl, scale_omega_vc, 
                                      scale_omega_q, scale_omega_vp}; 
                                      
  array[2] real scale_sigma = {scale_sigma_p, scale_sigma_a}; 
  
}
parameters{  
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = 0> TVQ;
  real<lower = 0> TVVP;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  vector<lower = 0>[2] sigma;
  cholesky_factor_corr[2] L_Sigma;
  
  matrix[n_random, n_subjects] Z;
  
}
transformed parameters{
  
  vector[n_subjects] eta_cl;
  vector[n_subjects] eta_vc;
  vector[n_subjects] eta_q;
  vector[n_subjects] eta_vp;
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] Q;
  vector[n_subjects] VP;

  vector[n_obs] ipred;
  
  matrix[2, 2] R_Sigma;
  matrix[2, 2] Sigma;

  {
    row_vector[n_random] typical_values = 
      to_row_vector({TVCL, TVVC, TVQ, TVVP});

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';
    
    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
                          
    R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    Sigma = quad_form_diag(R_Sigma, sigma);

    vector[n_total] dv_ipred;
    matrix[n_total, n_cmt] x_ipred;

    eta_cl = col(eta, 1);
    eta_vc = col(eta, 2);
    eta_q = col(eta, 3);
    eta_vp = col(eta, 4);

    CL = col(theta, 1);
    VC = col(theta, 2);
    Q = col(theta, 3);
    VP = col(theta, 4);

    for(j in 1:n_subjects){

      array[n_random + 1] real theta_params = {CL[j], Q[j], VC[j], VP[j], 0}; // KA = 0. Skip absorption
      
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_twocpt(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         theta_params)';
                      
                      
      dv_ipred[subj_start[j]:subj_end[j]] = 
        x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];

    }

    ipred = dv_ipred[i_obs];

  }
  
}

model{ 
  
  // Priors
  TVCL ~ cauchy(0, scale_tvcl);
  TVVC ~ cauchy(0, scale_tvvc);
  TVQ ~ cauchy(0, scale_tvq);
  TVVP ~ cauchy(0, scale_tvvp);

  omega ~ normal(0, scale_omega);
  L ~ lkj_corr_cholesky(lkj_df_omega);
  
  sigma ~ normal(0, scale_sigma);
  L_Sigma ~ lkj_corr_cholesky(lkj_df_sigma);
  
  to_vector(Z) ~ std_normal();
  
  // Likelihood
  for(i in 1:n_obs){
    dv_obs[i] ~ normal(ipred[i], sqrt(square(ipred[i]) * Sigma[1, 1] + 
                                             Sigma[2, 2] + 
                                             2*ipred[i]*Sigma[2, 1]));
  }

}
generated quantities{

  real<lower = 0> sigma_p = sigma[1];
  real<lower = 0> sigma_a = sigma[2];
  
  real<lower = 0> sigma_sq_p = square(sigma_p);
  real<lower = 0> sigma_sq_a = square(sigma_a); 

  real<lower = 0> omega_cl = omega[1];
  real<lower = 0> omega_vc = omega[2];
  real<lower = 0> omega_q = omega[3];
  real<lower = 0> omega_vp = omega[4];

  real<lower = 0> omega_sq_cl = square(omega_cl);
  real<lower = 0> omega_sq_vc = square(omega_vc);
  real<lower = 0> omega_sq_q = square(omega_q);
  real<lower = 0> omega_sq_vp = square(omega_vp);

  real cor_cl_vc;
  real cor_cl_q;
  real cor_cl_vp;
  real cor_vc_q;
  real cor_vc_vp;
  real cor_q_vp;
  
  real omega_cl_vc;
  real omega_cl_q;
  real omega_cl_vp;
  real omega_vc_q;
  real omega_vc_vp;
  real omega_q_vp;

  real cor_p_a;
  real sigma_p_a;
  
  vector[n_obs] pred;
  vector[n_obs] dv_ppc;
  vector[n_obs] log_lik;
  vector[n_obs] res;
  vector[n_obs] wres;
  vector[n_obs] ires;
  vector[n_obs] iwres;


  {

    matrix[n_random, n_random] R = multiply_lower_tri_self_transpose(L);
    matrix[n_random, n_random] Omega = quad_form_diag(R, omega);

    vector[n_total] dv_pred;
    matrix[n_total, n_cmt] x_pred;

    array[n_random + 1] real theta_params_tv = {TVCL, TVQ, TVVC, TVVP, 0}; // KA = 0. Skip absorption

    cor_cl_vc = R[1, 2];
    cor_cl_q = R[1, 3];
    cor_cl_vp = R[1, 4];
    cor_vc_q = R[2, 3];
    cor_vc_vp = R[2, 4];
    cor_q_vp = R[3, 4];
    
    omega_cl_vc = Omega[1, 2];
    omega_cl_q = Omega[1, 3];
    omega_cl_vp = Omega[1, 4];
    omega_vc_q = Omega[2, 3];
    omega_vc_vp = Omega[2, 4];
    omega_q_vp = Omega[3, 4];
    
    cor_p_a = R_Sigma[2, 1];
    sigma_p_a = Sigma[2, 1];

    for(j in 1:n_subjects){

      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_twocpt(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         theta_params_tv)';

      dv_pred[subj_start[j]:subj_end[j]] = 
        x_pred[subj_start[j]:subj_end[j], 2] ./ TVVC;
    }

    pred = dv_pred[i_obs];

  }

  res = dv_obs - pred;                                                                             
  ires = dv_obs - ipred;

  for(i in 1:n_obs){
    real ipred_tmp = ipred[i];
    real sigma_tmp = sqrt(square(ipred_tmp) * sigma_sq_p + sigma_sq_a + 
                          2*ipred_tmp*sigma_p_a);
    dv_ppc[i] = normal_rng(ipred_tmp, sigma_tmp);
    log_lik[i] = normal_lpdf(dv_obs[i] | ipred_tmp, sigma_tmp);
    wres[i] = res[i]/sigma_tmp;
    iwres[i] = ires[i]/sigma_tmp;
  }

}

