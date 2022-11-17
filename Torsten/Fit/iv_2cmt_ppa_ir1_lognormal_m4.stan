// IV Infusion
// Two-compartment PK Model with Indirect Response 1 PD Model
// IIV on CL, VC, Q, VP (full covariance matrix) for PK
// IIV on KIN, KOUT, IC50 (full covariance matrix) for PD. IMAX is fixed to be 1
// proportional plus additive error for PK - DV = CP(1 + eps_p_pk) + eps_a_pk
// lognormal error for PD - DV = f(.)*exp(eps_pd)
// General ODE solution using Torsten
// Implements threading for within-chain parallelization
// Deals with BLOQ values by the "CDF trick"
// Since we have a normal distribution on the error for PK, but the DV for PK 
//   must be > 0, it truncates the likelihood below at 0
// For PPC for PK, it generates values from a normal that is truncated below at 
//   0.


functions{

  array[] int sequence(int start, int end) { 
    array[end - start + 1] int seq;
    for (n in 1:num_elements(seq)) {
      seq[n] = n + start - 1;
    }
    return seq; 
  } 
  
  int num_between(int lb, int ub, array[] int y){
    
    int n = 0;
    for(i in 1:num_elements(y)){
      if(y[i] >= lb && y[i] <= ub)
                    n = n + 1;
    }
    return n;
    
  }
  
  array[] int find_between(int lb, int ub, array[] int y) {
    // vector[num_between(lb, ub, y)] result;
    array[num_between(lb, ub, y)] int result;
    int n = 1;
    for (i in 1:num_elements(y)) {
      if (y[i] >= lb && y[i] <= ub) {
        result[n] = y[i];
        n = n + 1;
      }
    }
    return result;
  }
  
  vector find_between_vec(int lb, int ub, array[] int idx, vector y) {
    
    vector[num_between(lb, ub, idx)] result;
    int n = 1;
    if(num_elements(idx) != num_elements(y)) reject("illegal input");
    for (i in 1:rows(y)) {
      if (idx[i] >= lb && idx[i] <= ub) {
        result[n] = y[i];
        n = n + 1;
      }
    }
    return result;
  }
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
  vector two_cmt_ir1_ode(real t, vector y, array[] real params, 
                         array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real q = params[2];
    real vc = params[3];
    real vp = params[4];
    real ka = params[5];
    real kin = params[6];
    real kout = params[7];
    real ic50 = params[8];
    real imax = params[9];  // It's fixed to 1 in this particular model
    real hill = params[10]; // It's fixed to 1 in this particular model
    real r_0 = params[11];
    
    real ke = cl/vc;
    real k_cp = q/vc;
    real k_pc = q/vp;
    
    real conc = y[2]/vc;
    
    real inh = (imax*pow(conc, hill))/(pow(ic50, hill) + pow(conc, hill));
    real response = y[4] + r_0;
    
    
    vector[4] dydt;

    dydt[1] = -ka*y[1];
    dydt[2] = ka*y[1] - (ke + k_cp)*y[2] + k_pc*y[3];
    dydt[3] = k_cp*y[2] - k_pc*y[3];
    dydt[4] = kin*(1 - inh) - kout*response;
    
    return dydt;
  }
  
  real partial_sum_lpmf(array[] int seq_subj, int start, int end,
                        vector dv_obs, array[] int dv_obs_id, array[] int i_obs,
                        array[] real amt, array[] int cmt, array[] int evid, 
                        array[] real time, array[] real rate, array[] real ii, 
                        array[] int addl, array[] int ss,
                        array[] int subj_start, array[] int subj_end, 
                        real TVCL, real TVVC, real TVQ, real TVVP, 
                        vector omega, matrix L, matrix Z, 
                        vector sigma, matrix L_Sigma,
                        real TVKIN, real TVKOUT, real TVIC50, real imax, real hill,
                        vector omega_pd, matrix L_pd, matrix Z_pd, 
                        real sigma_pd,
                        vector lloq, array[] int bloq,
                        int n_cmt,
                        int n_random, int n_subjects, int n_total, 
                        int n_random_pd){
                           
    real ptarget = 0;
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVQ, TVVP});
    
    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
                          
    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
    
    row_vector[n_random_pd] typical_values_pd = to_row_vector({TVKIN, TVKOUT, 
                                                               TVIC50});
    
    matrix[n_subjects, n_random_pd] eta_pd = diag_pre_multiply(omega_pd, 
                                                               L_pd * Z_pd)';

    matrix[n_subjects, n_random_pd] theta_pd =
                          (rep_matrix(typical_values_pd, n_subjects) .* 
                              exp(eta_pd));
    
                              
    int N = end - start + 1;        // number of subjects in this slice  
    vector[n_total] dv_ipred;       
    matrix[n_total, n_cmt] x_ipred; 
  

    int n_obs_slice = num_between(subj_start[start], subj_end[end], i_obs);
    array[n_obs_slice] int i_obs_slice = find_between(subj_start[start], 
                                                      subj_end[end], i_obs);
                                                
    vector[n_obs_slice] dv_obs_slice = find_between_vec(start, end, 
                                                        dv_obs_id, dv_obs);
    
    vector[n_obs_slice] ipred_slice;
    
    vector[n_obs_slice] lloq_slice = lloq[i_obs_slice];
    array[n_obs_slice] int bloq_slice = bloq[i_obs_slice];
    
    array[n_obs_slice] int cmt_slice = cmt[i_obs_slice];
    
    
    for(n in 1:N){            // loop over subjects in this slice
    
      int nn = n + start - 1; // nn is the ID of the current subject
      
      row_vector[n_random] theta_nn = theta[nn]; // access the PK parameters for subject nn
      real cl = theta_nn[1];
      real vc = theta_nn[2];
      real q = theta_nn[3];
      real vp = theta_nn[4];
      
      row_vector[n_random_pd] theta_nn_pd = theta_pd[nn]; // access the PD parameters for subject nn
      real kin = theta_nn_pd[1];
      real kout = theta_nn_pd[2];
      real ic50 = theta_nn_pd[3];
      real r_0 = kin/kout;
      
      array[11] real params = {cl, q, vc, vp, 0,  // The 0 is for KA. Skip the absorption
                               kin, kout, ic50, imax, hill, r_0};
      
      x_ipred[subj_start[nn]:subj_end[nn], ] =
        pmx_solve_rk45(two_cmt_ir1_ode,
                       n_cmt,
                       time[subj_start[nn]:subj_end[nn]],
                       amt[subj_start[nn]:subj_end[nn]],
                       rate[subj_start[nn]:subj_end[nn]],
                       ii[subj_start[nn]:subj_end[nn]],
                       evid[subj_start[nn]:subj_end[nn]],
                       cmt[subj_start[nn]:subj_end[nn]],
                       addl[subj_start[nn]:subj_end[nn]],
                       ss[subj_start[nn]:subj_end[nn]],
                       params)';
                      
      for(k in subj_start[nn]:subj_end[nn]){
        if(cmt[k] == 2){
          dv_ipred[k] = x_ipred[k, 2] / vc;
        }else if(cmt[k] == 4){
          dv_ipred[k] = x_ipred[k, 4] + r_0;
        }
      }       
    
    }
  
    ipred_slice = dv_ipred[i_obs_slice];
    
    for(i in 1:n_obs_slice){
      if(cmt_slice[i] == 2){
        real sigma_tmp = sqrt(square(ipred_slice[i]) * Sigma[1, 1] +
                                           Sigma[2, 2] + 
                                           2*ipred_slice[i]*Sigma[2, 1]);
        if(bloq_slice[i] == 1){
          ptarget += log_diff_exp(normal_lcdf(lloq_slice[i] | ipred_slice[i], 
                                                              sigma_tmp),
                                  normal_lcdf(0.0 | ipred_slice[i], sigma_tmp)) -
                     normal_lccdf(0.0 | ipred_slice[i], sigma_tmp); 
        }else{
          ptarget += normal_lpdf(dv_obs_slice[i] | ipred_slice[i], sigma_tmp) -
                     normal_lccdf(0.0 | ipred_slice[i], sigma_tmp);
        }
      }else if(cmt_slice[i] == 4){
        if(bloq_slice[i] == 1){
          ptarget += lognormal_lcdf(lloq_slice[i] | log(ipred_slice[i]), sigma_pd);
        }else{
          ptarget += lognormal_lpdf(dv_obs_slice[i] | log(ipred_slice[i]), sigma_pd);
        }
      }                                         
    }                          
    return ptarget;
                           
  }
  
}
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
  vector<lower = 0>[n_total] dv;
  array[n_subjects] int subj_start;
  array[n_subjects] int subj_end;
  vector[n_total] lloq;
  array[n_total] int bloq;
  
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
 
  real<lower = 0> scale_tvkin;    // Prior Scale parameter for CL
  real<lower = 0> scale_tvkout;   // Prior Scale parameter for VC
  real<lower = 0> scale_tvic50;   // Prior Scale parameter for Q
  
  real<lower = 0> scale_omega_kin;  // Prior scale parameter for omega_kin
  real<lower = 0> scale_omega_kout; // Prior scale parameter for omega_kout
  real<lower = 0> scale_omega_ic50; // Prior scale parameter for omega_ic50
  
  real<lower = 0> lkj_df_omega_pd;  // Prior degrees of freedom for omega cor mat for PD
  
  real<lower = 0> scale_sigma_pd;   // Prior Scale parameter for lognormal error for PD
 
}
transformed data{ 
  
  int grainsize = 1;
  
  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  vector[n_obs] lloq_obs = lloq[i_obs];
  array[n_obs] int bloq_obs = bloq[i_obs];
  
  int n_random = 4;                    // Number of random effects for PK
  int n_random_pd = 3;                 // Number of random effects for PD
  
  array[n_random] real scale_omega = {scale_omega_cl, scale_omega_vc, 
                                      scale_omega_q, scale_omega_vp}; 
  array[2] real scale_sigma = {scale_sigma_p, scale_sigma_a};
  
  array[n_random_pd] real scale_omega_pd = {scale_omega_kin, scale_omega_kout, 
                                            scale_omega_ic50}; 
  
  array[n_subjects] int seq_subj = sequence(1, n_subjects); // reduce_sum over subjects
  
  int n_cmt = 4;
  
  real imax = 1.0;
  real hill = 1.0;
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
  
  real<lower = 0> TVKIN;       
  real<lower = 0> TVKOUT; 
  real<lower = 0> TVIC50;
  
  vector<lower = 0>[n_random_pd] omega_pd;
  cholesky_factor_corr[n_random_pd] L_pd;
  
  real<lower = 0> sigma_pd;
  
  matrix[n_random_pd, n_subjects] Z_pd;
  
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
  
  TVKIN ~ cauchy(0, scale_tvkin);
  TVKOUT ~ cauchy(0, scale_tvkout);
  TVIC50 ~ cauchy(0, scale_tvic50);
  
  omega_pd ~ normal(0, scale_omega_pd);
  L_pd ~ lkj_corr_cholesky(lkj_df_omega_pd);
  
  sigma_pd ~ normal(0, scale_sigma_pd);
  
  to_vector(Z_pd) ~ std_normal();
  
  // Likelihood
  target += reduce_sum(partial_sum_lupmf, seq_subj, grainsize,
                       dv_obs, dv_obs_id, i_obs,
                       amt, cmt, evid, time, 
                       rate, ii, addl, ss, subj_start, subj_end, 
                       TVCL, TVVC, TVQ, TVVP, omega, L, Z,
                       sigma, L_Sigma, 
                       TVKIN, TVKOUT, TVIC50, imax, hill,
                       omega_pd, L_pd, Z_pd, 
                       sigma_pd,
                       lloq, bloq, n_cmt,
                       n_random, n_subjects, n_total, n_random_pd);
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

  real<lower = 0> omega_kin = omega_pd[1];
  real<lower = 0> omega_kout = omega_pd[2];
  real<lower = 0> omega_ic50 = omega_pd[3];

  real<lower = 0> omega_sq_kin = square(omega_kin);
  real<lower = 0> omega_sq_kout = square(omega_kout);
  real<lower = 0> omega_sq_ic50 = square(omega_ic50);

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

  vector[n_subjects] eta_cl;
  vector[n_subjects] eta_vc;
  vector[n_subjects] eta_q;
  vector[n_subjects] eta_vp;
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] Q;
  vector[n_subjects] VP;
  
  real cor_kin_kout;
  real cor_kin_ic50;
  real cor_kout_ic50;
  real omega_kin_kout;
  real omega_kin_ic50;
  real omega_kout_ic50;
  
  vector[n_subjects] eta_kin;
  vector[n_subjects] eta_kout;
  vector[n_subjects] eta_ic50;
  vector[n_subjects] KIN;
  vector[n_subjects] KOUT;
  vector[n_subjects] IC50;
  
  vector[n_obs] ipred;
  vector[n_obs] pred;
  vector[n_obs] dv_ppc;
  vector[n_obs] log_lik;
  vector[n_obs] res;
  vector[n_obs] wres;
  vector[n_obs] ires;
  vector[n_obs] iwres;
 
  {
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVQ, TVVP});

    matrix[n_random, n_random] R = multiply_lower_tri_self_transpose(L);
    matrix[n_random, n_random] Omega = quad_form_diag(R, omega);

    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
                          
    row_vector[n_random_pd] typical_values_pd = to_row_vector({TVKIN, TVKOUT, 
                                                               TVIC50});
    
    matrix[n_random_pd, n_random_pd] R_pd = 
                                        multiply_lower_tri_self_transpose(L_pd);
    matrix[n_random_pd, n_random_pd] Omega_pd = quad_form_diag(R_pd, omega_pd);
    
    matrix[n_subjects, n_random_pd] eta_pd = diag_pre_multiply(omega_pd, 
                                                               L_pd * Z_pd)';

    matrix[n_subjects, n_random_pd] theta_pd =
                          (rep_matrix(typical_values_pd, n_subjects) .* 
                              exp(eta_pd));

    vector[n_total] dv_pred;
    matrix[n_total, n_cmt] x_pred;
    vector[n_total] dv_ipred;
    matrix[n_total, n_cmt] x_ipred;
    
    real r_0_tv = TVKIN/TVKOUT;
    
    array[11] real params_tv = {TVCL, TVQ, TVVC, TVVP, 0,  // 0 is for KA to skip absorption
                                TVKIN, TVKOUT, TVIC50, imax, hill,
                                r_0_tv};
    
    eta_cl = col(eta, 1);
    eta_vc = col(eta, 2);
    eta_q = col(eta, 3);
    eta_vp = col(eta, 4);
    
    CL = col(theta, 1);
    VC = col(theta, 2);
    Q = col(theta, 3);
    VP = col(theta, 4);

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
    
    eta_kin = col(eta_pd, 1);
    eta_kout = col(eta_pd, 2);
    eta_ic50 = col(eta_pd, 3);
    
    KIN = col(theta_pd, 1);
    KOUT = col(theta_pd, 2);
    IC50 = col(theta_pd, 3);

    cor_kin_kout = R_pd[1, 2];
    cor_kin_ic50 = R_pd[1, 3];
    cor_kout_ic50 = R_pd[2, 3];
    
    omega_kin_kout = Omega_pd[1, 2];
    omega_kin_ic50 = Omega_pd[1, 3];
    omega_kout_ic50 = Omega_pd[2, 3];
    
    
    for(j in 1:n_subjects){
      
      real r_0 = KIN[j]/KOUT[j];
      array[11] real params = {CL[j], Q[j], VC[j], VP[j], 0,
                               KIN[j], KOUT[j], IC50[j], imax, hill,
                               r_0};
      
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(two_cmt_ir1_ode,
                       n_cmt,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       params)';
                      
      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){
          dv_ipred[k] = x_ipred[k, 2] / VC[j];
        }else if(cmt[k] == 4){
          dv_ipred[k] = x_ipred[k, 4] + r_0;
        }
      }    


      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(two_cmt_ir1_ode, 
                       n_cmt, 
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       params_tv)';
      
      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){
          dv_pred[k] = x_pred[k, 2] / TVVC;
        }else if(cmt[k] == 4){
          dv_pred[k] = x_pred[k, 4] + r_0_tv;
        }
      }  
    }

    pred = dv_pred[i_obs];
    ipred = dv_ipred[i_obs];

  }
  
  for(i in 1:n_obs){
    if(cmt[i_obs[i]] == 2){
      real ipred_tmp = ipred[i];
      real sigma_tmp = sqrt(square(ipred_tmp) * sigma_sq_p + sigma_sq_a + 
                            2*ipred_tmp*sigma_p_a);
      dv_ppc[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
      if(bloq_obs[i] == 1){
        log_lik[i] = log_diff_exp(normal_lcdf(lloq_obs[i] | ipred_tmp, sigma_tmp),
                                  normal_lcdf(0.0 | ipred_tmp, sigma_tmp)) -
                     normal_lccdf(0.0 | ipred_tmp, sigma_tmp); 
      }else{
        log_lik[i] = normal_lpdf(dv_obs[i] | ipred_tmp, sigma_tmp) - 
                     normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
      }
      res[i] = dv_obs[i] - pred[i];                                                                             
      ires[i] = dv_obs[i] - ipred[i];
      wres[i] = res[i]/sigma_tmp;
      iwres[i] = ires[i]/sigma_tmp;
    }else if(cmt[i_obs[i]] == 4){
      dv_ppc[i] = lognormal_rng(log(ipred[i]), sigma_pd);
      if(bloq_obs[i] == 1){
        log_lik[i] = lognormal_lcdf(lloq_obs[i] | log(ipred[i]), sigma_pd);
      }else{
        log_lik[i] = lognormal_lpdf(dv_obs[i] | log(ipred[i]), sigma_pd);
      }
      res[i] = log(dv_obs[i]) - log(pred[i]);
      ires[i] = log(dv_obs[i]) - log(ipred[i]);
      wres[i] = res[i]/sigma_pd;
      iwres[i] = ires[i]/sigma_pd;
    }
  }
  
}

