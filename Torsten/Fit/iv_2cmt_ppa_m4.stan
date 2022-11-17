// IV Infusion
// Two-compartment PK Model
// IIV on CL, VC, Q, VP (full covariance matrix)
// proportional plus additive error - DV = CP(1 + eps_p) + eps_a
// Closed form solution using Torsten
// Implements threading for within-chain parallelization
// Deals with BLOQ values by the "CDF trick" (M4)
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   truncates the likelihood below at 0
// For PPC, it generates values from a normal that is truncated below at 0


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
  
  real partial_sum_lpmf(array[] int seq_subj, int start, int end,
                        vector dv_obs, array[] int dv_obs_id, array[] int i_obs,
                        array[] real amt, array[] int cmt, array[] int evid, 
                        array[] real time, array[] real rate, array[] real ii, 
                        array[] int addl, array[] int ss,
                        array[] int subj_start, array[] int subj_end, 
                        real TVCL, real TVVC, real TVQ, real TVVP, 
                        vector omega, matrix L, matrix Z, 
                        vector sigma, matrix L_Sigma,
                        vector lloq, array[] int bloq,
                        int n_random, int n_subjects, int n_total){
                           
    real ptarget = 0;
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVQ, TVVP});
    
    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
                          
    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
    
                              
    int N = end - start + 1;    // number of subjects in this slice  
    vector[n_total] dv_ipred;   
    matrix[n_total, 3] x_ipred; 
  

    int n_obs_slice = num_between(subj_start[start], subj_end[end], i_obs);
    array[n_obs_slice] int i_obs_slice = find_between(subj_start[start], 
                                                      subj_end[end], i_obs);
                                                
    vector[n_obs_slice] dv_obs_slice = find_between_vec(start, end, 
                                                        dv_obs_id, dv_obs);
    
    vector[n_obs_slice] ipred_slice;
    
    vector[n_obs_slice] lloq_slice = lloq[i_obs_slice];
    array[n_obs_slice] int bloq_slice = bloq[i_obs_slice];
    
    
    for(n in 1:N){            // loop over subjects in this slice
    
      int nn = n + start - 1; // nn is the ID of the current subject
      
      row_vector[n_random] theta_nn = theta[nn]; // access the parameters for subject nn
      real cl = theta_nn[1];
      real vc = theta_nn[2];
      real q = theta_nn[3];
      real vp = theta_nn[4];
      array[5] real theta_params = {cl, q, vc, vp, 0}; // The 0 is for KA. Skip the absorption
      
      x_ipred[subj_start[nn]:subj_end[nn], ] =
        pmx_solve_twocpt(time[subj_start[nn]:subj_end[nn]],
                         amt[subj_start[nn]:subj_end[nn]],
                         rate[subj_start[nn]:subj_end[nn]],
                         ii[subj_start[nn]:subj_end[nn]],
                         evid[subj_start[nn]:subj_end[nn]],
                         cmt[subj_start[nn]:subj_end[nn]],
                         addl[subj_start[nn]:subj_end[nn]],
                         ss[subj_start[nn]:subj_end[nn]],
                         theta_params)';
                      
      dv_ipred[subj_start[nn]:subj_end[nn]] = 
        x_ipred[subj_start[nn]:subj_end[nn], 2] ./ vc;
    
    }
  
    
    ipred_slice = dv_ipred[i_obs_slice];
    
    for(i in 1:n_obs_slice){
      real sigma_tmp = sqrt(square(ipred_slice[i]) * Sigma[1, 1] +
                                         Sigma[2, 2] + 
                                         2*ipred_slice[i]*Sigma[2, 1]);
      if(bloq_slice[i] == 1){
        // ptarget += log(normal_cdf(lloq_slice[i] | ipred_slice[i], sigma_tmp) -
        //                normal_cdf(0.0 | ipred_slice[i], sigma_tmp)) -
        //            normal_lccdf(0.0 | ipred_slice[i], sigma_tmp);  
        ptarget += log_diff_exp(normal_lcdf(lloq_slice[i] | ipred_slice[i], 
                                                            sigma_tmp),
                                normal_lcdf(0.0 | ipred_slice[i], sigma_tmp)) -
                   normal_lccdf(0.0 | ipred_slice[i], sigma_tmp); 
      }else{
        ptarget += normal_lpdf(dv_obs_slice[i] | ipred_slice[i], sigma_tmp) -
                   normal_lccdf(0.0 | ipred_slice[i], sigma_tmp);
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
 
}
transformed data{ 
  
  int grainsize = 1;
  
  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  vector[n_obs] lloq_obs = lloq[i_obs];
  array[n_obs] int bloq_obs = bloq[i_obs];
  
  int n_random = 4;                    // Number of random effects
  
  array[n_random] real scale_omega = {scale_omega_cl, scale_omega_vc, 
                                      scale_omega_q, scale_omega_vp}; 
  array[2] real scale_sigma = {scale_sigma_p, scale_sigma_a};
  
  array[n_subjects] int seq_subj = sequence(1, n_subjects); // reduce_sum over subjects
  
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
  target += reduce_sum(partial_sum_lupmf, seq_subj, grainsize,
                       dv_obs, dv_obs_id, i_obs,
                       amt, cmt, evid, time, 
                       rate, ii, addl, ss, subj_start, subj_end, 
                       TVCL, TVVC, TVQ, TVVP, omega, L, Z,
                       sigma, L_Sigma, 
                       lloq, bloq,
                       n_random, n_subjects, n_total);
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

  vector[n_subjects] eta_cl;
  vector[n_subjects] eta_vc;
  vector[n_subjects] eta_q;
  vector[n_subjects] eta_vp;
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] Q;
  vector[n_subjects] VP;
  
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

    vector[n_total] dv_pred;
    matrix[n_total, 3] x_pred;
    vector[n_total] dv_ipred;
    matrix[n_total, 3] x_ipred;
    
    array[5] real theta_params_tv = {TVCL, TVQ, TVVC, TVVP, 0}; // 0 is for KA to skip absorption
    
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

    for(j in 1:n_subjects){

      array[5] real theta_params = {CL[j], Q[j], VC[j], VP[j], 0};
      
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
    ipred = dv_ipred[i_obs];

  }

  res = dv_obs - pred;                                                                             
  ires = dv_obs - ipred;

  for(i in 1:n_obs){
    real ipred_tmp = ipred[i];
    real sigma_tmp = sqrt(square(ipred_tmp) * sigma_sq_p + sigma_sq_a + 
                          2*ipred_tmp*sigma_p_a);
    // dv_ppc[i] = normal_rng(ipred_tmp, sigma_tmp);
    dv_ppc[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
    if(bloq_obs[i] == 1){
      // log_lik[i] = log(normal_cdf(lloq_obs[i] | ipred_tmp, sigma_tmp) - 
      //                  normal_cdf(0.0 | ipred_tmp, sigma_tmp)) -
      //              normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
      log_lik[i] = log_diff_exp(normal_lcdf(lloq_obs[i] | ipred_tmp, sigma_tmp),
                                normal_lcdf(0.0 | ipred_tmp, sigma_tmp)) -
                   normal_lccdf(0.0 | ipred_tmp, sigma_tmp); 
    }else{
      log_lik[i] = normal_lpdf(dv_obs[i] | ipred_tmp, sigma_tmp) - 
                   normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
    }
    wres[i] = res[i]/sigma_tmp;
    iwres[i] = ires[i]/sigma_tmp;
  }
  
}

