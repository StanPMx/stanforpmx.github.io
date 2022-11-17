// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model
// IIV on CL. V. and Ka (full covariance matrix)
// proportional plus additive error - DV = CP(1 + eps_p) + eps_a
// Linear ODE solution using Torsten
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
                        real TVCL, real TVV, real TVKA, 
                        vector omega, matrix L, matrix Z, 
                        vector sigma, matrix L_Sigma,
                        vector lloq, array[] int bloq, 
                        int n_cmt, array[] real bioav, array[] real tlag, 
                        int n_random, int n_subjects, int n_total){
                           
    real ptarget = 0;
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVV, TVKA});
    
    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
                          
    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
    
                              
    int N = end - start + 1;    // number of subjects in this slice  
    vector[n_total] dv_ipred;   
    matrix[n_total, 2] x_ipred; 
  
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
      real v = theta_nn[2];
      real ka = theta_nn[3];
      
      matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
      K[1, 1] = -ka;
      K[2, 1] = ka;
      K[2, 2] = -cl/v;
      
      x_ipred[subj_start[nn]:subj_end[nn], ] =
        pmx_solve_linode(time[subj_start[nn]:subj_end[nn]],
                         amt[subj_start[nn]:subj_end[nn]],
                         rate[subj_start[nn]:subj_end[nn]],
                         ii[subj_start[nn]:subj_end[nn]],
                         evid[subj_start[nn]:subj_end[nn]],
                         cmt[subj_start[nn]:subj_end[nn]],
                         addl[subj_start[nn]:subj_end[nn]],
                         ss[subj_start[nn]:subj_end[nn]],
                         K, bioav, tlag)';
                      
      dv_ipred[subj_start[nn]:subj_end[nn]] = 
        x_ipred[subj_start[nn]:subj_end[nn], 2] ./ v;
    
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
  
  real<lower = 0> scale_tvcl;     // Prior Scale parameter for CL
  real<lower = 0> scale_tvv;      // Prior Scale parameter for V
  real<lower = 0> scale_tvka;     // Prior Scale parameter for KA
  
  real<lower = 0> scale_omega_cl; // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_v;  // Prior scale parameter for omega_v
  real<lower = 0> scale_omega_ka; // Prior scale parameter for omega_ka
  
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
  
  int n_random = 3;                    // Number of random effects
  
  array[n_random] real scale_omega = {scale_omega_cl, scale_omega_v, 
                                      scale_omega_ka}; 
  array[2] real scale_sigma = {scale_sigma_p, scale_sigma_a};
  
  array[n_subjects] int seq_subj = sequence(1, n_subjects); // reduce_sum over subjects
  
  int n_cmt = 2;
  array[n_cmt] real bioav = rep_array(1.0, n_cmt);
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVV; 
  real<lower = TVCL/TVV> TVKA;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  vector<lower = 0>[2] sigma;
  cholesky_factor_corr[2] L_Sigma;
  
  matrix[n_random, n_subjects] Z;
  
}

model{ 
  
  // Priors
  TVCL ~ cauchy(0, scale_tvcl);
  TVV ~ cauchy(0, scale_tvv);
  TVKA ~ normal(0, scale_tvka) T[TVCL/TVV, ];

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
                       TVCL, TVV, TVKA, omega, L, Z,
                       sigma, L_Sigma, 
                       lloq, bloq,
                       n_cmt, bioav, tlag,
                       n_random, n_subjects, n_total);
}
generated quantities{
  
  real<lower = 0> sigma_p = sigma[1];
  real<lower = 0> sigma_a = sigma[2];
  
  real<lower = 0> sigma_sq_p = square(sigma_p);
  real<lower = 0> sigma_sq_a = square(sigma_a); 

  real<lower = 0> omega_cl = omega[1];
  real<lower = 0> omega_v = omega[2];
  real<lower = 0> omega_ka = omega[3];

  real<lower = 0> omega_sq_cl = square(omega_cl);
  real<lower = 0> omega_sq_v = square(omega_v);
  real<lower = 0> omega_sq_ka = square(omega_ka);

  real cor_cl_v;
  real cor_cl_ka;
  real cor_v_ka;
  real omega_cl_v;
  real omega_cl_ka;
  real omega_v_ka;
  
  real cor_p_a;
  real sigma_p_a;

  vector[n_subjects] eta_cl;
  vector[n_subjects] eta_v;
  vector[n_subjects] eta_ka;
  vector[n_subjects] CL;
  vector[n_subjects] V;
  vector[n_subjects] KA;
  vector[n_subjects] KE;

  vector[n_obs] ipred;
  vector[n_obs] pred;
  vector[n_obs] dv_ppc;
  vector[n_obs] log_lik;
  vector[n_obs] res;
  vector[n_obs] wres;
  vector[n_obs] ires;
  vector[n_obs] iwres;
 
  {
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVV, TVKA});

    matrix[n_random, n_random] R = multiply_lower_tri_self_transpose(L);
    matrix[n_random, n_random] Omega = quad_form_diag(R, omega);

    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));

    vector[n_total] dv_pred;
    matrix[n_total, 2] x_pred;
    vector[n_total] dv_ipred;
    matrix[n_total, 2] x_ipred;

    matrix[n_cmt, n_cmt] K_tv = rep_matrix(0, n_cmt, n_cmt);
    K_tv[1, 1] = -TVKA;
    K_tv[2, 1] = TVKA;
    K_tv[2, 2] = -TVCL/TVV;

    eta_cl = col(eta, 1);
    eta_v = col(eta, 2);
    eta_ka = col(eta, 3);

    CL = col(theta, 1);
    V = col(theta, 2);
    KA = col(theta, 3);
    KE = CL ./ V;

    cor_cl_v = R[1, 2];
    cor_cl_ka = R[1, 3];
    cor_v_ka = R[2, 3];

    omega_cl_v = Omega[1, 2];
    omega_cl_ka = Omega[1, 3];
    omega_v_ka = Omega[2, 3];
    
    cor_p_a = R_Sigma[2, 1];
    sigma_p_a = Sigma[2, 1];


    for(j in 1:n_subjects){

      matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
      K[1, 1] = -KA[j];
      K[2, 1] = KA[j];
      K[2, 2] = -CL[j]/V[j];
      
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         K, bioav, tlag)';
                      
                      
      dv_ipred[subj_start[j]:subj_end[j]] = 
        x_ipred[subj_start[j]:subj_end[j], 2] ./ V[j];

      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         K_tv, bioav, tlag)';

      dv_pred[subj_start[j]:subj_end[j]] = 
        x_pred[subj_start[j]:subj_end[j], 2] ./ TVV;
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

