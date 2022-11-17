// IV Infusion
// Two-compartment PK Model
// IIV on CL, VC, Q, VP (full covariance matrix)
// proportional plus additive error - DV = CP(1 + eps_p) + eps_a
// Closed form solution using Torsten
// Predictions are generated from a normal that is truncated below at 0

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
}

data{
  
  int n_subjects;
  int n_subjects_new;
  int n_time_new;
  array[n_time_new] real time;
  array[n_time_new] real amt;
  array[n_time_new] int cmt;
  array[n_time_new] int evid;
  array[n_time_new] real rate;
  array[n_time_new] real ii;
  array[n_time_new] int addl;
  array[n_time_new] int ss;
  array[n_subjects_new] int subj_start;
  array[n_subjects_new] int subj_end;
  
}
transformed data{ // done
  
  int n_random = 4;                    // Number of random effects
  
}
parameters{  // done
  
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
generated quantities{

  vector[n_time_new] cp; // concentration with no residual error
  vector[n_time_new] dv; // concentration with residual error

  {
    vector[n_random] typical_values = to_vector({TVCL, TVVC, TVQ, TVVP});
    
    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
    
    matrix[n_random, n_subjects_new] eta_new;
    matrix[n_subjects_new, n_random] theta_new;
    vector[n_subjects_new] CL;
    vector[n_subjects_new] VC;
    vector[n_subjects_new] Q;
    vector[n_subjects_new] VP;
    matrix[n_time_new, 3] x_cp;
    
    for(i in 1:n_subjects_new){
      eta_new[, i] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L));
    }
    theta_new = (rep_matrix(typical_values, n_subjects_new) .* exp(eta_new))';
    
    CL = col(theta_new, 1);
    VC = col(theta_new, 2);
    Q = col(theta_new, 3);
    VP = col(theta_new, 4);
    
    for(j in 1:n_subjects_new){
      
      array[5] real theta_params = {CL[j], Q[j], VC[j], VP[j], 0};
      
      x_cp[subj_start[j]:subj_end[j],] =
        pmx_solve_twocpt(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         theta_params)';

      cp[subj_start[j]:subj_end[j]] = x_cp[subj_start[j]:subj_end[j], 2] ./ VC[j];
    
    }


    for(i in 1:n_time_new){
      if(cp[i] == 0){
        dv[i] = 0;
      }else{
        real cp_tmp = cp[i];
        real sigma_tmp = sqrt(square(cp_tmp) * Sigma[1, 1] + Sigma[2, 2] + 
                              2*cp_tmp*Sigma[2, 1]);
        dv[i] = normal_lb_rng(cp_tmp, sigma_tmp, 0.0);
        
      }
    }
  
  }

}


