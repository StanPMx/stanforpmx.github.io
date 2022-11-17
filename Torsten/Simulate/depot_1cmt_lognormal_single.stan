// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model
// Single subject
// lognormal error - DV = CP*exp(eps)
// Closed form solution using a Torsten function

data{
  
  int n_obs;
  array[n_obs] real amt;
  array[n_obs] int cmt;
  array[n_obs] int evid;
  array[n_obs] real rate;
  array[n_obs] real ii;
  array[n_obs] int addl;
  array[n_obs] int ss;
  array[n_obs] real time;
  
  real<lower = 0> CL;     
  real<lower = 0> V;      
  real<lower = CL/V> KA;     
  
  real<lower = 0> sigma;  
 
}
model{ 

}
generated quantities{

  vector[n_obs] cp;
  vector[n_obs] dv;

  {
    
    matrix[n_obs, 2] x_cp;
    array[3] real theta_params = {CL, V, KA}; 
    
    x_cp = pmx_solve_onecpt(time,
                            amt,
                            rate,
                            ii,
                            evid,
                            cmt,
                            addl,
                            ss,
                            theta_params)';
                               
    cp = x_cp[, 2] ./ V;
    
  }

  for(i in 1:n_obs){
    if(cp[i] == 0){
      dv[i] = 0;
    }else{
      dv[i] = lognormal_rng(log(cp[i]), sigma);  
    }
    
  }
}


     