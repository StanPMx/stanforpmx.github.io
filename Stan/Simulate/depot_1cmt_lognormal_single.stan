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
  
  real<lower = 0> CL;     
  real<lower = 0> V;      
  real<lower = CL/V> KA;     
  
  real<lower = 0> sigma;  
 
}
transformed data{ 
  
  vector[n_obs] time_since_dose = to_vector(time) - time_of_first_dose;
  
}
model{ 

}
generated quantities{

  vector[n_obs] cp;
  vector[n_obs] dv;

  for(i in 1:n_obs){
    
    if(time_since_dose[i] <= 0){
      cp[i] = 0;
      dv[i] = 0;
    }else{
      cp[i] = depot_1cmt(dose, CL, V, KA, time_since_dose[i]);
      dv[i] = lognormal_rng(log(cp[i]), sigma);
    }
    
  }
}

