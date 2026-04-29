// sample from posterior of tiv model with n cell types with differing betas.
// includes lognormally distributed observation noise and observation threshold
// on viral load, and sampling error for cell numbers (no misclassification error)
// fit beta only.

// ODE model
functions {
  vector two_beta_innate_adaptive(real t,
  vector y,
  vector beta1,
  real delta,
  real p1,
  real c1,
  real T_0,
  int n_C, 
  real q, 
  real phi, 
  real rho, 
  real d, 
  real m, 
  real k1, 
  real A_max) {
    vector [(3*n_C+3)]dydt; // LHS of ODEs
    int V1 = 3*n_C + 1; 
    int F1 = 3*n_C + 2; 
    int A1 = 3*n_C + 3; 
    
    // because we divide all compartments by T_0, we need to transform beta, phi and k1
    // model equations
    dydt[1:n_C] = -beta1 .* y[1:n_C] * y[V1] * T_0 - phi * y[1:n_C] * y[F1] * T_0 + rho * y[(n_C+1):(2*n_C)]; // dT_i/dt
    dydt[(n_C+1):(2*n_C)] = phi * y[1:n_C] * y[F1] * T_0 - rho * y[(n_C+1):(2*n_C)]; // dR_i/dt
    dydt[(2*n_C+1):(3*n_C)] = beta1 .* y[1:n_C] * y[V1] * T_0 - delta * y[(2*n_C+1):(3*n_C)]; // dI_i/dt
    dydt[V1] = p1 * sum(y[(2*n_C+1):(3*n_C)]) - (c1 + k1 * T_0 * y[A1])*y[V1];
    dydt[F1] =q* sum(y[(2*n_C+1):(3*n_C)]) - d * y[F1];
    dydt[A1] = m * y[A1] * (1 - y[A1] / A_max);
    return dydt;
  }
  
  // calculates solution
  array[] vector calc_sol(int n_C,
  data real c1,
  data real delta,
  data real T_0,
  data real p,
  data real V_0,
  data real alpha,
  data real A_0, 
  data real q, 
  data real phi, 
  data real rho, 
  data real d, 
  data real m, 
  data real k1, 
  data real A_max, 
  data int T_S, // number of sampling times excluding 0
  data array[] real ts, // unique solving times
  vector beta1) {
    array[T_S] vector[(3*n_C+3)] sol; // solution of ODEs which corresponds to viral load and cell type numbers
    vector[3*n_C+3] y0; // vector of initial conditions for ODEs
    
    // set initial values for ODEs
    // note that we divide all compartments by T_0 to improve numerical stability
    // we will multiply all compartments by T_0 again at the end
    
    y0[1] = alpha;
    y0[2] = 1 - alpha;
    for(i in 1:(2*n_C)) {
      y0[n_C+i] = 0;
    }
    y0[3*n_C+1] = V_0 / T_0;
    y0[3*n_C+2] = 0;
    y0[3*n_C+3] = A_0 / T_0;
    
    // solve ODEs
    sol = ode_bdf(two_beta_innate_adaptive, y0, 0, ts, beta1, delta, p, c1, T_0, n_C, q, phi, rho, d, m, k1, A_max);
    for(i in 1:T_S) {
      sol[i] = sol[i] * T_0;
    }
    
    return sol;
  }
  
}

// data and fixed parameter values
data {
  int<lower=0> N_samples;
  real<lower=0.> c1; // rate of loss of infectivity of infectious virus
  real<lower=0.> delta; // decay rate of infected cells
  real<lower=0.> T_0; // initial number of cells in each well
  real<lower=0.> V_0; // initial viral load
  real<lower=0.> A_0;
  real<lower=0.> p; // production rate of virus from infected cells
  real<lower=0.> alpha; // initial proportion of cells in cell type 1
  real<lower=0.> q;
  real<lower=0.> phi;
  real<lower=0.> rho;
  real<lower=0.> d;
  real<lower=0.> m;
  real<lower=0.> k1;
  real<lower=0.> A_max;
  
  int<lower=1> T_S; // number of sampling times excluding 0
  array[T_S] real ts; // unique solving times
  
  array[N_samples, 2] real beta1;
}


// relates sampled parameters to model prediction of viral load
transformed data {
  
  array[T_S] vector[9] sol;
  
  array[N_samples, T_S] real sol_V1;
  
  // for each sample
  for(j in 1:N_samples) {

    // calculate the solution
    sol = calc_sol(2, c1, delta, T_0, p, V_0, alpha, A_0, q, phi, rho, d, m, k1, A_max, T_S, ts, to_vector(beta1[j]));
    
    // save the viral load
    for(i in 1:(T_S)) {
      sol_V1[j,i] = sol[i,7]; 
    }
  }
}

generated quantities {
  // save the transformed data
  array[N_samples, T_S] real sol_V = sol_V1;
}
