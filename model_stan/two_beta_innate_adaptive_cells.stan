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
  
  // calculates log likelihood
  real calc_log_lik(data int n_C,
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
  data real sigma_V,
  data real lod_V,
  data real rho0,
  data int T_V, // total number of samples for viral load
  data int T_C,
  data int T_S, // number of sampling times excluding 0
  data array[] int t_V_idx, // indices linking sample number to sampling time for viral load
  data array[] int t_C_idx,
  data vector obs_V, // vector of viral load data
  data array[,] int obs_C,
  data array[] real ts, // unique solving times
  vector beta1,
  int evaluate_likelihood) {
    //     log_lik = calc_log_lik(2, c1, delta, T_0, p, V_0, alpha, A_0, q, phi, rho, d, m, k1, A_max, sigma_V,
    // lod_V, T_V, T_S, t_V_idx, obs_V, ts, beta1, evaluate_likelihood);
    real log_lik = 0;
    
    if(evaluate_likelihood) {
      array[T_S] vector[(3*n_C+3)] sol; // solution of ODEs which corresponds to viral load and cell type numbers
      vector[3*n_C+3] y0; // vector of initial conditions for ODEs
      real alpha0 = (1 - rho0)/rho0;
      vector[2] temp1;
      
      temp1[2] = 0;
      
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
        for(j in 1:(2*n_C+1)) { // ensure positivity of solutions to construct simplex
        temp1[1] = sol[i,j];
        sol[i,j] = max(temp1);
        }
      }
      
      for(i in 1:T_V) {
        vector[2] y_pos; // temporary variable for rounding of viral load
        real LL;
        real sigma_V_nat_log; // sigma_V on natural log scale
        //(as Stan provides lognormal function on natural log scale)
        
        y_pos[2] = 1e-8; // temporary variable for rounding of viral load
        sigma_V_nat_log = sigma_V * 10 / exp(1);
        if(t_V_idx[i] == 0) {
          y_pos[1] = V_0;
        } else {
          y_pos[1] = sol[t_V_idx[i],(3*n_C+1)];
        }
        
        if(obs_V[i] > lod_V) {
          LL = lognormal_lpdf(obs_V[i] | log(max(y_pos)), sigma_V_nat_log);
        } else {
          LL = lognormal_lcdf(lod_V | log(max(y_pos)), sigma_V_nat_log);
        }
        log_lik += LL;
      }
      
      for(i in 1:T_C) {
        vector[2*n_C] true_n_cells_vec;
        vector[2*n_C] alpha_vec;
        real n_cells;
        // get proportion of cells in each compartment
        // combine T and R
        for(j in 1:n_C) {
          if(t_C_idx[i] == 0) {
            true_n_cells_vec[j] = y0[j] + y0[n_C + j];
            true_n_cells_vec[n_C + j] = y0[2*n_C + j];
          } else {
            true_n_cells_vec[j] = sol[t_C_idx[i],j] + sol[t_C_idx[i],n_C + j];
            true_n_cells_vec[n_C + j] = sol[t_C_idx[i],2*n_C + j];
          }
        }
        
        n_cells = sum(true_n_cells_vec);
        // and multiply by alpha0 which accounts for the overdispersion
        alpha_vec = true_n_cells_vec / n_cells * alpha0;
        for(j in 1:(2*n_C)) {
          if(alpha_vec[j] < 1e-8)
          alpha_vec[j] = 1e-8;
        }
        
        log_lik += dirichlet_multinomial_lpmf(obs_C[,i] | alpha_vec);
      }
    }
    return log_lik;
  }
}

// data and fixed parameter values
data {
  real<lower=0.> c1; // rate of loss of infectivity of infectious virus
  real<lower=0.> delta; // decay rate of infected cells
  real<lower=0.> T_0; // initial number of cells in each well
  real<lower=0.> V_0; // initial viral load
  real<lower=0.> A_0;
  real<lower=0.> p; // production rate of virus from infected cells
  real<lower=0., upper=1.> alpha; // initial proportion of cells in cell type 1
  real<lower=0.> q;
  real<lower=0.> phi;
  real<lower=0.> rho;
  real<lower=0.> d;
  real<lower=0.> m;
  real<lower=0.> k1;
  real<lower=0.> A_max;
  
  real<lower=0.> sigma_V; // standard deviation for observation model for viral load
  real<lower=0.> lod_V; // lower limit of detection for viral load
  real<lower=0., upper=1.> rho0; // overdispersion parameter
  
  int<lower=1> T_V; // total number of samples for viral load
  int<lower=1> T_C; // total number of samples for cell type
  int<lower=1> T_S; // number of sampling times excluding 0
  
  array[T_V] int t_V_idx; // indices linking sample number to sampling time for viral load
  array[T_C] int t_C_idx; // indices linking sample number to sampling time for cell counts
  vector[T_V] obs_V; // vector of viral load data
  array[2*2,T_C] int obs_C; // matrix of cell count data. columns are samples,
  //rows are T[1], T[2], ..., T[n_C], I[1], I[2], ..., I[n_C]
  array[T_S] real ts; // unique solving times
  
  int evaluate_likelihood; // if 1, sample from posterior, if 0, sample from prior
}

// model parameters
parameters {
  // ranges are ranges of uniform priors
  real<lower=-8.,upper=-4.> log10_beta; // mean value
  real<lower=-1., upper=1.> Delta_scaled;
}

transformed parameters {
  real log_lik;
  vector[2] beta1;
  real Delta;
  real mean_beta1;
  
  mean_beta1 = 10^log10_beta;
  
  if(Delta_scaled < 0) {
    Delta = Delta_scaled * (1 - alpha);
  } else {
    Delta = Delta_scaled * alpha;
  }
  
  beta1[1] =  mean_beta1 * (1 - Delta / alpha);
  beta1[2] =  mean_beta1 * (1 + Delta / (1 - alpha));
  
  log_lik = calc_log_lik(2, c1, delta, T_0, p, V_0, alpha, A_0, q, phi, rho, d, m, k1, A_max, sigma_V,
  lod_V, rho0, T_V, T_C, T_S,t_V_idx, t_C_idx, obs_V, obs_C, ts, beta1, evaluate_likelihood);
}

model {
  target += log_lik;
}
