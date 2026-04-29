#' simulate dirichlet multinomial observations
#' @param true_n_cells_vec vector of length 2 * n_C, where n_C is the number of cell types.
#' concatenated vector of number of uninfected cells of each type and
#' number of infected cells of each type in the well
#' @param rho overdispersion parameter between 0 and 1
#' @param sample_size sample size
#' @return n_sims x 2 * n_C matrix, where n_C is the number of cell types.
#' concatenated vector of number of uninfected cells of each type and
#' number of infected cells of each type in the well.
#' if n_sims == 1, returns vector instead
sim_dirichlet_multinomial <- function(true_n_cells_vec, rho, sample_size, n_sims = 1) {
  # https://en.wikipedia.org/wiki/Dirichlet-multinomial_distribution
  # alpha0 <- (1 - rho)/rho
  # normalise(true_n_cells_vec) * alpha0

  if(missing(sample_size)) sample_size <- sum(true_n_cells_vec)
  sims <- dirmult::simPop(J=n_sims,
                  n = sample_size,
                  pi = normalise(true_n_cells_vec),
                  theta = rho)
  sims <- sims$data
  if(n_sims == 1) {
    return(as.numeric(sims))
  }
  sims
}

#' solve TIV model with heterogeneous betas for different values of sigma
#' @param pars list of parameter values from get_two_compartment_pars
#' @param sigma scalar
#' @param solve_TIV_two_beta function to solve TIV model with heterogeneous betas
#' @return tibble with model solution
#' @export
change_sigma_wrapper <- function(solve_TIV_two_beta, pars, sigma) {
  pars$beta1 <- mean(pars$beta1) + sigma * c(-1, 1)
  solve_TIV_two_beta(pars) %>%
    mutate(sigma = sigma)
}

#' simulate data from model with heterogeneous beta
#' @param sol tibble with columns t and V -- solution for viral load
#' @param sampling_times vector of sampling times -- must be in sol$t
#' @param sampling_pars list if sampling parameters including lod_V (limit of detection) and sigma_V (sd of lognormal noise)
#' @param sim_filename file in which to store simulation
#' @param cells logical: whether to simulate observations for cells
#' @return tibble with model solution
#' @export
simulate_data <- function(sol, sampling_times, sampling_pars, cells, sim_filename) {
  
  time_tibble <- tibble(t = rep(sampling_times$V, each = sampling_pars$n_reps))

  obs_V <- sol %>%
    filter(t %in% sampling_times$V) %>%
    full_join(time_tibble) %>%
    mutate(obs_V = add_log10normal_noise(V, sampling_pars$sigma_V), # lognormal noise
           obs_V = ifelse(obs_V < sampling_pars$lod_V, 1, obs_V)) %>% # set values below detection limit to 1 
    select(t, V, obs_V)
  
  if(cells) {
    time_tibble <- tibble(t = rep(sampling_times$C, each = sampling_pars$n_reps))
    n_C <- sum(grepl("T1", colnames(sol)))

    cell_counts <- sol %>%
      filter(t %in% sampling_times$C) %>%
      full_join(time_tibble) %>%
      as.matrix
    
    # check whether R1 exists, if so, combine T and R
    
    exist_resistant_cells <- "R[1]" %in% colnames(cell_counts)
    
    if(exist_resistant_cells) {
      cell_compartment_names <- c(paste0("T1[", seq_len(n_C), "]"),
                                  paste0("R[", seq_len(n_C), "]"),
                                  paste0("I1[", seq_len(n_C), "]"))
      
      cell_counts <- cell_counts[, cell_compartment_names]
      cell_counts[,paste0("T1[", seq_len(n_C), "]")] <- cell_counts[,paste0("T1[", seq_len(n_C), "]")] +
        cell_counts[,paste0("R[", seq_len(n_C), "]")]
    }
    
    cell_compartment_names <- c(paste0("T1[", seq_len(n_C), "]"),
                                paste0("I1[", seq_len(n_C), "]"))
    
    cell_counts <- cell_counts[, cell_compartment_names]

    obs_C <- t(apply(cell_counts, 1, 
                     function(x) sim_dirichlet_multinomial(x, 
                                                           sampling_pars$rho, 
                                                           sample_size = sampling_pars$cell_sample_size)))
    
    colnames(obs_C) <- cell_compartment_names

    obs_C <- obs_C %>%
      as_tibble %>%
      bind_cols(time_tibble) %>%
      select(t, everything())
  } else {
    obs_C <- NULL
  }
  
  sim_data <- list(obs_V = obs_V, obs_C = obs_C)
  saveRDS(sim_data, sim_filename)
  invisible(NULL)
}

add_log10normal_noise <- function(true_values,sd_noise){
  vnapply(true_values, function(x) rl10norm(1,log10(x),sd_noise))
}

#' Analogue for rlnorm but in base 10
#'
#' \code{rl10norm} is an analogue for rlnorm but in base 10
#' 
#' @param n a numeric vector of length 1: number of observations
#' @param mean a vector of means.
#' @param sd a vector of standard deviations.
#' @return a numeric vector of length n, distributed according to the log10normal distribution.
#' The log10normal distribution in log10 space has mu =  mean and sigma = sd.

rl10norm <- function(n, mean = 1, sd = 1){
  mean <- mean[!is.na(mean)]
  sd <- sd[!is.na(sd)]
  10^rnorm(n,mean,sd)
}

#' function which imposes an observation threshold on data.
#' 
#' \code{impose_observation_threshold} imposes an observation threshold on data.
#' 
#' Run this after imposing other noise (e.g. lognormal)
#' 
#' @param noisy_values numeric vector of noisy values on which to impose threshold
#' @param threshold below this threshold, set values to -1 to denote below threshold
#' @return numeric vector of values with imposed threshold
impose_observation_threshold <- function(noisy_values,threshold){
  noisy_values[noisy_values < threshold] <- -1
  noisy_values
}

#' calculate R_0
#' @param pars named vector or list of parameter values with elements "beta1", "T_0", "p", "delta", "c1"
#' if there is heterogeneity in the model, pars needs to be a list, and beta1, T_0, p can be vectors
#' @return scalar: value of R_0
#' @export
calc_R_0 <- function(pars) {
  if(!inherits(pars, "list")) pars <- as.list(pars)
  # check dimensions are compatible
  n_C <- length(pars$T_0)
  stopifnot(length(pars$beta1) %in% c(1, n_C))
  stopifnot(length(pars$p) %in% c(1, n_C))
  sum(pars$beta1 * pars$T_0 * pars$p) / pars$delta / pars$c1
}

#' calculate R_0 for homogeneous case (faster)
#' @param beta1 value of beta1
#' @param T_0 value of T_0
#' @param p value of p
#' @param delta value of delta
#' @param c1 value of c1
#' @return scalar: value of R_0
#' @export
calc_R_0_homogeneous <- function(beta1, T_0, p, delta, c1) {
  beta1 * T_0 * p / delta / c1
}

#' calculate r for homogeneous model
#' @inheritParams calc_R_0
#' @return scalar: value of r
#' @export
calc_r <- function(pars) {
  if(!inherits(pars, "list")) pars <- as.list(pars)
  stopifnot(length(pars$p) == 1) # this expression doesn't hold for heterogeneous p
  beta_T_0 <- sum(pars$beta1 * pars$T_0)
  -(pars$c1 + pars$delta) / 2 + 
    sqrt((pars$delta - pars$c1)^2 + 4 * beta_T_0 * pars$p) / 2
}

#' calculate r for homogeneous case (faster)
#' @inheritParams calc_R_0_homogeneous
#' @return scalar: value of R_0
#' @export
calc_r_homogeneous <- function(beta1, T_0, p, delta, c1) {
  -(c1 + delta) / 2 + sqrt((delta - c1)^2 + 4 * beta1 * T_0 * p) / 2
}

#' calculate burst size for homogeneous model
#' @inheritParams calc_R_0
#' @return scalar: value of r
#' @export
calc_burst_size_homogeneous <- function(pars) {
  if(!inherits(pars, "list")) pars <- as.list(pars)
  pars$p / pars$delta
}


