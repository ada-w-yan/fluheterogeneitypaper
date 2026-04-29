#' get default parameter values
#' @return named vector of parameter values
#' @export
get_default_pars <- function() {
  pars <- c(p = 0.35,
            c1 = 4, # reduced so antibodies have more of an effect
            q = 1e-4, # increased q as there is only one innate immune mechanism
            d = 2,
            rho = 0.05,
            phi = .14,
            beta1 = 5e-6, # reduce beta to keep R_0 the same when c is reduced
            T_0 = 7e7,
            V_0 = 1,
            A_0 = 0.34,
            A_max = 1e4,
            delta = 3,
            m = 1,
            k1 = 1)
  
  pars
}

#' get directory of git repository
#' @param spartan logical: running on cluster or not
#' @return directory name
#' @export
get_git_dir <- function(spartan) {
  if(spartan) {
    git_dir <- "/data/gpfs/projects/punim0053/ayan/git_repos/fluheterogeneousmodels/"
  } else {
    git_dir <- "~/git_repos/fluheterogeneousmodels/"
  }
  git_dir
}

#' get default solving times
#' @return vector of solving times
#' @export
get_solving_time <- function() {
  seq(0, 7, by = .1)
}

#' get two compartment parameters
#' @param pars named vector of parameter values
#' @return list of parameter values
#' @export
get_two_compartment_pars <- function(pars) {
  pars_two_compartment <- as.list(pars)
  pars_two_compartment$beta1 <- rep(pars[["beta1"]], 2)
  pars_two_compartment$alpha <- .5
  pars_two_compartment <- pars_two_compartment[names(pars_two_compartment) != "I_0"]
  pars_two_compartment
}

#' get default values for Delta beta
#' @return vector of values for Delta beta
#' @export
get_sigma_vec <- function() {
  1e-7*seq(0, 4)
}

#' get parameter values for sampling
#' @return list of parameter values
#' @export
get_sampling_pars <- function() {
  list(sigma_V = 0.5, # standard deviation of observation model for viral load in log10 space
       lod_V = 10, # limit of detection for viral load
       rho = 0.01, # overdispersion of cell type counts
       n_reps = 3, # number of replicates)
       cell_sample_size = 5000) # number of cells in sample 
}

#' get times for sampling
#' @return named list of sampling times for cells and viral load
#' @export
get_sampling_times <- function() {
  list(C = seq(0, 2),
       V = seq(0, 7))
}