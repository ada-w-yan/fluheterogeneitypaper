#' get model filename
#' @param model_name string
#' @return string: model filename
#' @export
get_model_filename <- function(model_name) {
  git_repo_dirname <- "../"
  paste0(git_repo_dirname, "model_stan/", model_name, ".stan")
}

#' append indices to parameter name
#'
#' @param par_name string: parameter name
#' @param idx numeric vector
#' @return string.  e.g. if par_name == "abc" and idx == c(1,2), returns "abc[1,2]"
#' @export
index_par_name <- function(par_name, idx) {
  idx <- paste(idx, collapse = ",")
  paste0(par_name, "[", idx, "]")
}

#' extract model predictions as matrix
#'
#' @param pred stanfit object
#' @param N_samples number of samples used to make model predictions
#' @param pred_name (optional) name of model prediction.  If not specified, return everything
#' @return matrix of model predictions
#' @importFrom dplyr %>%
#' @export
extract_pred <- function(pred, N_samples, pred_name) {
  pred <- as.matrix(pred)
  if(!missing(pred_name)) {
    pred <- pred[,grepl(pred_name, colnames(pred), fixed = TRUE)]
  } else {
    pred <- pred[1, -ncol(pred)]
  }
  pred <- matrix(pred, ncol = N_samples)
  pred
}

#' extract samples from fit object and thin
#'
#' @param fit stanfit object
#' @param N_samples number of samples to extract.  if 0 or greater than the number of available samples, extract all
#' @param par_names (optional) character vector: namse of parameters to extract.  if missing, extract all
#' @param drop if par_names is specified and is of length 1, and drop = TRUE, return vector;
#' otherwise return matrix
#' @param by_chain if TRUE, return list of matrices with one matrix per chain.  Then N_samples is interpreted as per
#' chain.  if FALSE, return one matrix for the whole chain.
#' @return matrix with named columns containing samples
#' @export
extract_fit <- function(fit, N_samples, par_names, drop = FALSE, by_chain = FALSE) {
  ext_fit <- as.matrix(fit)
  if(missing(par_names)) {
    par_names <- colnames(ext_fit)
  }
  if(by_chain) {
    N_chains <- length(fit@inits)
  } else {
    N_chains <- 1
  }
  chain_length <- nrow(ext_fit) / N_chains
  N_samples <- adjust_N_samples(N_samples, chain_length)
  thinned_samples <- round(seq(1, chain_length, length.out = N_samples))
  if(by_chain) {
    ext_fit <- lapply(seq_len(N_chains),
                      function(x) ext_fit[seq(chain_length * (x - 1) + 1, x * chain_length),])
    ext_fit <- lapply(ext_fit, function(x) x[thinned_samples, par_names, drop = drop])
  } else {
    ext_fit <- ext_fit[thinned_samples, par_names, drop = drop]
  }
  
  ext_fit
}

#' get the maximum likelihood parameters from a stan fit
#'
#' @param fit_obj a stanfit object, matrix containing draws from posterior, or output from cmdstanr
#' @return named vector with parameters and generated quantities of maximum likelihood
#' @export
get_max_LL_pars <- function(fit_obj) {
  if(inherits(fit_obj, "CmdStanFit")) {
    draws <- fit_obj$draws() # draws are stored in a 3d array
    # log_lik is stored in a 2d array
    # due to background transforms in Stan, lp__ is the log likelihood with
    # a function added to it, rather than just the long likelihood.
    # this is less of a problem if the data is informative.
    # ideally we computed the log likelihood and stored it in log_lik,
    # but if not we use lp__ as an approximation.
    
    log_lik <- tryCatch(draws[,,"log_lik"], error = function(e) draws[,,"lp__"])
    
    max_log_lik <- max(log_lik)
    which_max_log_lik <- which.max(log_lik)[1] # break ties
    which_max_log_lik_idx1 <- ceiling(which_max_log_lik / dim(log_lik)[1])
    which_max_log_lik_idx2 <- which_max_log_lik - (which_max_log_lik_idx1 - 1) * dim(log_lik)[1]
    
    max_LL_values <- draws[which_max_log_lik_idx2, which_max_log_lik_idx1,] %>%
      as.numeric()
    names(max_LL_values) <- dimnames(draws)$variable
    return(max_LL_values)
  }
  
  if(inherits(fit_obj, "stanfit")) {
    fit_obj <- as.matrix(fit_obj)
  }
  
  log_lik_str <- if("log_lik" %in% colnames(fit_obj)) "log_lik" else "lp__"
  max_LL_row <- fit_obj[which.max(fit_obj[,"lp__"])[1],]
  max_LL_pars <- as.numeric(max_LL_row)
  names(max_LL_pars) <- colnames(fit_obj)
  max_LL_pars
}

#' new function for extracting parameter values (to do: compare with old)
#'
#' @param fit fit object from RStan or cmdstanr
#' @param par_name parameter names in fit
#' @return tibble with columns iter, chain, value, par_name
extract_par_values <- function(fit, par_names) {
  get_chain <- Vectorize(\(x) strsplit(x, ".", fixed = TRUE)[[1]][[1]])
  get_par_name <- Vectorize(\(x) strsplit(x, ".", fixed = TRUE)[[1]][[2]])
  use_cmdstan <- inherits(fit, "CmdStanFit")
  if(use_cmdstan) {
    fit <- fit$draws(par_names)
  } else {
    fit <- as.array(fit$fit)
    fit <- fit[,,par_names]
  }
  pars <- fit %>%
    as_tibble() %>%
    mutate(iter = seq_len(n())) %>%
    pivot_longer(-iter)
  if(use_cmdstan) {
    pars <- pars %>%
      mutate(par_name = get_par_name(name),
             chain = get_chain(name)) %>%
      select(iter, chain, value, par_name)
  } else {
    pars <- pars %>%
      rename(par_name = name)
  }
  pars
  
}

#' new function for plotting histograms (to do: compare with old)
#'
#' @param fit fit object from RStan or cmdstanr
#' @param par_names vector of parameter names in fit to plot
#' @param by_chain logical. Fill histogram by chain or not.
#' @param true_values (optional) vector of true values.  Must be of same length as par_names
#' @param prior_min (optional) vector of x axis lower limits.  Must be of same length as par_names
#' @param prior_max (optional) vector of x axis upper limits.  Must be of same length as par_names
#' @return ggplot object
plot_hist_stan <- function(fit, par_names, by_chain = TRUE, true_values, prior_min, prior_max) {
  
  if(!missing(true_values)) {
    true_values = tibble(par_name = par_names,
                         true_value = true_values)
  }
  
  if(!missing(prior_min) & !missing (prior_max)) {
    prior <- tibble(par_name = par_names,
                    prior_min = prior_min,
                    prior_max = prior_max)
  }
  
  
  g <- lapply(par_names, extract_par_values, fit = fit) %>%
    bind_rows() %>%
    ggplot() +
    facet_wrap(~par_name, nrow = 1, scales = "free") +
    theme_bw()
  
  if(by_chain) {
    g <- g + geom_histogram(aes(x = value, fill = chain))
  } else {
    g <- g + geom_histogram(aes(x = value))
  }
  
  if(!missing(true_values)) {
    g <- g +
      geom_vline(data = true_values, aes(xintercept = true_value))
  }
  
  if(!missing(prior_min) & !missing (prior_max)) {
    g <- g +
      geom_blank(data = prior, aes(x = prior_min)) +
      geom_blank(data = prior, aes(x = prior_max))
  }
  g
}

#' #' make model predictions
#' #'
#' #' @param model stan model object
#' #' @param data_list list of data and model parameters to feed to model
#' #' @return stanfit object
#' #' @export
#' make_pred <- function(model, data_list) {
#'   if(is.character(model)) {
#'     model <- rstan::stan_model(get_model_filename(model))
#'   }
#'   pred <- rstan::sampling(model,
#'                           data = data_list,
#'                           chains = 1, iter = 1,
#'                           algorithm = "Fixed_param")
#'   pred
#' }
#' 

#' choose subset of 8 chains from MCMC so that there are no divergences and 
#' Rhat is below 1.01 (if possible) or minimised
#' 
#' @param idx index of simulation estimation parameter set
#' @param Delta_scaled value of delta_scaled used for simulation
#' @param fit_wrapper function to fit model to data (in practice used to retrive filename)
#' @return tibble with columns idx (from input); Delta_scaled (from input);
#' prop_divergent (proportion of divertent transitions in selected chains);
#' max_rhat (maximum rhat calculated using selected chains);
#' chosen_chain_combn (indices of chosen chains)
get_summaries_subset_chains <- function(idx, Delta_scaled, fit_wrapper) {
  
  get_subset_chains <- function(idx, Delta_scaled) {
    fit_filename <- fit_wrapper(idx, Delta_scaled, eval = FALSE)
    if(!file.exists(fit_filename)) {
      return(list(fit = NA, prop_divergent = NA, max_rhat = NA, chosen_chain_combn = NA))
    }
    fit <- readRDS(fit_filename)
    
    n_draws <- length(fit$lp())
    chain_divergences <- fit$diagnostic_summary()$num_divergent
    n_chains <- length(chain_divergences)
    prop_divergent <- sum(chain_divergences) / n_draws
    max_rhat <- fit$summary() %>% pull(rhat) %>% max()
    
    # subset chains if there are convergence issues
    fit <- as_draws(fit)
    exclude_chain <- 0
    
    max_rhat_thres <- 1.01
    if(prop_divergent > 0 || max_rhat > max_rhat_thres) {
      # try including all combinations of chains of at least length 3
      # first test chains without divergences
      
      # make all combinations of chains of at least length 3
      # if there are fewer than 3 chains with no divergences, return NAs
      n_no_divergences <- length(which(chain_divergences == 0))
      if(n_no_divergences < 3) {
        return(list(fit = NA, prop_divergent = NA, max_rhat = NA, chosen_chain_combn = NA))
      }
      chain_combns <- do.call("c", lapply(seq(3, n_no_divergences), function(i) combn(which(chain_divergences == 0), i, FUN = list)))
      # aim to include largest number of chains
      chain_combns <- rev(chain_combns)
      
      # test chain combinations from longest to shortest and stop if rhat is below the threshold
      get_max_rhat <- function(chain_combn) {
        subset_fit <- posterior::subset_draws(fit, chain = chain_combn)
        summary_subset <- posterior::summarise_draws(subset_fit)
        max_rhat <- max(summary_subset$rhat)
        max_rhat
      }
      
      chosen_chain_combn <- NULL
      current_min_rhat <- Inf
      current_min_rhat_idx <- 0
      for(i in seq_along(chain_combns)) {
        max_rhat <- get_max_rhat(chain_combns[[i]])
        if(max_rhat < current_min_rhat) {
          current_min_rhat <- max_rhat
          current_min_rhat_idx <- i
        }
        if(max_rhat < max_rhat_thres) {
          chosen_chain_combn <- chain_combns[[i]]
          break
        }
      }
      
      # if no chain combinations have rhat < 1.01, choose the one with the lowest max rhat
      if(is.null(chosen_chain_combn)) chosen_chain_combn <- chain_combns[[current_min_rhat_idx]]
      
      fit <- posterior::subset_draws(fit, chain = chosen_chain_combn)
      chain_divergences <- chain_divergences[chosen_chain_combn]
      prop_divergent <- sum(chain_divergences) / n_draws * length(chosen_chain_combn) / n_chains
      max_rhat <- current_min_rhat
    } else {
      chosen_chain_combn <- seq_len(n_chains)
    }
    
    list(fit = fit, prop_divergent = prop_divergent, max_rhat = max_rhat, chosen_chain_combn = chosen_chain_combn)
  }
  
  fit <- get_subset_chains(idx, Delta_scaled)
  
  tibble(idx = idx, Delta_scaled = Delta_scaled, 
         prop_divergent = fit$prop_divergent, 
         max_rhat = fit$max_rhat,
         chosen_chain_combn = fit$chosen_chain_combn)
}
