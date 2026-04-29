#' read in inferred values of beta
#' @param filename_fn function which
#' returns the filename of the fit
#' @param ... arguments to be passed to filename_fn
#' @return tibble with columns iter, par_name, value, chain, sigma where par_name is beta1[1] or beta1[2]
#' @export
get_betas_old <- function(filename_fn, sigma) {
  get_chain <- Vectorize(\(x) strsplit(x, ".", fixed = TRUE)[[1]][[1]])
  fit <- filename_fn(sigma) %>% readRDS()
  n_C <- 2
  par_names <- paste0("beta1[", seq_len(n_C), "]")
  
  # check if beta1 was stored, if not, reconstruct from log10_beta1 and delta_beta
  # returns TRUE if fit$draws(beta_names[1]) evaluates, FALSE otherwise
  beta_stored <- tryCatch(!is.null(fit$draws(par_names[1])), error = function(e) FALSE)
  if(!beta_stored) {
    par_names <- c("log10_beta", paste0("delta_beta[", seq_len(n_C), "]"))
  }
  
  pars <- lapply(par_names, extract_par_values, fit = fit) %>%
    bind_rows() %>%
    mutate(chain = get_chain(chain),
           sigma = !!sigma)
  if(!beta_stored) {
    pars <- pars %>%
      pivot_wider(names_from = "par_name", values_from = "value") %>%
      mutate(across(starts_with("delta_beta"),
                    ~ 10^log10_beta * 2. * .x)) %>%
      rename_with(~str_replace(.x, "delta_beta", "beta1"), starts_with("delta_beta")) %>%
      select(-log10_beta) %>%
      pivot_longer(starts_with("beta"), names_to = "par_name")
  }
  pars
}

#' read in inferred values of beta
#' @param filename_fn function which
#' returns the filename of the fit
#' @param ... arguments to be passed to filename_fn
#' @return tibble with columns iter, par_name, value, chain, sigma where par_name is beta1[1] or beta1[2]
#' @export
get_betas <- function(filename_fn, ...) {
  get_chain <- Vectorize(\(x) strsplit(x, ".", fixed = TRUE)[[1]][[1]])
  fit_filename <- filename_fn(...)
  if(!file.exists(fit_filename)) return(NULL)
  fit <- fit_filename %>% readRDS()
  n_C <- 2
  par_names <- paste0("beta1[", seq_len(n_C), "]")
  
  id_tibble <- as_tibble(list(...))
  
  pars <- lapply(par_names, extract_par_values, fit = fit) %>%
    bind_rows() %>%
    mutate(chain = get_chain(chain)) %>%
    cross_join(id_tibble) # append the values of the ... arguments
  
  pars
}

#' read in inferred values of parameters
#' @param filename_fn function which returns the filename of the fit
#' @param par_names parmaeter names
#' @param ... arguments to be passed to filename_fn
#' @return tibble with columns iter, par_name, value, chain, sigma where par_name is beta1[1] or beta1[2]
#' @export
get_pars <- function(filename_fn, par_names, ...) {
  get_chain <- Vectorize(\(x) strsplit(x, ".", fixed = TRUE)[[1]][[1]])
  fit_filename <- filename_fn(...)
  if(!file.exists(fit_filename)) return(NULL)
  fit <- fit_filename %>% readRDS()
  id_tibble <- as_tibble(list(...))
  
  pars <- lapply(par_names, extract_par_values, fit = fit) %>%
    bind_rows() %>%
    mutate(chain = get_chain(chain)) %>%
    cross_join(id_tibble) # append the values of the ... arguments
  
  pars
}

#' take tibble of inferred beta values and enforce beta1 > beta2
#' @param betas tibble with columns iter, par_name, value, chain, sigma where par_name is beta1[1] or beta1[2]
#' @return tibble with columns iter, par_name, value, chain, sigma
#' @export
enforce_beta_larger <- function(betas) {
  grouping_names <- colnames(betas) %>% setdiff(c("par_name", "value"))
  betas %>%
    # group_by(iter, chain, sigma) %>%
    group_by(pick(grouping_names)) %>%
    summarise(beta1_temp = min(value),
              beta2_temp = max(value)) %>%
    ungroup() %>%
    rename(`beta1[1]` = beta1_temp,
           `beta1[2]` = beta2_temp) %>%
    pivot_longer(starts_with("beta"), names_to = "par_name")
}

#' calculate inferred mean value of beta
#' @param betas tibble with columns iter, par_name, value, chain, ... where par_name is beta1[1] or beta1[2]
#' @return tibble with columns iter, chain, ..., mean_beta
#' @export
calc_mean_beta <- function(betas) {
  grouping_names <- colnames(betas) %>% setdiff(c("par_name", "value"))
  betas %>%
    group_by(pick(grouping_names)) %>%
    summarise(mean_beta = mean(value)) %>%
    ungroup()
}

#' calculate inferred weighted mean value of beta (weighted by $T_0i$)
#' @param betas tibble with columns iter, par_name, value, chain, ... where par_name is beta1[1] or beta1[2]
#' @param samples tibble with columns alpha, ... where ... are the same indexing variables as in betas
#' @return tibble with columns iter, chain, weighted_mean_beta, ...
#' @export
calc_weighted_mean_beta <- function(betas, samples) {
  indexing_variables <- colnames(betas) %>% setdiff(c("iter", "chain", "par_name", "value"))
  samples <- samples %>% select(any_of(c(indexing_variables, "alpha")))

  betas %>%
    pivot_wider(names_from = "par_name", values_from = "value") %>%
    full_join(samples) %>%
    mutate(weighted_mean_beta = `beta1[1]` * alpha + `beta1[2]` * (1 - alpha))
}

#' sample from prior distribution for the mean of beta
#' @param n_C number of cell types
#' @return tibble with columns mean_beta (contains the samples)
#' @export
gen_prior_mean_beta <- function(n_C = 2, n_samples_prior = 1e4) {
  10^runif(n_samples_prior * n_C, -8, -4) %>%
    matrix(ncol = n_C) %>%
    rowMeans() %>%
    enframe(name = NULL, value = "mean_beta")
}

#' sample from prior distribution for the ratio of individual beta to beta
#' @param n_C number of cell types
#' @return tibble with columns par_name, value and sigma = "prior"
#' @export
gen_prior_indiv_beta <- function(n_C = 2, n_samples_prior = 1e4) {
  prior_samples <- 10^runif(n_samples_prior * n_C, -8, -4) %>%
    matrix(ncol = n_C)
  tibble(par_name = rep(paste0("beta1[", seq_len(n_C), "]"), n_samples_prior / 2),
         value = prior_samples[,1] / rowMeans(prior_samples),
         sigma = "prior")
}

#' plot mean betas
#' @param mean_betas tibble with columns iter, chain, sigma, mean_beta
#' @param true_value true_value of mean beta
#' @return plot
#' @export
plot_mean_betas_histogram <- function(mean_betas, true_value) {

  mean_betas <- mean_betas %>%
    select(-iter, -chain) %>%
    select(mean_beta, everything())
  
  n_dim <- dim(mean_betas)[2] - 1 # number of parameters being varied
  colnames(mean_betas) <- c("mean_beta", as.character(seq_len(n_dim)))
  
  prior_samples <- gen_prior_mean_beta() %>%
    cross_join(select(mean_betas, -mean_beta))
  
  g <- mean_betas %>%
    bind_rows(prior_samples) %>%
    ggplot() +
    geom_vline(xintercept = true_value) +
    geom_histogram(aes(x = mean_beta)) +
    theme_bw() +
    scale_x_log10("Mean value of beta", guide = guide_axis(angle = 90))
  
  if(n_dim == 1) {
    g <- g + facet_wrap(~`1`, nrow = 1)
  } else {
    g <- g + facet_grid(`1`~`2`)
  }
  
  g
  
}

#' plot ratio of individual beta to mean beta
#' @param betas tibble with columns iter, par_name, value, chain, sigma where par_name is beta1[1] or beta1[2]
#' @param mean_betas tibble with columns iter, chain, sigma, mean_beta
#' @param true_value true values of mean beta
#' @return plot
#' @export
plot_beta_ratios <- function(betas, mean_betas, true_values) {
  betas <- betas %>%
    full_join(mean_betas) %>%
    mutate(value = value / mean_beta) %>%
    select(par_name, value, sigma, chain) %>%
    mutate(sigma = as.character(sigma))
  
  subplot <- function(idx) {
    true_values <- true_values %>%
      filter(par_name == paste0("beta1[", idx, "]"))
    x_lim <- if(idx == 1) c(0, 1) else c(1, 2)
    
    g <- betas %>%
      filter(par_name == paste0("beta1[", idx, "]")) %>%
      ggplot() +
      facet_grid(par_name~sigma, scales = "free") +
      geom_histogram(aes(x = value)) +
      geom_vline(data = true_values, aes(xintercept = value)) +
      theme_bw() + 
      scale_x_continuous(limits = x_lim, "Ratio of individual beta to mean beta", guide = guide_axis(angle = 90))
    g
  }
  
  lapply(seq_len(2), subplot) %>%
    gridExtra::grid.arrange(grobs = ., ncol = 1)
}

#' make a tibble with the true values of beta1[1] and beta1[2]
#' @param beta1 scalar: value of beta1
#' @param sigma_vec: vector of values of sigma
#' @return tibble with columns par_name, value and sigma
#' @export
get_true_values_beta_tibble <- function(beta1, sigma_vec, n_C = 2) {
  tibble(par_name = rep(paste0("beta1[", seq_len(n_C), "]"), each = length(sigma_vec)),
         value = as.numeric((beta1 + outer(sigma_vec, c(-1, 1))) / beta1),
         sigma = rep(sigma_vec, n_C))
}

#' get model predictions for viral load
#'
#' @param dir_name directory name
#' @param filename filename
#' @param data_list data list for stan fitting
#' @param used_pred if TRUE, predictions were made after model fitting; if FALSE, predictions
#' were made alongside model fitting
#' @return tibble with columns idx, condition, t, V
get_sol_V <- function(filename, pred_times, used_pred) {

  pred <- readRDS(filename)
  pred <- pred$draws()
  pred_names <- dimnames(pred)$variable
  pred_values <- as.numeric(pred)
  if(used_pred) {
    names(pred_values) <- pred_names
  } else {
    names(pred_values) <- rep(pred_names, each = length(pred_values) / length(pred_names))
  }
  
  rm(pred)
  rm(pred_names)
  
  sol_V_idx <- grep("sol_V", names(pred_values))
  sol_V <- pred_values[sol_V_idx]
  sol_V <- tibble(name = names(sol_V),
                  value = sol_V) %>%
    # pivot_longer(everything()) %>%
    mutate(name = str_split_i(name, fixed("["), 2),
           name = str_replace_all(name, fixed("]"), ""),
           # the idx is stored as the first string if used_pred = TRUE, 
           # so it moves the string positions along by 1
           t_idx = str_split_i(name, fixed(","), 1 + used_pred))
  
  # the idx is stored as the first string if used_pred = TRUE
  if(used_pred) {
    sol_V <- sol_V %>%
      mutate(idx = as.numeric(str_split_i(name, fixed(","), 1)))
  } else {
    # otherwise recreate idx
    sol_V <- sol_V %>%
      group_by(t_idx) %>%
      mutate(idx = seq_len(n())) %>%
      ungroup()
  }
  
  sol_V %>%
    mutate(t_idx = as.numeric(t_idx),
           t = pred_times[t_idx]) %>%
    select(idx, t, value) %>%
    rename(V = value)
}

#' plot viral load from data and model fits
#' 
#' @param sol_V model prediction for viral load
#' @param sim_data simulated data
#' @return ggplot object
plot_single_sol <- function(sol_V, sim_data) {
  n_draws_to_plot <- 10
  n_draws <- max(sol_V$idx)
  n_draws_to_plot <- max(n_draws_to_plot, n_draws)
  plot_draws_idx <- round(seq(1, n_draws, length.out = n_draws_to_plot))
  
  quantiles <- sol_V %>%
    group_by(t) %>%
    summarise(lower = quantile(V, prob = 0.025),
              median = quantile(V, prob = 0.5),
              upper = quantile(V, prob = 0.975)) %>%
    ungroup()
  
  sol_V <- sol_V %>%
    filter(idx %in% plot_draws_idx)

  grey <- "#D8D8D8"
  ggplot(quantiles) +
    geom_ribbon(aes(x = t, ymin = lower, ymax = upper), fill = grey) +
    geom_point(data = sim_data, aes(x = t, y = obs_V)) +
    geom_line(data = sol_V, aes(x = t, y = V, group = idx), color = "darkgrey") +
    geom_line(data = sim_data, aes(x = t, y = V)) +
    theme_bw() +
    scale_y_log10("Viral load (pfu)", limits = c(1, 1e10)) +
    xlab("Time (days)")
}

#' plot model predictions for viral load
#' @param pred tibble containing model predictions.
#' contains columns t, lower (2.5% quantile of prediction), 
#' upper (97.5% quantile of prediction), max_LL (maximum likelihood prediction), sigma
#' @param all_sim_data tibble with simulated data.
#' contains columns t, obs_V, sigma
#' @param lod_V lower limit of detection for viral load
#' @return ggplot object
#' @export
plot_pred_V <- function(pred, all_sim_data, lod_V) {
  pred %>%
    filter(compartment_name == "V") %>%
    ggplot(aes(x = t)) +
    facet_wrap(~sigma, nrow = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#D8D8D8") +
    geom_line(aes(y = max_LL)) +
    geom_point(data = all_sim_data, aes(y = obs_V)) +
    geom_hline(yintercept = sampling_pars$lod_V, linetype = "dotted") +
    theme_bw() +
    scale_y_log10("Viral load (pfu/mL)", limits = c(1, 1e9)) +
    xlab("Time (days)")
}

#' plot model predictions for cell compartments
#' @param pred tibble containing model predictions.
#' contains columns t, lower (2.5% quantile of prediction), 
#' upper (97.5% quantile of prediction), max_LL (maximum likelihood prediction), sigma
#' @param all_sim_data tibble with simulated data.
#' contains columns t, obs_V, sigma
#' @param lod_V lower limit of detection for viral load
#' @return ggplot object
#' @export
plot_pred_C <- function(pred, all_sim_data, plot_proportion) {
  if(plot_proportion) {
    pred <- pred %>%
      filter(grepl("prop", compartment_name)) %>%
      mutate(compartment_name = gsub("prop_", "", compartment_name))
    all_sim_data <- all_sim_data %>%
      select(-value) %>%
      rename(value = prop_value)
  } else {
    pred <- pred %>%
      filter(!grepl("prop", compartment_name), compartment_name != "V")
  }
  
  g <- pred %>%
    ggplot(aes(x = t)) +
    facet_grid(compartment_name~sigma, scales = "free") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#D8D8D8") +
    geom_line(aes(y = max_LL)) +
    geom_point(data = all_sim_data, aes(y = value)) +
    theme_bw() +
    xlab("Time (days)")
  if(plot_proportion) {
    g <- g +
      scale_y_log10("Proportion of sampled cells") +
      coord_cartesian(xlim = c(0, 2), ylim = c(1e-5, 1))
  } else {
    g <- g +
      scale_y_log10("Number of cells") +
      coord_cartesian(xlim = c(0, 2), ylim = c(1, 1e8))
  }
  g
}

#' make model predictions
#' @param filename_fn function to get filename where parameter values are stored.  Takes argument sigma.
#' @param make_pred_fn function to make model predictions. Has arguments idx, beta1 and beta2,
#' and returns a tibble with columns idx, t and compartment_name and value.
#' @param sigma scalar value of sigma.
#' @return pred tibble containing model predictions.
#' contains columns t, lower (2.5% quantile of prediction),
#' upper (97.5% quantile of prediction), max_LL (maximum likelihood prediction), sigma
#' @export
make_pred_beta <- function(filename_fn, make_pred_fn, sigma) {
  fit <- filename_fn(sigma) %>% readRDS()
  beta_wide <- get_betas(sigma, filename_fn) %>%
    mutate(idx = paste(iter, chain, sigma)) %>%
    select(idx, par_name, value) %>%
    mutate(par_name = gsub("1[", "", par_name, fixed = TRUE),
           par_name = gsub("]", "", par_name, fixed = TRUE)) %>%
    pivot_wider(names_from = par_name, values_from = value)
  
  pred <- apply_named_args(beta_wide, 1, make_pred_fn) %>%
    bind_rows()
  
  # calculate proportion of cells in each compartment at each timepoint
  calc_prop_wrapper <- function(pred) {
    pred %>%
      filter(grepl("]", compartment_name, fixed = TRUE)) %>%
      rename(replicate = idx) %>%
      mutate(sigma = !!sigma) %>%
      calc_prop() %>%
      select(-sigma, -value) %>%
      mutate(compartment_name = paste0("prop_", compartment_name)) %>%
      rename(value = prop_value)
  }
  
  prop_cells <- calc_prop_wrapper(pred)
  
  pred_CI <- pred %>%
    bind_rows(prop_cells) %>%
    group_by(t, compartment_name) %>%
    summarise(lower = quantile(value, probs = .025),
              upper = quantile(value, probs = .975)) %>%
    ungroup() %>%
    arrange(sigma, t, compartment_name)
  
  max_LL_pars <- get_max_LL_pars(fit)
  if("log10_beta1[1]" %in% names(max_LL_pars)) {
    pred_max_LL <- make_pred_fn(NA,
                                beta1 = 10^max_LL_pars[["log10_beta1[1]"]],
                                beta2 = 10^max_LL_pars[["log10_beta1[2]"]])
  } else if("delta_beta[1]" %in% names(max_LL_pars)) {
    pred_max_LL <- make_pred_fn(NA,
                                beta1 = 10^max_LL_pars[["log10_beta"]] * 2. * max_LL_pars[["delta_beta[1]"]],
                                beta2 = 10^max_LL_pars[["log10_beta"]] * 2. * max_LL_pars[["delta_beta[2]"]])
  }
  
  
  prop_cells <- calc_prop_wrapper(pred_max_LL)
  
  pred_max_LL <- pred_max_LL %>%
    bind_rows(prop_cells) %>%
    arrange(sigma, t, compartment_name) %>%
    select(value) %>%
    rename(max_LL = value)
  
  pred <- bind_cols(pred_CI, pred_max_LL) %>%
    mutate(sigma = !!sigma)
  pred
}

#' make traceplot of the log likelihood
#' @param filename_fn function to get filename where parameter values are stored.  Takes argument sigma.
#' @param sigma scalar value of sigma.
#' @return ggplot object
#' @export
plot_trace_log_lik <- function(filename_fn, sigma) {
  fit <- filename_fn(sigma) %>%
    readRDS()
  fit$draws() %>%
    mcmc_trace(pars = "log_lik") +
    ggtitle(sigma) +
    theme(legend.position = "none") +
    scale_x_continuous(guide = guide_axis(angle = 90))
}

#' calculate proportion of each cell type and infection status for each time point and replicate
#' @param sol tibble with columns t, sigma, replicate, compartment_name and value
#' @return tibble with columns t, sigma, replicate, compartment_name and value
#' @export
calc_prop <- function(sol) {
  sol %>%
    group_by(t, sigma, replicate) %>%
    mutate(prop_value = value / sum(value)) %>%
    ungroup()
}
