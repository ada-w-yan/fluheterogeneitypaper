#' fit model to data using cmdstan
#' 
#' @param beta_T_0 scalar
#' @param p scalar
#' @param k scalar
#' @param delta scalar
#' @param c1 scalar
calc_r_teiv <- function(beta_T_0, p, k, delta, c1) {
  # linearise around DFE
  mat <- matrix(c(-k, 0, beta_T_0, k, -delta, 0, 0, p, -c1), nrow = 3)
  # find eigenvalues
  eigenvalues = eigen(mat)$values
  # find largest real eigenvalue
  max(Re(eigenvalues[abs(Im(eigenvalues)) < 1e-6])) # you have to choose the precision you like here
}

#' get model filename
#' @param model_name string
#' @return string: model filename
#' @export
get_model_filename <- function(model_name) {
  git_repo_dirname <- "../"
  paste0(git_repo_dirname, "model_stan/", model_name, ".stan")
}

#' get the number of iterations post warm-up for each chain in a stanfit object
#'
#' @param fit stanfit object
#' @return vector: iterations post warm-up for each chain
#' @export
get_n_iter <- function(fit) {
  vnapply(fit@stan_args, function(x) x$iter - x$warmup)
}

#' Change fig.width and fig.height Dynamically Within an R Markdown Chunk
#' http://michaeljw.com/blog/post/subchunkify/
#'
#' @param g a ggplot object
#' @param fig_height scalar: figure height
#' @param fig_width scalar: figure width
#' @return gpglot object with correct height and width
subchunkify <- function(g, fig_height=7, fig_width=5) {
  g_deparsed <- paste0(deparse(
    function() {g}
  ), collapse = '')

  sub_chunk <- paste0("
  `","``{r sub_chunk_", floor(runif(1) * 10000), ", fig.height=",
                      fig_height, ", fig.width=", fig_width, ", echo=FALSE}",
                      "\n(",
                      g_deparsed
                      , ")()",
                      "\n`","``
  ")

  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
}

#' adjust N_samples for predictions
#'
#' @param N_samples original value for N_samples
#' @param chain_length length of NUTS chain
#' @return integer: chain_length if N_samples == 0, min(N_samples, chain_length) otherwise
#' @export
adjust_N_samples <- function(N_samples, chain_length) {
  if(N_samples == 0 || (N_samples > chain_length)) {
    if(N_samples > chain_length) {
      warning("number of samples requested exceeds length of chain, using length of chain instead")
    }
    N_samples <- chain_length
  }
  N_samples
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
