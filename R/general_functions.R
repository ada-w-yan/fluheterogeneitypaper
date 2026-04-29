# complete doc

#' checks if file exists in rds or csv format, and reads rds file preferentially
#'
#' checks if file exists in rds or csv format, and reads from original extension if given;
#' else reads rds file preferentially
#'
#' @param filename character vector of length 1: with or without extension
#' @return the read data frame
read_csv_or_rds <- function(filename,...){
  # check if extension given; if so, remove and keep track of original extension
  original_ext <- NULL
  if(as.logical(match(substr(filename,nchar(filename)-3,nchar(filename)),
                      c(".rds", ".csv"), nomatch = 0))){
    original_ext <- substr(filename,nchar(filename)-3,nchar(filename))
    filename <- substr(filename,1,nchar(filename)-4)
  }

  if(file.exists(paste0(filename,".rds"))){
    if(file.exists(paste0(filename,".csv"))){
      # if file exists in both extensions and original extension was not given, read rds
      if(is.null(original_ext)){
        warning(paste0(filename, " exists in both rds and csv -- reading rds"))
        try_read(filename,rds = TRUE,...)
      } else if(original_ext == ".rds"){
        # if file exists in both extensions and original extension was  given,
        # read original extension
        try_read(filename,rds = TRUE,...)
      } else {
        try_read(filename,rds = FALSE,...)
      }
    } else {
      # if file only exists in one extension, read that extension
      try_read(filename,rds = TRUE,...)
    }
  } else if (file.exists(paste0(filename,".csv"))){
    try_read(filename,rds = FALSE,...)
    # throw error if file doesn't exist in either extension
  } else {
    stop(paste0(filename, " not found"))
  }
}

#' reads rds or csv with tryCatch
#'
#' reads rds or csv with tryCatch
#'
#' rationale for existence: if it's not known beforehand whether we're reading
#'  an rds or csv, may pass nonsensical arguments, so ignore if these are causing an error
#'
#' @param filename character vector of length 1: with or without extension
#' @return the read data frame
try_read <- function(filename, rds, ...){
  tryCatch({
    if(rds){
      readRDS(paste0(filename,".rds"),...)
    } else {
      try_fread(paste0(filename,".csv"),...)
    }
  }, error = function(c){
    if(rds){
      readRDS(paste0(filename,".rds"))
    } else {
      try_fread(paste0(filename,".csv"))
    }
    warning(paste0("Additional arguments ignored when reading ", filename))
  })
}

#' data.table::fread with tryCatch
#'
#' data.table::fread sometimes crashes.  If crashing, try read.csv then converting.
#'
#' @param filename character vector of length 1 ending in .csv
#' @return the read data frame
try_fread <- function(filename, ...) {
  tryCatch ({
    data.table::fread(filename, ...)
  }, error = function(c) {
    data_table <- read.csv(filename, ...)
    data.table::as.data.table(data_table)
  })
}

#' thins a csv or rds file, saves it with new name, optionally deletes old file
#'
#' thins a csv or rds file, saves it with new name, optionally deletes old file
#'
#' @param filename character vector of length 1: with .csv or .rds extension
#' @param thin numeric vector of length 1 (must be in integer if not integer format):
#' thin every this many iterations
#' @param save_as_rds logical vector of length 1: if TRUE, save thinned chain as
#' rds regardless of extension of filename; if FALSE, save as original extension
#' @param remove_after logical vector of length 1: remove old file if TRUE
#' @return NULL
thin_csv_or_rds <- function(filename, thin, save_as_rds, remove_after){

  # check if thin is an integer
  if(!all.equal(thin, as.integer(thin))){
    stop(paste0("thin = ", as.character(thin), ": should be integer"))
  }

  # read file
  filename_wo_ext <- substr(filename, 1, nchar(filename)-4)
  ext <- substr(filename, nchar(filename)-3, nchar(filename))
  if(ext == ".rds"){
    table <- readRDS(filename)
  } else if(ext == ".csv"){
    table <- data.table:: fread(filename)
  } else {
    stop(paste0(filename, " is neither csv nor rds"))
  }

  #thin
  table <- table[seq(1, nrow(table), by = thin),]

  #save thinned file
  if(save_as_rds || ext == ".rds"){
    saveRDS(table, paste0(filename_wo_ext, "_thin_", as.character(thin),".rds"))
  } else{
    write.table(table, paste0(filename_wo_ext, "_thin_", as.character(thin),ext),
                row.names=FALSE,col.names=TRUE,sep=",",append=FALSE)
  }

  # remove old file
  if(remove_after){
    file.remove(filename)
  }
  invisible(NULL)
}

#' saves all variables in parent environment
#'
#' same as save.image except for parent environment rather than global
#'
#' @param filename character vector of length 1: no extension
#' @return NULL
save_all <- function(filename){
  save(list = ls(all.names = TRUE, envir = parent.frame()), file = paste0(filename,".RData"), envir = parent.frame())
}

#' makes a directory from a filename
#'
#' makes a directory from a filename
#'
#' @param filename character vector of length 1
#' @return NULL
mkdir_from_filename <- function(filename){
  dir_name <- dirname(filename)
  if(!dir.exists(dir_name)){
    dir.create(dir_name,recursive = TRUE)
  }
  invisible(NULL)
}

#' calculate geometric mean of a vector
#'
#' @param x numeric vector
#' @param na.rm logical vector of length 1
#' @return numeric vector of length 1: geometric mean
gm_mean <- function(x, na.rm=TRUE){
  if(na.rm) {
    denom <- sum(!is.na(x))
  } else {
    denom <- length(x)
  }
  exp(sum(log(x[x > 0]), na.rm=na.rm) / denom)
}

#' a version of Map() where we specify the format of the output
#'
#' @param FUN function
#' @param FUN.VALUE prespecified type of return value
#' @return vector with prespecified type of return value
Map_vapply <- function(FUN, FUN.VALUE, ...){
  Map_output <- Map(FUN, ...)
  vapply(Map_output, identity, FUN.VALUE)
}

#' solve an ODE with error handling
#'
#' if an error occurs, solve with a smaller time interval
#' for solution elements except the time vector, set to tol if the solution is below tol
#'
#' @param mod an ode_system object generated by odin
#' @param solving_time times at which to solve the ODE
#' @param tol numeric vector of length 1: tolerance
#' @return a deSolve object: solution of the ODE
solve_ODE_error_handling <- function(mod, solving_time, tol = 1e-6) {

  ## maximum number of times to try to solve with smaller time interval
  max_fail_iter <- 5

  sol <- tryCatch(mod$run(solving_time), error = function(e) "error") # solve ODEs

  ## if error occurred while solving ODEs, try smaller solving interval
  if(!is.matrix(sol)) {
    fail_iter <- 0
    n_div <- 10
    solving_time_temp <- solving_time
    while(fail_iter < max_fail_iter && !is.matrix(sol)) {
      # subdivide each time step into n_div steps
      solving_time_temp <- interpolate_vector(solving_time_temp, n_div)
      # resolve
      sol <- tryCatch(mod$run(solving_time_temp), error = function(e) "error")
      fail_iter <- fail_iter + 1
    }
    # retrive solution for original time steps
    sol <- sol[solving_time_temp %in% solving_time,]
  }

  # if ODE solver exits early, pad with last value
  sol <- sol[sol[,1] %in% solving_time,]
  if(nrow(sol) < length(solving_time)) {
    last_row <- sol[nrow(sol),]
    rep_last_rows <- rep(list(last_row), length(solving_time) - nrow(sol))
    sol <- do.call(rbind, c(list(sol), rep_last_rows))
    sol[,1] <- solving_time
  }

  # for solution elements except the time vector, set to tol if the solution is below tol
  non_time_elements <- sol[,-1]
  non_time_elements[non_time_elements < tol] <- tol
  sol[,-1] <- non_time_elements

  sol
}

#' interpolate ordered numeric vector into n_div intervals betwen successive values
#'
#' @param vec numeric vector to interpolate
#' @param n_div numeric vector of length 1: number of intervals into which to subdivide
#' successive values of the vector
#' @return numeric vector of length (length(vec) - 1) * n_div + 1: interpolated vector
interpolate_vector <- function(vec, n_div) {
  x <- seq(1, length(vec), by = 1/n_div)
  temp <- approx(seq_along(vec), vec, xout = x)
  vec_out <- temp$y
  vec_out
}

#' a version of lapply() that supplies FUN with both the name and the value of each component
#'
#' @param X see documentation for lapply
#' @param FUN function which takes a character vector of length 1 as a first argument and
#' something else as a second argument
#' @return see documentation for lapply
lapply_w_name <- function(X, FUN){
  Map(FUN, names(X), unname(X))
}

#' a version of a function that suppresses warnings
#'
#' @param f a function
#' @return the same function without warnings
no_warn <- function(f) {
  function(...) {
    suppressWarnings(f(...))
  }
}

#' make string from number in format acceptable for filename
#'
#' \code{make_filename_from_number} makes a string from a number, replacing
#' periods and plus signs
#' @param x numeric vector of length 1
#' @param decimal_points numeric vector of length 1. number of decimal points to
#'   round to
#' @param scientific logical vector of length 1.  whether to use scientific
#'   notation
#' @return character vector of lenth 1: string suitable for use as filename
make_filename_from_number <- function(x, decimal_points = 3, scientific = FALSE)
{
  if(scientific) {
    formatted_x <- format(signif(x, decimal_points + 1), scientific = TRUE)
  } else {
    formatted_x <- round(x, digits = decimal_points)
  }
  formatted_x <- sub(".", "point", formatted_x, fixed = TRUE)
  formatted_x <- sub("+", "", formatted_x, fixed = TRUE)
  formatted_x
}


#' load contents of RData file into a list
#'
#' @param filename character vector of length 1: RData filename
#' @param var_names optional argument specifying variables to load.  If missing,
#' load all variables
#' @return list containing variables with names var_names
load_into_list <- function(filename, var_names) {
  load_env <- new.env()
  load(filename, envir = load_env)
  output <- as.list(load_env)
  if(!missing(var_names)) {
    output <- get_vars_from_list_with_check(output, var_names)
  }
  output
}

#' gather variables with given names from an environment into a list
#'
#' @param var_names character vector of variable names to gather
#' @param envir environment in which to find parameters
#' @return list of variables with names var_names
list_vars_from_environment <- function(var_names, envir = parent.frame()) {
  env_vars <- as.list(envir)
  env_vars <- get_vars_from_list_with_check(env_vars, var_names)
  env_vars
}

#' dump variables from a list to the parent frame
#'
#' @param x list containing variables
#' @param var_names character vector containing names of variables to dump into
#' parent frame
#' @param overwrite optional parameter controlling the behaviour if any variables
#' with the same name(s) alraedy exist in the parent frame. If set to "warn",
#' throws a warning; if set to "error", throws an error
#' @return NULL
list2here <- function(x, var_names, overwrite) {
  if(!missing(var_names)) {
    x <- get_vars_from_list_with_check(x, var_names)
  }


  if(!missing(overwrite)) {
    parent_frame_vars <- ls(parent.frame)
    overwriting_vars <- intersect(names(x), parent_frame_vars)
    if(length(overwriting_vars) > 0) {
      if(overwrite == "warn") {
        lapply(overwriting_vars, function(x) warning(cat("overwriting", x)))
      } else if(overwrite == "error") {
        stop(cat("Attempting to overwrite", x[1]))
      }
    }
  }

  list2env(x, envir = parent.frame())
  invisible(NULL)
}

#' extract variables from list, throwing an error if they are not found
#'
#' @param x list of variables
#' @param var_names character vector containing names of variables to extract
#' @return list of selected variables
get_vars_from_list_with_check <- function(x, var_names) {
  missing_vars <- var_names[!(var_names %in% names (x))]
  if(length(missing_vars) > 0) {
    stop(cat("variables missing from list: ", paste0(missing_vars, collapse = " ")))
  }
  x <- x[var_names]
  x
}

#' make factor out of a numeric vector
#'
#' @param vec numeric vector
#' @return a character version of the vector, with levels in numeric order
make_factor <- function(vec) {
  factor(as.character(vec), levels = unique(as.character(sort(vec))))
}

suppress_final_line_warning <- function(w) {
  if(any(grepl("incomplete final line found on", w, fixed = TRUE))) {
    invokeRestart("muffleWarning")
  }
}

read_code_unspace <- function(filename) {
  code <- withCallingHandlers(readLines(filename), warning = suppress_final_line_warning)
  # remove spaces
  code <- gsub(" ", "", code)
}

#' return names of all functions in a .R file
#'
#' @param filename location of .R file
#' @return character vector where each entry is the name of a function in the file
get_func_names <- function(filename) {
  # read the .R file, suppressing warnings about incomplete final line
  code <- read_code_unspace(filename)
  # identify lines which define functions
  func_str <- "<-function("
  potential_func_lines <- grep(func_str, code, fixed = TRUE)

  # determine whether each function is an outermost function
  is_outside_func <- function(code, potential_func_line) {
    # assume functions on first line are outside functions
    if(potential_func_line == 1) {
      return(TRUE)
    }
    # paste all code up to name of potential function
    pasted_code <- paste0(code[seq_len(potential_func_line - 1)], collapse = "")
    # potential_func_subline <- strsplit(code, split = func_str, fixed = TRUE)
    # potential_func_subline <- potential_func_subline[1]
    # pasted_code <- paste0(pasted_code, potential_func_subline, collapse = "")
    # count number of open and close curly brackets up to potential function name
    count_characters <- function(pasted_code, char_in) {
      n <- gregexpr(char_in, pasted_code, fixed = TRUE)
      length(n[[1]])
    }
    n_brackets <- vapply(c("{", "}"), function(x) count_characters(pasted_code, x), numeric(1))
    # the functino is an outermost function if the number of open and close brackets is the same
    n_brackets[1] == n_brackets[2]
  }

  func_lines <- potential_func_lines[vapply(potential_func_lines,
                                            function(x) is_outside_func(code, x),
                                            logical(1))]
  # split off the function names
  func_lines <- strsplit(code[func_lines], split = func_str, fixed = TRUE)
  func_lines <- vapply(func_lines, function(x) x[1], character(1))
  func_lines
}

#' a version of apply. if MARGIN = 1, the values in each column of X are passed to FUN
#' as named arguments according to the column name. if MARGIN = 2, the values in
#' each row of X are passed to FUN
#' as named arguments according to the row name
#' @param X see arguments for apply.  If MARGIN = 1, colnames(X) must match the named
#' arguments of FUN.
#' @param MARGIN see arguments for apply
#' @param FUN see arguments for apply.  Must have named arguments
#' @return see apply
apply_named_args <- function(X, MARGIN = 1, FUN, keep_class = TRUE) {
  if(MARGIN == 2 | !keep_class) {
    return(apply(X, MARGIN, function(X) do.call(FUN, as.list(X))))
  }
  args <- formalArgs(FUN)
  X <- c(list(FUN), as.list(X[,args]))
  do.call(Map, X)
}

#' make a default list of random number generator seeds for MCMC
#'
#' @param n_replicates the number of parallel chains to run to assess convergence
#' @param n_temperatures the number of temperatures
#' @param seed_start integer to start seed sequence
#' @return an n_replicates by n_temperatures of integers
make_seed <- function(n_replicates, n_temperatures, seed_start = 1) {
    starting_point_seed <- lapply(seq_len(n_replicates),
                                  function(x) (x - 1) * n_temperatures +
                                      seq_len(n_temperatures) + seed_start - 1)
}

#' find the indices at which a vector changes sign
#'
#' @param x a numeric vector
#' @return a logical vector of length length(x) - 1.  TRUE at element i indicates
#' that x[i] has a different sign from x[i+1].
find_sign_change_ind <- function(x) {
  which(abs(diff(sign(x))) == 2)
}

#' short for \code{vapply(X, FUN, logical(1), ...)}
#'
#' short for \code{vapply(X, FUN, logical(1), ...)}
#'
#' @param X first input to vapply
#' @param FUN second input to vapply
#' @return output of \code{vapply(X, FUN, logical(1), ...)}: logical vector of
#' same length as X
vlapply <- function(X, FUN, ...) {
  vapply(X, FUN, logical(1), ...)
}

#' short for \code{vapply(X, FUN, integer(1), ...)}
#'
#' short for \code{vapply(X, FUN, integer(1), ...)}
#'
#' @param X first input to vapply
#' @param FUN second input to vapply
#' @return output of \code{vapply(X, FUN, integer(1), ...)}: integer vector of
#' same length as X
viapply <- function(X, FUN, ...) {
  vapply(X, FUN, integer(1), ...)
}

#' short for \code{vapply(X, FUN, numeric(1), ...)}
#'
#' short for \code{vapply(X, FUN, numeric(1), ...)}
#'
#' @param X first input to vapply
#' @param FUN second input to vapply
#' @return output of \code{vapply(X, FUN, numeric(1), ...)}: numeric vector of
#' same length as X
vnapply <- function(X, FUN, ...) {
  vapply(X, FUN, numeric(1), ...)
}

#' short for \code{vapply(X, FUN, character(1), ...)}
#'
#' short for \code{vapply(X, FUN, character(1), ...)}
#'
#' @param X first input to vapply
#' @param FUN second input to vapply
#' @return output of \code{vapply(X, FUN, character(1), ...)}: character vector
vcapply <- function(X, FUN, ...) {
  vapply(X, FUN, character(1), ...)
}

#' check if x is within \code{tol} of \code{x}
#'
#' @param x numeric vector of length 1
#' @param tol tolerance
#' @return logical vector of length 1
is_integer_like <- function(x, tol = sqrt(.Machine$double.eps)) {
  is.integer(x) || (is.numeric(x) && abs(x - round(x)) < tol)
}

#' find the file containing a function
#'
#' @param func_name name of function
#' @param dir_name directory in which to search (only searches in that directory
#' + subdirectories one level down)
#' @return the name of the file which has the function
find_func_file <- function(func_name, dir_name = ".") {
  # list current directory + directories one level down
  dirs <- list.dirs(path = dir_name, recursive = FALSE)
  dirs <- c(dirs, dir_name)
  # exclude Rproj.user -- not a relevant directory
  dirs <- dirs[!grepl("Rproj.user", dirs, fixed = TRUE)]
  # list files in those directories
  filenames <- lapply(dirs, function(x) list.files(x, pattern="\\.R$", full.names=TRUE))
  # list functions in those files
  func_names <- lapply(filenames, function(x) lapply(x, get_func_names))
  # find index of file containing target function
  in_file <- lapply(func_names, function(x) vlapply(x, function(y) any(func_name == y)))
  # find directory of that file
  in_dir <- vlapply(in_file, any)
  # return filename

  if(length(filenames[in_dir]) == 0) {
    stop("function not found")
  }
  func_file <- filenames[in_dir][[1]][in_file[in_dir][[1]]]
  # open file in Rstudio
  file.edit(func_file)
  func_file
}

#' loads a character vector of packages
#'
#' loads a character vector of packages
#'
#' @param package_vector: character vector to load
#' @return NULL
load_packages <- function(package_vector){
  lapply(package_vector,function(x) do.call(library,list(x)))
  invisible(NULL)
}

#' sources a character vector of files
#'
#' sources a character vector of files
#'
#' @param source_vector: character vector of files to source
#' @return NULL
source_files <- function(source_vector){
  lapply(source_vector,function(x) source(x))
  invisible(NULL)
}

#' calculate the nth moment from a probability mass function
#'
#' @param x numeric vector. x-values of probability mass function.  Must include
#' all values where the probability mass function is non-zero.
#' @param y numeric vector of the same length as x.  Values of probability mass
#' function at the points in x.
#' @param n integer >= 0.  Calculate nth moment.
#' @param central logical with default value FALSE.  If TRUE, calculate central
#' moment, otherwise calculate raw moment.
#' @return numeric vector of length 1: the nth moment
calc_moment_from_pmf <- function(x, y, n, central = FALSE) {
  stopifnot(n >= 0 && is_integer_like(n))
  stopifnot(length(x) == length(y))
  stopifnot(all.equal(sum(y), 1))
  if(n == 0) {
    return(0)
  }

  if(central) {
    first_moment <- calc_moment_from_pmf(x, y, n = 1, central = FALSE)
    adjustment <- first_moment
  } else {
    adjustment <- 0
  }

  moment <- sum(y * (x - adjustment)^n)
  moment
}

#' calculate percentiles for fitted parameters
#'
#' calculate percentiles for fitted parameters
#'
#' @param chain the MCMC chain: a data frame with n columns, where n is the number
#' of fitted parameters
#' @param prctiles a numeric vector of length m containing the percentiles to be calculated
#' (between 0 and 1)
#' @param par_names_plot parameter names for the output data frame
#' @return a data frame with n rows and m columns containing the percentiles
#' @export
print_prctiles <- function(chain, prctiles = c(.025,.5,.975), par_names_plot){
  prctile_table <- lapply(chain,function(x) quantile(x,prctiles, na.rm = TRUE))
  prctile_table <- t(as.data.frame(prctile_table))
  rownames(prctile_table) <- par_names_plot
  col_names <- as.character(prctiles*100)
  col_names <- trimws(format(col_names, digits = 3))
  col_names <- paste0(col_names,"\\%")
  colnames(prctile_table) <- col_names
  prctile_table
}

#' compile an odin model
#'
#' @param mod_path character vector of length 1: path to ODEs
#' @return compiled model
#' @export
compile_model <- function(mod_path){
  dir_name <- "../model_odin/"
  if(substr(mod_path, 1, nchar(dir_name)) != dir_name) {
    mod_path <- paste0(dir_name, mod_path)
  }
  odin::odin(mod_path, verbose=FALSE) # compile model
}

#' loads a character vector of packages
#'
#' loads a character vector of packages
#'
#' @param package_vector: character vector to load
#' @return NULL
load_packages <- function(package_vector){
  lapply(package_vector,function(x) do.call(library,list(x)))
  invisible(NULL)
}

#' sources a character vector of files
#'
#' sources a character vector of files
#'
#' @param source_vector: character vector of files to source
#' @return NULL
source_files <- function(source_vector){
  lapply(source_vector,function(x) source(x))
  invisible(NULL)
}

#' format a numeric vector of breaks nicely for printing
#'
#' @param breaks numeric vector
#' @param limits numeric vector of length 2. If present, breaks must be within limits
#' @return list with two elements:
#' breaks: numeric vector of same length
#' break_names: character vector of same length
#' order_magnitude: highest order of magnitude of breaks

gen_nice_breaks <- function(breaks, limits){

  # find highest order of magnitude
  order_magnitude <- floor(log10(max(abs(breaks))))
  # divide breaks by that order of magnitude
  breaks <- breaks / (10 ^ order_magnitude)
  # find minimum difference between breaks
  min_diff <- min(abs(diff(breaks)))
  # find order of minimum difference between breaks
  order_min_diff <- floor(log10(min_diff))
  # round breaks accordingly
  breaks <- round(breaks / 10 ^ (order_min_diff - 1)) * 10 ^ (order_min_diff - 1)

  # keep breaks within limits
  if(!missing(limits)){
    # divide limits by that order of magnitude
    limits <- limits / (10 ^ order_magnitude)
    limits_temp <- double(2)

    while(limits_temp[1] == limits_temp[2]){
      limits_temp[1] <- ceiling(limits[1] / 10 ^ order_min_diff) * 10 ^ order_min_diff
      limits_temp[2] <- floor(limits[2] / 10 ^ order_min_diff) * 10 ^ order_min_diff
      order_min_diff <- order_min_diff - 1
    }

    limits <- limits_temp

    replace_breaks_under_limit <- function(breaks, limit){
      if(any(breaks < limit)){
        breaks <- breaks[breaks >= limit]
        breaks <- c(breaks,limit)
      }
      breaks
    }

    breaks <- replace_breaks_under_limit(breaks,limits[1])
    breaks <- -replace_breaks_under_limit(-breaks,-limits[2])

  }

  break_names <- as.character(breaks)
  breaks <- breaks * (10 ^ order_magnitude)
  list(breaks = breaks, break_names = break_names,
       order_magnitude = order_magnitude)
}

#' wrapper function to feed arguments to ggsave
#'
#' @param g either a ggplot object, or a list of ggplot objects of length m,
#' or a nested list of ggplot objects (length n, each of length m)
#' @param filename_fn function which makes a filename by taking n arguments and
#' outputting a character vector of length 1. If g is a ggplot object or a list
#' of ggplot objects of length m, then n = 1
#' @param width numeric vector of length 1. width of figure in cm
#' @param height numeric vector of length 1. height of figure in cm
#' @param filename_args if g is a ggplot object, a charater vector of length 1.
#' if g is a list of ggplot objects of length m, a character vector of length m.
#' if g is a nested list of ggplot objects, a list of length 2, containing a
#' character vector of length n and a character vector of length m respectively.
#' @return NULL
#' @export
ggsave_wrapper <- function(g, filename_fn, width, height, filename_args) {

  ggsave_dims <- function(x, y) ggsave_wch(filename_fn(x), y, width = width, height = height, units = "cm")
  ggsave_dims2 <- function(x1, x2, y) ggsave_wch(filename_fn(x1, x2), y, width = width, height = height, units = "cm")

  if(is.ggplot(g[[1]][[1]])) {
    idx_matrix <- expand.grid(seq_along(g), seq_along(g[[1]]))
    apply(idx_matrix, 1, function(x) ggsave_dims2(filename_args[[1]][x[1]], filename_args[[2]][x[2]], g[[x[1]]][[x[2]]]))
  } else if(is.ggplot(g)){
    ggsave_dims(filename_args[[1]][1], g)
  } else {
    Map(ggsave_dims, filename_args[[1]], g)
  }
  invisible(NULL)
}

#' version of ggsave which suppresses warnings that rows are removed from the data frame,
#' and about transformations introducing infinite values
#' @export
ggsave_wch <- function(...) {
  h <- function(w) {
    if(any(grepl("Removed", w) | grepl("Transformation", w))) {
      invokeRestart("muffleWarning")
    }
  }
  withCallingHandlers(ggsave(...), warning = h)
}

#' normalises a vector so that it sums to 1
#'
#' @param vec numeric vector to be normalised
#' @return normalised vector of same length as vec
normalise <- function(vec) {

  if(all(vec == 0)) {
    return(vec)
  } else if(any(vec < 0)) {
    stop("vector to be normalised must have all non-negative elements")
  }

  vec / sum(vec)

}

get_mod_info <- function(mod_filename) {
    # read the .R file, suppressing warnings about incomplete final line
    mod_filename <- paste0("model_odin/", mod_filename)
    code <- read_code_unspace(mod_filename)

    get_pars <- function() {
      # identify lines which define parameters
      par_str <- "<-user("
      par_lines <- grep(par_str, code, fixed = TRUE)
      code <- code[par_lines] %>%
        lapply(., trimws)
      # exclude commented lines
      code <- code[substr(code, 1, 1) != "#"]

      extract_par_name_from_line <- function(par_str) {
        par_str <- strsplit(par_str, "<-user()", fixed = TRUE)[[1]][1]
        strsplit(par_str, "[", fixed = TRUE)[[1]][1]
      }

      vcapply(code, extract_par_name_from_line) %>% unique
    }

    get_compartments <- function() {
      compartment_str <- "initial("
      compartment_lines <- grep(compartment_str, code, fixed = TRUE)
      code <- code[compartment_lines] %>%
        lapply(., trimws)
      # exclude commented lines
      code <- code[substr(code, 1, 1) != "#"]

      extract_compartment_name_from_line <- function(compartment_str) {
        compartment_str <- strsplit(compartment_str, "initial(", fixed = TRUE)[[1]][2]
        compartment_str <- strsplit(compartment_str, "[", fixed = TRUE)[[1]][1]
        strsplit(compartment_str, ")", fixed = TRUE)[[1]][1]
      }

      vcapply(code, extract_compartment_name_from_line) %>% unique
    }

    list(pars = get_pars(), compartments = get_compartments())

}

order_cols <- function(my_df) {
  my_df[,order(colnames(my_df))]
}

is.error <- function(x) inherits(x, "try-error")

expect_equal_tol <- function(x1, x2, FUN = identity, my_message = NULL, ...)  {
  a <- try(testthat::expect_equal(FUN(x1), FUN(x2), tol = 1e-5, ...))

  if(is.error(a)) {
    message(my_message)
    return(list(FUN(x1), FUN(x2)))
  }
  invisible(a)
}

#' sum time series of odin solution across compartments defined as vectors
#'
#' @param df data frame containing odin solution
#' @param compartment_names character vector. Names of compartments to sum over
#' @return a data frame with one column for each element of compartment_names,
#' containing the solution summed over that vector compartment
sum_stages <- function(df, compartment_names){
  extract_and_sum_stages <- function(df, compartment_name) {
    extracted_df <- df[,grep(compartment_name, colnames(df)), drop = FALSE]
    if(ncol(extracted_df) == 0) {
      double(nrow(df))
    } else {
      rowSums(extracted_df)
    }
  }

  summed_df <- vapply(compartment_names, extract_and_sum_stages, numeric(nrow(df)), df = df)
  colnames(summed_df) <- compartment_names
  summed_df
}

grepl_any <- function(patterns, x, ...) {
  temp <- vapply(patterns, function(y) grepl(y, x, ...), logical(length(x)))
  apply(temp, 1, any)
}

grepl_all <- function(patterns, x, ...) {
  temp <- vapply(patterns, function(y) grepl(y, x, ...), logical(length(x)))
  apply(temp, 1, all)
}

## a useful function: rev() for strings
strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")


make_clean_dir <- function(my_dir) {
  system(paste0("rm -r ", my_dir))
  system(paste0("mkdir ", my_dir))
}

#' Emulates the ggplot colour palette.
#'
#' \code{gg_color_hue} emulates the ggplot colour palette.
#'
#' @param n numeric vector of length 1. The number of colours required.
#' @return a character vector of length \code{n}. Each element is a hex colour code.
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' format a p-value for plotting
#'
#' @param p-value numeric vector of length 1
#' @return character vector of length 1
make_p_value_str <- function(p_value) {
  if(is.na(p_value)) {
    ""
  } else if(p_value < 0.001) {
    "< 0.001"
  } else {
    as.character(round(p_value, 3))
  }
}

#' probabilistic rounding
#'
#' @param x scalar
#' @return x rounded up with probability x - floor(x)
prob_round <- function(x) {
  floor(x) + rbinom(1, 1, x - floor(x))
}

#' rounds a numeric vector probabilistically such that its sum (an integer) is preserved
#'
#' @param vec numeric vector which sums to an integer
#' @return numeric vector of the same length as vec: rounded vector
round_preserve_sum <- function(vec) {

  # if vector already integers, do nothing
  # old implementation -- all.equal(vec - round(vec)) -- fails for large values in vec
  if(isTRUE(all.equal(vec - round(vec), numeric(length(vec))))) {
    return(vec)
  }

  # if vector doesn't sum to integer,
  # round down so it does sum to an integer (ensuring positivity)
  sum_vec <- sum(vec)
  rounded_sum <- round(sum_vec)
  fractional_part <- rounded_sum - sum_vec
  if(!isTRUE(all.equal(rounded_sum - sum_vec, 0))) {
    larger_idx <- which(vec > fractional_part)
    # subtract fractional_part from one of those groups, chosen randomly
    subtract_idx <- sample(larger_idx, 1)
    vec[subtract_idx] <- vec[subtract_idx] - fractional_part
  }

  n_round_up <- rounded_sum - sum(floor(vec))
  round_up_idx <- sample.int(length(vec), size = n_round_up, replace = FALSE,
                             prob = vec - floor(vec))
  vec_out <- floor(vec)
  vec_out[round_up_idx] <- ceiling(vec[round_up_idx])
  vec_out
}

simple_solve <- function(gen, pars, solving_time) {
  mod <- gen$new(user = pars)
  sol <- solve_ODE_error_handling(mod, solving_time) %>%
    as_tibble()
  sol
}

#' make dir for Rmd outputs
#' @return string: directory name
#' @importFrom dplyr %>%
#' @export
make_dir_Rmd <- function(input_filename) {
  dir_name <- input_filename %>%
    sub(".Rmd", "_files1/", ., fixed = TRUE) %>%
    paste0("output/", .)
  if(!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
  dir_name
  
}
