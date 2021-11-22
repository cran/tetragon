#' tetragon
#'
#' @param df A data frame with time features on columns.
#' @param seq_len Positive integer. Time-step number of the projected sequence.
#' @param deriv Integer vector. Number of differentiation operations to perform for each original time feature. 0 = no change; 1: one diff; 2: two diff.
#' @param ci Confidence interval. Default: 0.8.
#' @param method String. Distance method for calculating distance matrix among sequences. Possibile options are: "euclidean", "manhattan", "chebyshev", "gower", "lorentzian", "jaccard", "dice", "squared_euclidean", "divergence", "clark", "avg". For further information, please make reference to philentropy package. Default: NULL (random selection among all possible options).
#' @param distr String. DIstribution used to expand the distance matrix. Possible options are: "norm", "cauchy", "logis", "t", "exp". Default: NULL (random selection among all possible options).
#' @param n_windows Positive integer. Number of validation tests to measure/sample error. Default: 10.
#' @param n_sample Positive integer. Number of samples for grid or random search. Default: 30.
#' @param search String. Two option available: "grid", "random". Default: "random".
#' @param fixed Logical. Setting to TRUE, calculate a single model (if the variables passed to the functions are completed). Default: FALSE.
#' @param dates Date. Vector with dates for time features.
#' @param seed Positive integer. Random seed. Default: 42.
#'
#' @author Giancarlo Vercellino \email{giancarlo.vercellino@gmail.com}
#'
#' @return This function returns a list including:
#' \itemize{
#' \item exploration: list of all not-null models, complete with predictions, test metrics, prediction stats and plot
#' \item history: a table with the sampled models, hyper-parameters, validation errors, weighted average rank
#' \item best: results for the best model in term of weighted average rank, including:
#' \itemize{
#' \item wt_avg_best: hyper-parameters of the best model selected through grid/random search
#' \item predictions: min, max, q25, q50, q75, quantiles at selected ci, mean, sd, mode, kurtosis, skewness for each time feature
#' \item testing_errors: testing errors for each time feature (ME, MAE, MSE, MPE, MAPE, sCE, MAPE, sMAE, sMSE, MASE, RMSSE, rMAE, rRMSE, rAME). For further information, make reference to greybox package.
#' \item pred_stats: for each predicted time feature, IQR to range, KL-divergence, risk ratio, upside probability, averaged across prediction time-points and at the terminal points.
#' }
#' \item time_log
#' }
#'
#' @export
#'
#' @import purrr
#' @import tictoc
#' @importFrom readr parse_number
#' @importFrom lubridate seconds_to_period is.Date as.duration
#' @importFrom stats weighted.mean
#' @importFrom imputeTS na_kalman


#'@examples
#'\donttest{
#'tetragon(covid_in_europe[, c(2, 4)], seq_len = 40, n_sample = 2, deriv = c(1, 2))
#'}

tetragon <- function(df, seq_len = NULL, deriv = NULL, ci = 0.8, method = NULL, distr = NULL, n_windows = 10, n_sample = 30, search = "random", fixed = FALSE, dates = NULL, seed = 42)
{

  hood <- function(df, seq_len, deriv, ci, method, distr, n_windows, dates)
  {

    n_ts <- nrow(df)
    feat_names <- colnames(df)

    test_index <- unique(round(seq(2 * seq_len, n_ts, length.out = n_windows)))
    if(length(test_index) < n_windows){message("not enough data, testing on ", length(test_index), " possible windows")}

    results <- map(test_index, ~ tryCatch(engine(df[1:.x,, drop = F], seq_len, deriv, ci, method, distr, measure_error = TRUE, n_sim = 30, dates), error = function(e) NA))
    not_na <- !is.na(results)
    results <- results[not_na]
    window_errors <- map(results, ~ .x$errors)
    testing_errors <- Reduce("+", window_errors)/length(test_index)
    model <- engine(df, seq_len, deriv, ci, method, distr, measure_error = FALSE, n_sim = 100, dates)

    outcome <- list(model = model, testing_errors = testing_errors)

    return(outcome)
  }


  set.seed(seed)

  tic.clearlog()
  tic("time")

  if(!is.data.frame(df)){stop("time features must be in dataframe format")}
  if(anyNA(df)){df <- as.data.frame(na_kalman(df))}

  n_ts <- nrow(df)
  n_feat <- ncol(df)
  max_limit <- floor(n_ts/3)

  if(!is.null(method)){method <- match.arg(arg = method, choices = c("euclidean", "manhattan", "chebyshev", "gower", "lorentzian", "jaccard", "dice", "squared_euclidean", "divergence", "clark", "avg"), several.ok = TRUE)}
  if(!is.null(distr)){distr <- match.arg(arg = distr, choices = c("norm", "t", "cauchy", "exp", "logis"), several.ok = TRUE)}

  if(!is.null(seq_len) && any(seq_len < 2 | seq_len > max_limit)){seq_len[seq_len < 2] <- 2; seq_len[seq_len > max_limit] <- max_limit; message("setting seq_len value between 2 and 1/3 time-feature length")}
  if(!is.null(deriv) && any(deriv > 2 | deriv < 0)){deriv[deriv > 2] <- 2; deriv[deriv < 0] <- 0; message("setting deriv value between 0 and 2")}

  if(fixed == TRUE && (is.null(seq_len) | is.null(deriv) | is.null(method)) && length(deriv) == n_feat){stop("fixed flag requires non-null values for seq_len, method, distr, and non-null values for deriv with length equal to the number of time features")}

  if(fixed == FALSE)
  {
    if(is.null(seq_len)){sl_sample_field <- 2:floor(sqrt(n_ts))} else {sl_sample_field <- seq_len}
    if(is.null(deriv)){d_sample_field <- 0:2} else {d_sample_field <- deriv}

    if(is.null(method)){m_sample_field <- c("euclidean", "manhattan", "chebyshev", "gower", "lorentzian",
                                            "jaccard", "dice", "squared_euclidean", "divergence", "clark", "avg")} else {m_sample_field <- method} ###

    if(is.null(distr)){dsr_sample_field <- c("norm", "t", "cauchy", "exp", "logis")} else {dsr_sample_field <- distr}

    if(search == "random")
    {
      if(length(sl_sample_field) == 1){sl_sample <- rep(sl_sample_field, n_sample)} else {sl_sample <- sample(sl_sample_field, n_sample, replace = TRUE)}
      if(length(d_sample_field) == 1){d_sample <- replicate(n_sample, rep(d_sample_field, n_feat), simplify = F)}
      if(length(d_sample_field) > 1){d_sample <- replicate(n_sample, sample(d_sample_field, n_feat, replace = TRUE), simplify = F)}
      #if(length(d_sample_field) > 1 & length(d_sample_field) != n_feat){d_sample <- replicate(n_sample, sample(d_sample_field, n_feat, replace = TRUE), simplify = F)}
      #if(length(d_sample_field) > 1 & length(d_sample_field) == n_feat){d_sample <- replicate(n_sample, d_sample_field, simplify = F)}
      if(length(m_sample_field) == 1){m_sample <- replicate(n_sample, m_sample_field, simplify = F)}
      if(length(m_sample_field) > 1){m_sample <- replicate(n_sample, sample(m_sample_field, 1, replace = TRUE), simplify = F)}
      if(length(dsr_sample_field) == 1){dsr_sample <- replicate(n_sample, dsr_sample_field, simplify = F)}
      if(length(dsr_sample_field) > 1){dsr_sample <- replicate(n_sample, sample(dsr_sample_field, 1, replace = TRUE), simplify = F)}
    }

    if(search == "grid")
    {
      d_perm <- narray::split(expand.grid(replicate(n_feat, d_sample_field, simplify = F)), along = 1)
      grid_list <- narray::split(expand.grid(sl_sample_field, d_perm, m_sample_field, dsr_sample_field), 2)
      sl_sample <- grid_list[[1]]
      d_sample <- map(grid_list[[2]], ~ unlist(.x))
      m_sample <- grid_list[[3]]
      dsr_sample <- grid_list[[4]]
    }

    sl_sample <- pmap(list(sl_sample, d_sample), ~ (..1 > (n_ts - max(..2) - n_windows)) * floor((n_ts - max(..2) - n_windows)) +  (..1 <= (n_ts - max(..2) - n_windows)) * ..1)
    sl_sample <- map2(sl_sample, d_sample, ~ ((.x - max(.y)) == 0) * (.x + 2) + ((.x - max(.y)) == 1) * (.x + 1) + ((.x - max(.y)) > 1) * .x)

    if(length(seq_len)==1 && any(sl_sample != seq_len)){message("fixing seq_len for available data")}

    hyper_params <- list(sl_sample, d_sample, m_sample, dsr_sample)
  }

  if(fixed == TRUE){hyper_params <- list(sl_sample = list(seq_len), d_sample = list(deriv), m_sample = list(method), dsr_sample = list(distr))}

  exploration <- suppressMessages(pmap(hyper_params, ~ tryCatch(hood(df, seq_len = ..1, deriv = ..2, ci, method = ..3, distr = ..4, n_windows, dates), error = function(e) NA)))

  not_na <- !is.na(exploration)
  exploration <- exploration[not_na]

  if(fixed == FALSE)
  {
    aggr_errors <- t(as.data.frame(map(exploration, ~ tryCatch(apply(.x$testing_errors, 2, function(x) x[which.max(abs(x))]), error = function(e) rep(NA, 13)))))
    colnames(aggr_errors) <- paste0("max_", colnames(aggr_errors))
    weights <- apply(aggr_errors, 2, function(x) {abs(sd(x[is.finite(x)], na.rm = TRUE)/mean(x[is.finite(x)], na.rm = TRUE))})
    finite_w <- is.finite(weights)
    history <- as.data.frame(Reduce(cbind, list(sl_sample[not_na], d_sample[not_na])))
    history$method <- unlist(m_sample[not_na])
    history$distr <- unlist(dsr_sample[not_na])
    colnames(history) <- c("seq_len", "deriv", "method", "distr")
    history <- data.frame(history, aggr_errors)
    rownames(history) <- NULL

    history$wgt_avg_rank <- round(apply(apply(abs(aggr_errors[, finite_w, drop = FALSE]), 2, rank), 1, weighted.mean, w = weights[finite_w]), 2)
    history <- history[order(history$wgt_avg_rank),]

    best_index <- as.numeric(rownames(history[1,]))
    best <- list(wt_avg_best = history[1, c(1:4, 18)], predictions = exploration[[best_index]][[1]]$predictions, testing_errors = exploration[[best_index]]$testing_errors, pred_stats = exploration[[best_index]][[1]]$pred_stats, plot = exploration[[best_index]][[1]]$plot)
  }

  if(fixed == TRUE)
  {
    history <- NULL
    best <- exploration[[1]]
  }

  toc(log = TRUE)
  time_log <- seconds_to_period(round(parse_number(unlist(tic.log())), 0))

  outcome <- list(exploration = exploration, history = history, best = best, time_log = time_log)

  return(outcome)
}
