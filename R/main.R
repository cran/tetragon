#' tetragon
#'
#' @param df A data frame with time features as columns. They could be continuous variables or not.
#' @param seq_len Positive integer. Time-step number of the projected sequence. Default: NULL (random selection between maximum boundaries).
#' @param ci Confidence interval. Default: 0.8.
#' @param method String. Distance method for calculating distance matrix among sequences. Options are: "euclidean", "manhattan", "maximum", "minkowski". Default: NULL (random selection among all possible options).
#' @param distr String. Distribution used to expand the distance matrix. Options are: "norm", "logis", "t", "exp", "chisq". Default: NULL (random selection among all possible options).
#' @param n_windows Positive integer. Number of validation tests to measure/sample error. Default: 3 (but a larger value is strongly suggested to really understand your accuracy).
#' @param n_sample Positive integer. Number of samples for random search. Default: 30.
#' @param dates Date. Vector with dates for time features.
#' @param error_scale String. Scale for the scaled error metrics (only for continuous variables). Two options: "naive" (average of naive one-step absolute error for the historical series) or "deviation" (standard error of the historical series). Default: "naive".
#' @param error_benchmark String. Benchmark for the relative error metrics (only for continuous variables). Two options: "naive" (sequential extension of last value) or "average" (mean value of true sequence). Default: "naive".
#' @param seed Positive integer. Random seed. Default: 42.
#'
#' @author Giancarlo Vercellino \email{giancarlo.vercellino@gmail.com}
#'
#' @return This function returns a list including:
#' \itemize{
#' \item exploration: list of all explored models, complete with predictions, testing metrics and plots
#' \item history: a table with the sampled models, hyper-parameters, validation errors
#' \item best: results for the best model including:
#' \itemize{
#' \item predictions: min, max, q25, q50, q75, quantiles at selected ci, and a bunch of specific measures for each point fo predicted sequences
#' \item testing_errors: testing errors for one-step and sequence for each ts feature
#' \item plots: confidence interval plot for each time feature
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
#' @importFrom stats weighted.mean ecdf lm rnorm na.omit
#' @importFrom imputeTS na_kalman
#' @importFrom narray split
#' @importFrom utils head tail
#' @importFrom modeest mlv1 mfv1
#' @importFrom moments skewness kurtosis
#' @importFrom stats dcauchy dexp dlogis dnorm dt pcauchy pexp plogis pnorm pt median sd density qcauchy qlogis qnorm qt quantile runif cor
#' @import abind
#' @import philentropy
#' @import abind
#' @import ggplot2
#' @import greybox
#' @import stringr
#' @importFrom dqrng dqsample
#' @importFrom scales number
#' @importFrom entropy entropy
#' @importFrom Rfast Dist


#'@examples
#'\donttest{
#'tetragon(covid_in_europe[, c(2, 4)], seq_len = 40, n_sample = 2)
#'}

tetragon <- function(df, seq_len = NULL, ci = 0.8, method = NULL, distr = NULL, n_windows = 3, n_sample = 30, dates = NULL, error_scale = "naive", error_benchmark = "naive", seed = 42)
{
  set.seed(seed)

  tic.clearlog()
  tic("time")

  if(!is.data.frame(df)){stop("time features must be in dataframe format")}

  n_ts <- nrow(df)
  n_feat <- ncol(df)
  feat_names <- colnames(df)

  feat_n_class <- NULL
  feat_level_names <- NULL

  all_classes <- all(map_lgl(df, ~ is.factor(.x) | is.character(.x)))
  all_numbers <- all(map_lgl(df, ~ is.numeric(.x) | is.integer(.x)))

  if(!all_classes & !all_numbers){stop("only numbers or factors, but not both")}

  event_class <- all_classes & !all_numbers
  if(anyNA(df) & event_class == FALSE){df <- as.data.frame(na_kalman(df))}

  if(event_class == TRUE)
  {
    df <- as.data.frame(map(df, ~ as.factor(.x)))
    feat_level_names <- map(df, ~ levels(.x))
    feat_n_class <- map_dbl(feat_level_names, ~ length(.x))
    df <- as.data.frame(data.matrix(df)-1)
    deriv <- rep(0, n_feat)
  }

  if(event_class == FALSE){deriv <- map_dbl(df, ~ best_deriv(.x))}

  max_limit <- floor((n_ts - max(deriv))/(4 + n_windows + 1))
  if(max_limit < 3){stop("not enough data for the validation windows")}

  if(all(c(length(seq_len) == 1, length(method) == n_feat, length(distr) == n_feat, n_sample == 1))){fixed <- TRUE} else {fixed <- FALSE}
  if(!is.null(method)){method <- match.arg(arg = method, choices = c("euclidean", "manhattan", "maximum", "minkowski", "kullback_leibler"), several.ok = TRUE)}
  if(!is.null(distr)){distr <- match.arg(arg = distr, choices = c("norm", "t", "exp", "logis", "chisq"), several.ok = TRUE)}
  if(!is.null(seq_len) && any(seq_len < 3 | seq_len > max_limit)){suppressWarnings(seq_len[seq_len < 3] <- 3); suppressWarnings(seq_len[seq_len > max_limit] <- max_limit); message(paste0("setting seq_len within boundaries (min 3, max ", max_limit,")"))}

  ###OLD: c("euclidean", "manhattan", "chebyshev", "gower", "lorentzian", "jaccard", "dice", "squared_euclidean", "divergence", "clark", "avg")

  if(fixed == FALSE)
  {
    sl_sample <- sampler(seq_len, n_samp = n_sample, range = c(3, max_limit), integer = TRUE)
    m_sample <- sampler(method, n_samp = n_sample, range = c("euclidean", "manhattan", "maximum", "minkowski"), multi = n_feat)
    dsr_sample <- sampler(distr, n_samp = n_sample, range = c("norm", "t", "exp", "logis", "chisq"), multi = n_feat)

    if(event_class == TRUE){sl_sample <- map2(sl_sample, replicate(n_sample, deriv, simplify = FALSE), ~ ((.x - max(.y)) == 0) * (.x + 2) + ((.x - max(.y)) == 1) * (.x + 1) + ((.x - max(.y)) > 1) * .x)
    if(length(seq_len)==1 && any(sl_sample != seq_len)){message("fixing seq_len for available data")}}

    hyper_params <- list(sl_sample, m_sample, dsr_sample)
  }

  if(fixed == TRUE){hyper_params <- list(sl_sample = list(seq_len), m_sample = list(method), dsr_sample = list(distr))}

  exploration <- lapply(1:n_feat, function(ft) pmap(hyper_params, ~ tryCatch(windower(df[, ft, drop = TRUE], seq_len = ..1, deriv[ft], ci, method = ..2[ft], distr = ..3[ft], n_windows, dates, error_scale, error_benchmark, feat_names[ft], feat_n_class[ft], feat_level_names[[ft]], seed), error = function(e) NA)))

  #failures <- map_lgl(exploration, ~ anyNA(.x))
  #failed <- hyper_params[failures]
  #return(failed)

  exploration <- transpose(exploration)

  testing_errors <- map(map_depth(exploration, 2, ~ .x$testing_errors), ~ Reduce(rbind, .x))

  if(fixed == FALSE)
  {
    aggr_errors <- t(as.data.frame(map(testing_errors, ~ apply(.x, 2, function(x){mean(x[is.finite(x)], na.omit = TRUE)}))))
    history <- data.frame(seq_len = unlist(sl_sample))
    history$method <- m_sample
    history$distr <- dsr_sample
    history <- data.frame(history, aggr_errors)
    rownames(history) <- NULL

    if(event_class == FALSE){history <- ranker(history, focus = -c(1, 2, 3), inverse = NULL, absolute = c("me", "mpe", "sce"), reverse = FALSE)}
    if(event_class == TRUE){history <- ranker(history, focus = -c(1, 2, 3), inverse = NULL, absolute = NULL, reverse = FALSE)}

    best_index <- as.numeric(rownames(history[1,]))
  }

  if(fixed == TRUE)
  {
    history <- NULL
    best_index <- 1
  }


  pred_scores <- map_depth(exploration, 2, ~ .x$pred_scores)[[best_index]]
  predictions <- map(exploration[[best_index]], ~ .x$model$predictions)
  predictions <- map2(predictions, pred_scores, ~ cbind(.x, pred_scores = .y))

  names(predictions) <- feat_names
  plots <- map(exploration[[best_index]], ~ .x$model$plot)
  names(plots) <- feat_names

  best <- list(predictions = predictions, testing_errors = testing_errors[[best_index]], plots = plots)

  toc(log = TRUE)
  time_log <- seconds_to_period(round(parse_number(unlist(tic.log())), 0))

  outcome <- list(exploration = exploration, history = history, best = best, time_log = time_log)

  return(outcome)
}


###SUPPORT FUNCTIONS
windower <- function(ts, seq_len, deriv, ci, method, distr, n_windows, dates, error_scale, error_benchmark, feat_name, n_class, level_names, seed)
{
  n_ts <- length(ts)
  test_index <- unique(round(seq(seq_len * 4, n_ts, length.out = n_windows + 1)))
  n_test <- length(test_index)
  if((n_test - 1) < n_windows){message("not enough data, testing on ", (n_test - 1), " possible windows")}

  results <- map(1:n_test, ~ engine(ts[1:test_index[.x]], seq_len, deriv, ci, method, distr, measure_error = ifelse(.x < n_test, TRUE, FALSE), n_sim = 100, dates, error_scale, error_benchmark, feat_name, n_class, level_names, seed))

  window_errors <- head(map(results, ~ .x$errors), -1)
  testing_errors <- Reduce("+", window_errors)/length(window_errors)
  window_pred_scores <- head(map(results, ~ .x$pred_scores), -1)
  pred_scores <- round(Reduce("+", window_pred_scores)/length(window_pred_scores), 4)
  model <- results[[n_test]]

  outcome <- list(model = model, testing_errors = testing_errors, pred_scores = pred_scores)

  return(outcome)
}


####
engine <- function(ts, seq_len, deriv, ci, method, distr, measure_error, n_sim, dates, error_scale, error_benchmark, feat_name, n_class, level_names, seed)
{

  actual <- NULL
  if(measure_error == TRUE){actual <- tail(ts, seq_len); df <- head(ts, - seq_len)}
  n_ts <- length(ts)

  diff_model <- recursive_diff(ts, deriv)
  diff_vector <- diff_model$vector

  segmented <- as.data.frame(smart_reframer(diff_vector, seq_len, seq_len))
  dist <- Dist(segmented, method = method, p = 3)
  minmax <- function(x){(x - min(x))/(max(x) - min(x))}
  dist <- apply(dist, 2, minmax)

  if(distr == "norm")
  {
    mean_weights <- apply(dist, 1, mean, na.rm = TRUE)
    sd_weights <- apply(dist, 1, sd, na.rm = TRUE)
    weights <- replicate(n_sim, map2_dbl(mean_weights, sd_weights, ~ trunc_norm(1, a = 0, b = Inf, mean = .x, sd = .y)), simplify = FALSE)
  }

  if(distr == "t")
  {
    df_weights <- ncol(dist) - 1
    ncp_weights <- apply(dist, 1, mean, na.rm = TRUE)
    weights <- replicate(n_sim, map2_dbl(df_weights, ncp_weights, ~ trunc_t(1, a = 0, b = Inf, df = .x, ncp = .y)), simplify = FALSE)
  }

  if(distr == "cauchy")
  {
    location_weights <- apply(dist, 1, mean, na.rm = TRUE)
    scale_weights <- apply(dist, 1, sd, na.rm = TRUE)
    scale_weights[scale_weights <= 0] <- mean(scale_weights[scale_weights> 0])
    weights <- replicate(n_sim, map2_dbl(location_weights, scale_weights, ~ trunc_cauchy(1, a = 0, b = Inf, location = .x, scale = .y)), simplify = FALSE)
  }

  if(distr == "exp")
  {
    rate_weights <- 1/apply(dist, 1, sd, na.rm = TRUE)
    weights <- replicate(n_sim, map_dbl(rate_weights, ~ suppressWarnings(rexp(1, rate = .x))), simplify = FALSE)
  }

  if(distr == "pois")
  {
    lambda_weights <- apply(dist, 1, mean, na.rm = TRUE)
    weights <- replicate(n_sim, map_dbl(lambda_weights, ~ suppressWarnings(rpois(1, lambda = .x))), simplify = FALSE)
  }

  if(distr == "chisq")
  {
    df_weights <- ncol(dist) - 1
    ncp_weights <- apply(dist, 1, mean, na.rm = TRUE)
    weights <- replicate(n_sim, map2_dbl(df_weights, ncp_weights, ~ suppressWarnings(rchisq(1, df =.x, ncp = .y))), simplify = FALSE)
  }

  if(distr == "logis")
  {
    location_weights <- apply(dist, 1, mean, na.rm = TRUE)
    scale_weights <- apply(dist, 1, sd, na.rm = TRUE)
    weights <- replicate(n_sim, map2_dbl(location_weights, scale_weights, ~ trunc_logis(1, a = 0, b = Inf, location = .x, scale = .y)), simplify = FALSE)
  }

  if(distr == "empirical")
  {
    weights <- replicate(n_sim, matrix_expansion(dist), simplify = FALSE)
  }

  weights <- map(weights, ~ {.x[.x == Inf] <- suppressWarnings(max(.x[is.finite(.x)], na.rm = TRUE)); return(.x)})
  weights <- map(weights, ~ {.x[.x == -Inf] <- suppressWarnings(min(.x[is.finite(.x)], na.rm = TRUE)); return(.x)})
  weights <- map(weights, ~ {.x[is.nan(.x)] <- suppressWarnings(median(.x[is.finite(.x)], na.rm = TRUE)); return(.x)})

  inv_weights <- map(weights, ~ 1/.x)

  if(any(map_lgl(inv_weights, ~ anyNA(.x)))){print(seq_len); print(method); print(distr)}

  if(is.numeric(n_class)){inv_weights <- map(inv_weights, ~ round(minmax_transf(.x, min_r = 1, max_r = 10)))}

  #sample_pred_seeds <- t(mapply(function(wts) apply(segmented * inv_weights[[wts]], 2, sum, na.rm = TRUE)/sum(inv_weights[[wts]], na.rm = TRUE), wts = 1:length(inv_weights)))
  if(is.null(n_class))
  {
    sample_pred_seeds <- t(mapply(function(wts) apply(segmented * inv_weights[[wts]], 2, sum, na.rm = TRUE)/sum(inv_weights[[wts]], na.rm = TRUE), wts = 1:length(inv_weights)))
    heads <- diff_model$tail_value
    pred_seeds <- t(apply(sample_pred_seeds, 1, function(x) tail(invdiff(x, heads), seq_len)))
  }

  if(is.numeric(n_class)){pred_seeds <- t(mapply(function(wts) apply(expand(segmented, inv_weights[[wts]]), 2, function(x) most_frequent(unlist(x))), wts = 1:length(inv_weights)))}

  errors <- NULL
  pred_scores <- NULL
  if(measure_error == TRUE)
  {
    pred_seeds <- doxa_filter(ts = head(ts, - seq_len), pred_seeds, n_class)
    pred_scores <- prediction_score(pred_seeds, actual)
    errors <- apply(Reduce(rbind, map(split(pred_seeds, along = 1), ~ my_metrics(holdout = actual, forecast = .x, actual = head(ts, - seq_len), error_scale, error_benchmark, n_class))), 2, mean, na.rm = TRUE)
    errors <- t(as.data.frame(errors))
    rownames(errors) <- feat_name
  }

  plot <- NULL
  predictions <- NULL
  if(measure_error == FALSE)
  {
    predictions <- quantile_prediction(raw_pred = pred_seeds, seed_pred = NULL, raw_error = NULL, holdout = NULL, tail_value = NULL, ts, ci, doxa = TRUE,
    error_scale, error_benchmark, n_class, level_names, seed)$quant_pred

    if(is.Date(dates)){new_dates<- seq.Date(tail(dates, 1), tail(dates, 1) + seq_len * mean(diff(dates)), length.out = seq_len); rownames(predictions) <- as.character(new_dates)} else {new_dates <- 1:seq_len}
    if(is.Date(dates)){x_hist <- dates; x_forcat <- new_dates} else {x_hist <- 1:n_ts; x_forcat <- (n_ts + 1):(n_ts + seq_len); rownames(predictions) <- paste0("t", 1:seq_len)}
    x_lab <- paste0("Forecasting Horizon for sequence n = ", seq_len)
    y_lab <- paste0("Forecasting Values for ", feat_name)

    if(is.numeric(n_class)){ts <- level_names[ts + 1]}
    lower_b <- paste0((1-ci)/2 * 100, "%")
    upper_b <- paste0((ci+(1-ci)/2) * 100, "%")
    plot <- ts_graph(x_hist = x_hist, y_hist = ts, x_forcat = x_forcat, y_forcat = predictions[, "50%"], lower = predictions[, lower_b], upper = predictions[, upper_b], label_x = x_lab, label_y = y_lab)

   }

  outcome <- list(predictions = predictions, errors = errors, pred_scores = pred_scores, plot = plot)
  return(outcome)
}


######
smart_reframer <- function(ts, seq_len, stride)
{
  n_length <- length(ts)
  if(seq_len > n_length | stride > n_length){stop("vector too short for sequence length or stride")}
  if(n_length%%seq_len > 0){ts <- tail(ts, - (n_length%%seq_len))}
  n_length <- length(ts)
  idx <- base::seq(from = 1, to = (n_length - seq_len + 1), by = 1)
  reframed <- t(sapply(idx, function(x) ts[x:(x+seq_len-1)]))
  if(seq_len == 1){reframed <- t(reframed)}
  idx <- rev(base::seq(nrow(reframed), 1, - stride))
  reframed <- reframed[idx,,drop = FALSE]
  colnames(reframed) <- paste0("t", 1:seq_len)
  return(reframed)
}

###
recursive_diff <- function(vector, deriv)
{
  vector <- unlist(vector)
  head_value <- vector("numeric", deriv)
  tail_value <- vector("numeric", deriv)
  if(deriv==0){head_value = NULL; tail_value = NULL}
  if(deriv > 0){for(i in 1:deriv){head_value[i] <- head(vector, 1); tail_value[i] <- tail(vector, 1); vector <- diff(vector)}}
  outcome <- list(vector = vector, head_value = head_value, tail_value = tail_value)
  return(outcome)
}


###
invdiff <- function(vector, heads, add = FALSE)
{
  vector <- unlist(vector)
  if(is.null(heads)){return(vector)}
  for(d in length(heads):1){vector <- cumsum(c(heads[d], vector))}
  if(add == FALSE){return(vector[-c(1:length(heads))])} else {return(vector)}
}

###
matrix_expansion <- function(dist_mat)
{
  n_row <- nrow(dist_mat)
  n_col <- ncol(dist_mat)

  oth_cols <- mapply(function(r) dqsample(dist_mat[r, -r], size = 1), r = 1:n_row)
  pred <- sapply(1:n_row, function(x) dqsample(c(dist_mat[x, -x], oth_cols[-x]), size = 1))

  return(pred)
}

###
my_metrics <- function(holdout, forecast, actuals, error_scale = "naive", error_benchmark = "naive", n_class = NULL)
{
  if(is.null(n_class))
  {
    scale <- switch(error_scale, "deviation" = sd(actuals), "naive" = mean(abs(diff(actuals))))
    benchmark <- switch(error_benchmark, "average" = rep(mean(forecast), length(forecast)), "naive" = rep(tail(actuals, 1), length(forecast)))
    me <- ME(holdout, forecast, na.rm = TRUE)
    mae <- MAE(holdout, forecast, na.rm = TRUE)
    mse <- MSE(holdout, forecast, na.rm = TRUE)
    rmsse <- RMSSE(holdout, forecast, scale, na.rm = TRUE)
    mre <- MRE(holdout, forecast, na.rm = TRUE)
    mpe <- MPE(holdout, forecast, na.rm = TRUE)
    mape <- MAPE(holdout, forecast, na.rm = TRUE)
    rmae <- rMAE(holdout, forecast, benchmark, na.rm = TRUE)
    rrmse <- rRMSE(holdout, forecast, benchmark, na.rm = TRUE)
    rame <- rAME(holdout, forecast, benchmark, na.rm = TRUE)
    mase <- MASE(holdout, forecast, scale, na.rm = TRUE)
    smse <- sMSE(holdout, forecast, scale, na.rm = TRUE)
    sce <- sCE(holdout, forecast, scale, na.rm = TRUE)
    gmrae <- GMRAE(holdout, forecast, benchmark, na.rm = TRUE)
    out <- round(c(me = me, mae = mae, mse = mse, rmsse = rmsse, mpe = mpe, mape = mape, rmae = rmae, rrmse = rrmse, rame = rame, mase = mase, smse = smse, sce = sce, gmrae = gmrae), 3)
  }

  if(is.numeric(n_class))
  {
    avg <- suppressMessages(distance(rbind(holdout, forecast), method = "avg"))
    tanimoto <- suppressMessages(distance(rbind(holdout, forecast), method = "tanimoto"))
    hassebrook <- 1 - suppressMessages(distance(rbind(holdout, forecast), method = "hassebrook"))
    jaccard <- suppressMessages(distance(rbind(holdout, forecast), method = "jaccard"))
    taneja <- suppressMessages(distance(rbind(holdout, forecast), method = "taneja"))
    canberra <- suppressMessages(distance(rbind(holdout, forecast), method = "canberra"))
    gower <- suppressMessages(distance(rbind(holdout, forecast), method = "gower"))
    lorentzian <- suppressMessages(distance(rbind(holdout, forecast), method = "lorentzian"))
    clark <- suppressMessages(distance(rbind(holdout, forecast), method = "clark"))

    out <- round(c(avg, tanimoto, hassebrook, jaccard, taneja, canberra, gower, lorentzian, clark), 4)
  }

  return(out)
}


###
trunc_t <- function(n, a, b, df, ncp)
{
  min <- suppressWarnings(pt(a, df, ncp))
  max <- suppressWarnings(pt(b, df, ncp))
  u <- suppressWarnings(runif(n, min, max))
  out <- suppressWarnings(qt(u, df, ncp))
  return(out)
}

###
trunc_norm <- function(n, a, b, mean, sd)
{
  min <- suppressWarnings(pnorm(a, mean, sd))
  max <- suppressWarnings(pnorm(b, mean, sd))
  u <- suppressWarnings(runif(n, min, max))
  out <- suppressWarnings(qnorm(u, mean, sd))
  return(out)
}


###
trunc_cauchy <- function(n, a, b, location, scale)
{
  min <- suppressWarnings(pcauchy(a, location, scale))
  max <- suppressWarnings(pcauchy(b, location, scale))
  u <- suppressWarnings(runif(n, min, max))
  out <- suppressWarnings(qcauchy(u, location, scale))
  return(out)
}


###
trunc_logis <- function(n, a, b, location, scale)
{
  min <- suppressWarnings(plogis(a, location, scale))
  max <- suppressWarnings(plogis(b, location, scale))
  u <- suppressWarnings(runif(n, min, max))
  out <- suppressWarnings(qlogis(u, location, scale))
  return(out)
}


###
best_deriv <- function(ts, max_diff = 3, thresh = 0.001)
{
  pvalues <- vector(mode = "double", length = as.integer(max_diff))

  for(d in 1:(max_diff + 1))
  {
    model <- lm(ts ~ t, data.frame(ts, t = 1:length(ts)))
    pvalues[d] <- with(summary(model), pf(fstatistic[1], fstatistic[2], fstatistic[3],lower.tail=FALSE))
    ts <- diff(ts)
  }

  best <- tail(cumsum(pvalues < thresh), 1)

  return(best)
}

###
quantile_prediction <- function(raw_pred = NULL, seed_pred = NULL, raw_error = NULL, holdout = NULL, tail_value = NULL, ts, ci, doxa,
                                error_scale = "naive", error_benchmark = "naive", n_class = NULL, level_names = NULL, seed = 42)
{
  set.seed(seed)

  not_null<- c(!is.null(seed_pred), !is.null(raw_error))
  if(all(not_null))
  {
    pred_integration <- function(pred_seed, error_raw){as.data.frame(map2(pred_seed, error_raw, ~ .x + sample(.y, size = 1000, replace = TRUE)))}
    raw_pred <- pred_integration(seed_pred, raw_error)
    if(!is.null(tail_value)){raw_pred <- t(apply(raw_pred, 1, function(x) invdiff(x, heads = tail_value, add = FALSE)))}
  }

  if(!is.null(raw_pred))
  {
    if(doxa == TRUE){raw_pred <- doxa_filter(ts, raw_pred, n_class)}
    quants <- sort(unique(c((1-ci)/2, 0.25, 0.5, 0.75, ci+(1-ci)/2)))

    if(is.null(n_class))
    {
      p_stats <- function(x){c(min = suppressWarnings(min(x, na.rm = TRUE)), quantile(x, probs = quants, na.rm = TRUE), max = suppressWarnings(max(x, na.rm = TRUE)), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), mode = suppressWarnings(mlv1(x[is.finite(x)], method = "shorth")), kurtosis = suppressWarnings(kurtosis(x[is.finite(x)], na.rm = TRUE)), skewness = suppressWarnings(skewness(x[is.finite(x)], na.rm = TRUE)))}
      quant_pred <- as.data.frame(t(as.data.frame(apply(raw_pred, 2, p_stats))))
      p_value <- apply(raw_pred, 2, function(x) ecdf(x)(seq(min(raw_pred), max(raw_pred), length.out = 1000)))
      divergence <- c(NA, apply(p_value[,-1, drop = FALSE] - p_value[,-ncol(p_value), drop = FALSE], 2, function(x) abs(max(x, na.rm = TRUE))))
      upside_prob <- c(NA, colMeans(apply(raw_pred[,-1, drop = FALSE]/raw_pred[,-ncol(raw_pred), drop = FALSE], 2, function(x) x > 1)))
      iqr_to_range <- (quant_pred[, "75%"] - quant_pred[, "25%"])/(quant_pred[, "max"] - quant_pred[, "min"])
      median_range_ratio <- (quant_pred[, "max"] - quant_pred[, "50%"])/(quant_pred[, "50%"] - quant_pred[, "min"])
      quant_pred <- as.data.frame(round(cbind(quant_pred, iqr_to_range, median_range_ratio, upside_prob, divergence), 4))
    }

    if(is.numeric(n_class) & !is.null(level_names))
    {
      freq <- function(x) {sapply(0:(n_class-1), function(n) sum(x == n))}
      p_stats <- function(x){c(min = suppressWarnings(min(x, na.rm = TRUE)), quantile(x, probs = quants, na.rm = TRUE), max = suppressWarnings(max(x, na.rm = TRUE)), props = freq(x)/length(x), difformity = sd(freq(x)/length(x))/sd(c(1, rep(0, n_class - 1))), entropy = entropy(freq(x)))}
      quant_pred <- as.data.frame(t(as.data.frame(apply(raw_pred, 2, p_stats))))
      p_value <- apply(raw_pred, 2, function(x) ecdf(x)(0:(n_class - 1)))
      divergence <- c(NA, apply(p_value[,-1, drop = FALSE] - p_value[,-ncol(p_value), drop = FALSE], 2, function(x) abs(max(x, na.rm = TRUE))))
      upgrade_prob <- c(NA, colMeans(apply((raw_pred[,-1, drop = FALSE] + 1)/(raw_pred[,-ncol(raw_pred), drop = FALSE] + 1), 2, function(x) x > 1)))###ADDING ONE TO SOLVE ISSUE WITH CLASS ZERO
      quant_pred <- as.data.frame(round(cbind(quant_pred, upgrade_prob = upgrade_prob, divergence = divergence), 4))

      ###FIXING FACTOR LABELS
      n_quants <- length(quants) + 2###QUANTILES + MIN & MAX
      m_index <- as.matrix(quant_pred[, 1:n_quants] + 1)###INDEXING FROM O TO 1
      releveled <- as.data.frame(matrix(level_names[m_index], nrow(m_index), ncol(m_index)))
      colnames(releveled) <- c("min", paste0(quants * 100, "%"), "max")
      quant_pred <- as.data.frame(cbind(releveled, quant_pred[, c(paste0("props", 1:n_class), "difformity", "entropy", "upgrade_prob", "divergence")]))
      colnames(quant_pred) <- c(colnames(releveled), paste0("prop_", level_names), "difformity", "entropy", "upgrade_prob", "divergence")
    }
  }

  testing_error <- NULL
  if(!is.null(holdout))
  {
    if(is.null(n_class)){
      mean_pred <- colMeans(raw_pred)
      testing_error <- my_metrics(holdout, mean_pred, ts, error_scale, error_benchmark, n_class)}

    if(is.numeric(n_class)){
      testing_error <- apply(apply(raw_pred, 1, function(x) my_metrics(holdout, x, n_class = n_class)), 1, function(x) mean(x, na.rm = TRUE))
    }
  }

  outcome <- list(raw_pred = raw_pred, quant_pred = quant_pred, testing_error = testing_error)
  return(outcome)
}

###
doxa_filter <- function(ts, mat, n_class = NULL)
{
  discrete_check <- all(ts%%1 == 0)
  all_positive_check <- all(ts >= 0)
  all_negative_check <- all(ts <= 0)
  monotonic_increase_check <- all(diff(ts) >= 0)
  monotonic_decrease_check <- all(diff(ts) <= 0)

  monotonic_fixer <- function(x, mode)
  {
    model <- recursive_diff(x, 1)
    vect <- model$vector
    if(mode == 0){vect[vect < 0] <- 0; vect <- invdiff(vect, model$head_value, add = TRUE)}
    if(mode == 1){vect[vect > 0] <- 0; vect <- invdiff(vect, model$head_value, add = TRUE)}
    return(vect)
  }

  if(is.null(n_class))
  {
    if(all_positive_check){mat[mat < 0] <- 0}
    if(all_negative_check){mat[mat > 0] <- 0}
    if(discrete_check){mat <- round(mat)}
    #if(!is.null(bounded)){mat[mat > max(bounded)] <- max(bounded); mat[mat < min(bounded)] <- min(bounded)}
    if(monotonic_increase_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 0)))}
    if(monotonic_decrease_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 1)))}
  }

  if(is.numeric(n_class)){mat <- round(mat); mat[mat > (n_class - 1)] <- (n_class - 1); mat[mat < 1] <- 0}

  return(mat)
}


###
prediction_score <- function(integrated_preds, ground_truth)
{
  pfuns <- apply(integrated_preds, 2, function(x) ecdf(na.omit(x)))
  pvalues <- map2_dbl(pfuns, ground_truth, ~ .x(.y))
  scores <- 1 - 2 * abs(pvalues - 0.5)
  return(scores)
}

###
ranker <- function(df, focus, inverse = NULL, absolute = NULL, reverse = FALSE)
{
  rank_set <- df[, focus, drop = FALSE]
  if(!is.null(inverse)){rank_set[, inverse] <- - rank_set[, inverse]}###INVERSION BY COL NAMES
  if(!is.null(absolute)){rank_set[, absolute] <- abs(rank_set[, absolute])}###ABS BY COL NAMES
  index <- apply(scale(rank_set), 1, mean, na.rm = TRUE)
  if(reverse == FALSE){df <- df[order(index),]}
  if(reverse == TRUE){df <- df[order(-index),]}
  return(df)
}

###
ts_graph <- function(x_hist, y_hist, x_forcat, y_forcat, lower = NULL, upper = NULL, line_size = 1.3, label_size = 11,
                     forcat_band = "seagreen2", forcat_line = "seagreen4", hist_line = "gray43", label_x = "Horizon", label_y= "Forecasted Var", dbreak = NULL, date_format = "%b-%d-%Y")
{
  if(is.character(y_hist)){y_hist <- as.factor(y_hist)}
  if(is.character(y_forcat)){y_forcat <- factor(y_forcat, levels = levels(y_hist))}
  if(is.character(lower)){lower <- factor(lower, levels = levels(y_hist))}
  if(is.character(upper)){upper <- factor(upper, levels = levels(y_hist))}

  n_class <- NULL
  if(is.factor(y_hist)){class_levels <- levels(y_hist); n_class <- length(class_levels)}

  all_data <- data.frame(x_all = c(x_hist, x_forcat), y_all = as.numeric(c(y_hist, y_forcat)))
  forcat_data <- data.frame(x_forcat = x_forcat, y_forcat = as.numeric(y_forcat))

  if(!is.null(lower) & !is.null(upper)){forcat_data$lower <- as.numeric(lower); forcat_data$upper <- as.numeric(upper)}

  plot <- ggplot()+ geom_line(data = all_data, aes_string(x = "x_all", y = "y_all"), color = hist_line, size = line_size)
  if(!is.null(lower) & !is.null(upper)){plot <- plot + geom_ribbon(data = forcat_data, aes_string(x = "x_forcat", ymin = "lower", ymax = "upper"), alpha = 0.3, fill = forcat_band)}
  plot <- plot + geom_line(data = forcat_data, aes_string(x = "x_forcat", y = "y_forcat"), color = forcat_line, size = line_size)
  if(!is.null(dbreak)){plot <- plot + scale_x_date(name = paste0("\n", label_x), date_breaks = dbreak, date_labels = date_format)}
  if(is.null(dbreak)){plot <- plot + xlab(label_x)}
  if(is.null(n_class)){plot <- plot + scale_y_continuous(name = paste0(label_y, "\n"), labels = number)}
  if(is.numeric(n_class)){plot <- plot + scale_y_continuous(name = paste0(label_y, "\n"), breaks = 1:n_class, labels = class_levels)}
  plot <- plot + ylab(label_y)  + theme_bw()
  plot <- plot + theme(axis.text=element_text(size=label_size), axis.title=element_text(size=label_size + 2))

  return(plot)
}


###
sampler <- function(vect, n_samp, range = NULL, integer = FALSE, fun = NULL, multi = NULL)
{
  #if(is.character(vect)){vect <- match.arg(arg = vect, choices = range, several.ok = TRUE)}

  if(is.null(vect) & is.null(fun))
  {
    if(!is.character(range)){if(integer){set <- min(range):max(range)} else {set <- seq(min(range), max(range), length.out = 1000)}} else {set <- range}
    if(is.null(multi)){samp <- sample(set, n_samp, replace = TRUE)}
    if(is.numeric(multi)){samp <- replicate(n_samp, sample(set, multi, replace = TRUE), simplify = FALSE)}
  }

  if(is.null(vect) & !is.null(fun)){samp <- fun}

  if(is.null(multi)){
    if(length(vect)==1){samp <- rep(vect, n_samp)}
    if(length(vect) > 1){samp <- sample(vect, n_samp, replace = TRUE)}}

  if(is.numeric(multi)){
    if(length(vect)==1){samp <- replicate(n_samp, rep(vect, multi), simplify = FALSE)}
    if(length(vect) > 1){samp <- replicate(n_samp, sample(vect, multi, replace = TRUE), simplify = FALSE)}}

  return(samp)
}

###
expand <- function(mat, vect)
{
  unfiltered <- mapply(function(i) t(replicate(vect[i], mat[i,], simplify = T)), i = 1:length(vect))
  filtered <- keep(unfiltered, ~length(.x)>0)
  expanded <- as.matrix(Reduce(rbind, filtered))
  return(expanded)
}

###
minmax_transf <- function(vect, min_r, max_r)
{
  v_range <- diff(range(vect))
  rescaled <- (vect - min(vect))/v_range * (max_r - min_r) + min_r
  return(rescaled)
}

###
most_frequent <- function(vect)
{
  freqs <- table(vect)
  max_freq <- max(freqs)
  selected <- names(freqs[freqs == max_freq])
  most_freq <- sample(c(selected), 1)###IN CASE OF MULTIMODAL, RANDOM SELECTION AMONG SELECTED VALUES
  if(is.numeric(vect)){most_freq <- as.numeric(most_freq)}
  if(is.factor(vect)){most_freq <- factor(most_freq, levels = levels(vect))}
  return(most_freq)
}

