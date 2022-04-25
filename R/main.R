#' tetragon
#'
#' @param df A data frame with time features as columns.
#' @param seq_len Positive integer. Time-step number of the projected sequence. Default: NULL (random selection between maximum boundaries).
#' @param ci Confidence interval. Default: 0.8.
#' @param method String. Distance method for calculating distance matrix among sequences. Options are: "euclidean", "manhattan", "chebyshev", "gower", "lorentzian", "jaccard", "dice", "squared_euclidean", "divergence", "clark", "avg". Default: NULL (random selection among all possible options).
#' @param distr String. Distribution used to expand the distance matrix. Options are: "norm", "cauchy", "logis", "t", "exp", "empirical". Default: NULL (random selection among all possible options).
#' @param n_windows Positive integer. Number of validation tests to measure/sample error. Default: 3.
#' @param n_sample Positive integer. Number of samples for grid or random search. Default: 30.
#' @param dates Date. Vector with dates for time features.
#' @param error_scale String. Scale for the scaled error metrics. Two options: "naive" (average of naive one-step absolute error for the historical series) or "deviation" (standard error of the historical series). Default: "naive".
#' @param error_benchmark String. Benchmark for the relative error metrics. Two options: "naive" (sequential extension of last value) or "average" (mean value of true sequence). Default: "naive".
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
#' \item predictions: min, max, q25, q50, q75, quantiles at selected ci, mean, sd, mode, skewness, kurtosis, IQR to range, risk ratio, upside probability, divergence and entropy for each point fo predicted sequences
#' \item testing_errors: testing errors for one-step and sequence for each ts feature (me, mae, mse, rmsse, mpe, mape, rmae, rrmse, rame, mase, smse, sce, gmrae, pred_score)
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
#' @importFrom stats weighted.mean ecdf lm rnorm
#' @importFrom imputeTS na_kalman
#' @importFrom narray split
#' @importFrom utils head tail
#' @importFrom modeest mlv1
#' @importFrom moments skewness kurtosis
#' @importFrom stats dcauchy dexp dlogis dnorm dt pcauchy pexp plogis pnorm pt median sd density qcauchy qlogis qnorm qt quantile runif cor
#' @importFrom purrr map map2 pmap_dbl map_dbl map_lgl map2_dbl
#' @import abind
#' @import philentropy
#' @import abind
#' @import ggplot2
#' @import greybox
#' @import stringr
#' @importFrom dqrng dqsample
#' @importFrom scales number
#' @importFrom entropy entropy


#'@examples
#'\donttest{
#'tetragon(covid_in_europe[, c(2, 4)], seq_len = 40, n_sample = 2)
#'}

tetragon <- function(df, seq_len = NULL, ci = 0.8, method = NULL, distr = NULL, n_windows = 3, n_sample = 30, dates = NULL, error_scale = "naive", error_benchmark = "naive", seed = 42)
{

  windower <- function(ts, seq_len, deriv, ci, method, distr, n_windows, dates, error_scale, error_benchmark, feat_name)
  {
    ####SUPPORT ENGINE
    engine <- function(ts, seq_len, deriv, ci, method, distr, measure_error, n_sim, dates, error_scale, error_benchmark, collected_errors, feat_name)
    {
      segmenter <- function(vector, seq_len)
      {
        if(seq_len > length(vector)){stop("vector too short for sequence length")}
        n_segment <- floor(length(vector)/seq_len)
        residual <- length(vector)%%seq_len
        fixed_vector <- vector
        if(residual > 0){fixed_vector <- tail(vector, - residual)}
        segment_id <- sort(rep(1:n_segment, seq_len))
        seq_list <- narray::split(fixed_vector, along = 1, subsets = segment_id)
        segmented <- as.data.frame(Reduce(rbind, seq_list))
        rownames(segmented) <- NULL
        return(segmented)
      }

      ###
      prediction_score <- function(integrated_preds, ground_truth)
      {
        pfuns <- apply(integrated_preds, 2, ecdf)
        pvalues <- purrr::map2_dbl(pfuns, ground_truth, ~ .x(.y))
        scores <- mean(1 - 2 * abs(pvalues - 0.5))
        return(scores)
      }

      ####
      doxa_filter <- function(orig, mat, n_class = NULL)
      {
        discrete_check <- all(orig%%1 == 0)
        all_positive_check <- all(orig >= 0)
        all_negative_check <- all(orig <= 0)
        monotonic_increase_check <- all(diff(orig) >= 0)
        monotonic_decrease_check <- all(diff(orig) <= 0)
        class_check <- FALSE
        if(is.integer(n_class)){class_check <- length(unique(orig)) <= n_class}

        monotonic_fixer <- function(x, mode)
        {
          model <- recursive_diff(x, 1)
          vect <- model$vector
          if(mode == 0){vect[vect < 0] <- 0; vect <- invdiff(vect, model$head_value, add = T)}
          if(mode == 1){vect[vect > 0] <- 0; vect <- invdiff(vect, model$head_value, add = T)}
          return(vect)
        }

        if(all_positive_check){mat[mat < 0] <- 0}
        if(all_negative_check){mat[mat > 0] <- 0}
        if(discrete_check){mat <- floor(mat)}
        if(monotonic_increase_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 0)))}###SOMETHING SMARTER???
        if(monotonic_decrease_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 1)))}
        if(class_check){mat[!(mat %in% unique(orig))] <- ((mat[!(mat %in% unique(orig))] > max(unique(orig))) * max(unique(orig))) + ((mat[!(mat %in% unique(orig))] < min(unique(orig))) * min(unique(orig)))}

        return(mat)
      }

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
      my_metrics <- function(holdout, forecast, actual, error_scale, error_benchmark)
      {
        scale <- switch(error_scale, "deviation" = sd(actual), "naive" = mean(abs(diff(actual))))
        benchmark <- switch(error_benchmark, "average" = rep(mean(forecast), length(forecast)), "naive" = rep(tail(actual, 1), length(forecast)))

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
        return(out)
      }


      ts_graph <- function(x_hist, y_hist, x_forcat, y_forcat, lower = NULL, upper = NULL, line_size = 1.3, label_size = 11,
                           forcat_band = "darkorange", forcat_line = "darkorange", hist_line = "gray43", label_x = "Horizon", label_y= "Forecasted Var", dbreak = NULL, date_format = "%b-%d-%Y")
      {
        all_data <- data.frame(x_all = c(x_hist, x_forcat), y_all = c(y_hist, y_forcat))
        forcat_data <- data.frame(x_forcat = x_forcat, y_forcat = y_forcat)

        if(!is.null(lower) & !is.null(upper)){forcat_data$lower <- lower; forcat_data$upper <- upper}

        plot <- ggplot()+geom_line(data = all_data, aes_string(x = "x_all", y = "y_all"), color = hist_line, size = line_size)
        if(!is.null(lower) & !is.null(upper)){plot <- plot + geom_ribbon(data = forcat_data, aes(x = x_forcat, ymin = lower, ymax = upper), alpha = 0.3, fill = forcat_band)}
        plot <- plot + geom_line(data = forcat_data, aes(x = x_forcat, y = y_forcat), color = forcat_line, size = line_size)
        if(!is.null(dbreak)){plot <- plot + scale_x_date(name = paste0("\n", label_x), date_breaks = dbreak, date_labels = date_format)}
        if(is.null(dbreak)){plot <- plot + xlab(label_x)}
        plot <- plot + scale_y_continuous(name = paste0(label_y, "\n"), labels = number)
        plot <- plot + ylab(label_y)  + theme_bw()
        plot <- plot + theme(axis.text=element_text(size=label_size), axis.title=element_text(size=label_size + 2))

        return(plot)
      }

      trunc_t <- function(n, a, b, df, ncp)
      {
        min <- suppressWarnings(pt(a, df, ncp))
        max <- suppressWarnings(pt(b, df, ncp))
        u <- suppressWarnings(runif(n, min, max))
        out <- suppressWarnings(qt(u, df, ncp))
        return(out)
      }

      trunc_norm <- function(n, a, b, mean, sd)
      {
        min <- suppressWarnings(pnorm(a, mean, sd))
        max <- suppressWarnings(pnorm(b, mean, sd))
        u <- suppressWarnings(runif(n, min, max))
        out <- suppressWarnings(qnorm(u, mean, sd))
        return(out)
      }

      trunc_cauchy <- function(n, a, b, location, scale)
      {
        min <- suppressWarnings(pcauchy(a, location, scale))
        max <- suppressWarnings(pcauchy(b, location, scale))
        u <- suppressWarnings(runif(n, min, max))
        out <- suppressWarnings(qcauchy(u, location, scale))
        return(out)
      }

      trunc_logis <- function(n, a, b, location, scale)
      {
        min <- suppressWarnings(plogis(a, location, scale))
        max <- suppressWarnings(plogis(b, location, scale))
        u <- suppressWarnings(runif(n, min, max))
        out <- suppressWarnings(qlogis(u, location, scale))
        return(out)
      }

      actual <- NULL
      if(measure_error == TRUE){actual <- tail(ts, seq_len); df <- head(ts, - seq_len)}
      n_ts <- length(ts)

      diff_model <- recursive_diff(ts, deriv)
      diff_vector <- diff_model$vector

      segmented <- segmenter(diff_vector, seq_len)
      dist <- suppressMessages(distance(segmented, method = method))
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
      weights <- map(weights, ~ {.x[is.nan(.x)] <- suppressWarnings(mean(.x[is.finite(.x)], na.rm = TRUE)); return(.x)})

      inv_weights <- map(weights, ~ 1/.x)

      sample_pred_seeds <- t(mapply(function(wts) apply(segmented * inv_weights[[wts]], 2, sum, na.rm = TRUE)/sum(inv_weights[[wts]], na.rm = TRUE), wts = 1:length(inv_weights)))
      heads <- diff_model$tail_value
      pred_seeds <- t(apply(sample_pred_seeds, 1, function(x) tail(invdiff(x, heads), seq_len)))

      errors <- NULL
      raw_errors <- NULL
      pred_scores <- NULL
      if(measure_error == TRUE)
      {
        raw_errors <- t(mapply(function(i) actual - pred_seeds[i,],  i = 1:nrow(pred_seeds)))
        avg_pred_seeds <- colMeans(pred_seeds)

        if(!is.na(collected_errors))
        {
          collected_errors <- as.data.frame(collected_errors)
          sample_pred <- mapply(function(t) avg_pred_seeds[t] + sample(collected_errors[, t], size = 1000, replace = T), t = 1:seq_len)
          sample_pred <- doxa_filter(ts, sample_pred, n_class = NULL)
          pred_scores <- prediction_score(sample_pred, actual)
        }

        errors <- my_metrics(holdout = actual, forecast = avg_pred_seeds, actual = head(ts, - seq_len), error_scale, error_benchmark)
        errors <- t(as.data.frame(errors))
        rownames(errors) <- feat_name
      }


      plot <- NULL
      pred_stats <- NULL
      predictions <- NULL
      if(measure_error == FALSE)
      {
        collected_errors <- as.data.frame(collected_errors)
        avg_pred_seeds <- colMeans(pred_seeds)
        sample_pred <- mapply(function(t) avg_pred_seeds[t] + sample(collected_errors[, t], size = 1000, replace = T), t = 1: seq_len)
        sample_pred <- doxa_filter(ts, sample_pred, n_class = NULL)

        quants <- sort(unique(c((1-ci)/2, 0.25, 0.5, 0.75, ci+(1-ci)/2)))
        p_stats <- function(x){stats <- round(c(min = suppressWarnings(min(x, na.rm = T)), quantile(x, probs = quants, na.rm = T), max = suppressWarnings(max(x, na.rm = T)), mean = mean(x, na.rm = T), sd = sd(x, na.rm = T), mode = tryCatch(suppressWarnings(mlv1(x[is.finite(x)], method = "shorth")), error = function(e) NA), kurtosis = tryCatch(suppressWarnings(kurtosis(x[is.finite(x)], na.rm = T)), error = function(e) NA), skewness = tryCatch(suppressWarnings(skewness(x[is.finite(x)], na.rm = T)), error = function(e) NA)), 3); return(stats)}
        predictions <- as.data.frame(t(apply(sample_pred, 2, function(x) p_stats(x))))

        iqr_to_range <- (predictions[, "75%"] - predictions[, "25%"])/(predictions[, "max"] - predictions[, "min"])
        risk_ratio <- (predictions[, "max"] - predictions[, "50%"])/(predictions[, "50%"] - predictions[, "min"])
        growths <- mapply(function(m, s) rnorm(1000, m, s), m = predictions[, "mean"], s = predictions[, "sd"])
        upside_prob <- c(NA, colMeans(apply(growths[,-1]/growths[, - seq_len], 2, function(x) x > 1)))
        pvalues <- as.data.frame(map(as.data.frame(sample_pred), ~ ecdf(.x)(.x)))
        divergence <- c(NA, apply(pvalues[,-1] - pvalues[, - seq_len], 2, function(x) abs(max(x, na.rm = T))))
        entropy <- apply(sample_pred, 2, entropy)
        predictions <- cbind(predictions, iqr_to_range = iqr_to_range, risk_ratio = risk_ratio, upside_prob = upside_prob, divergence = divergence, entropy = entropy)
        if(is.Date(dates)){new_dates<- tail(dates, seq_len) + seq_len} else {new_dates <- 1:seq_len}
        predictions <- as.data.frame(predictions)
        rownames(predictions) <- new_dates

        if(is.Date(dates)){x_hist <- dates; x_forcat <- new_dates} else {x_hist <- 1:n_ts; x_forcat <- (n_ts + 1):(n_ts + seq_len)}

        plot <- ts_graph(x_hist = x_hist, y_hist = ts, x_forcat = x_forcat, y_forcat = predictions[, 3], lower = predictions[, 1], upper = predictions[, 5], label_y = feat_name)

      }

      outcome <- list(predictions = predictions, errors = errors, raw_errors = raw_errors, pred_scores = pred_scores, plot = plot)
      return(outcome)
    }


    ####
    n_ts <- length(ts)
    test_index <- unique(round(seq(seq_len * 4, n_ts, length.out = n_windows + 1)))
    n_test <- length(test_index)
    if((n_test - 1) < n_windows){message("not enough data, testing on ", (n_test - 1), " possible windows")}

    results <- vector("list", n_test)
    collected_errors <- vector("list", n_test - 1)
    for(i in 1:n_test)
      {
      results[[i]] <- engine(ts[1:test_index[i]], seq_len, deriv, ci, method, distr, measure_error = ifelse(i < n_test, TRUE, FALSE), n_sim = 100, dates, error_scale, error_benchmark, collected_errors = ifelse(i == 1, NA, collected_errors[i-1]), feat_name)
      if(i < n_test){collected_errors[[i]] <- results[[i]]$raw_errors}
      }

    window_errors <- head(map(results, ~ .x$errors), -1)
    testing_errors <- Reduce("+", window_errors)/length(window_errors)
    window_pred_scores <- head(tail(map(results, ~ .x$pred_scores), -1), -1)
    pred_scores <- Reduce("+", window_pred_scores)/length(window_pred_scores)
    testing_errors <- cbind(pred_scores, testing_errors)
    model <- results[[n_test]]

    outcome <- list(model = model, testing_errors = testing_errors)

    return(outcome)
  }

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

  set.seed(seed)

  tic.clearlog()
  tic("time")

  if(!is.data.frame(df)){stop("time features must be in dataframe format")}
  if(anyNA(df)){df <- as.data.frame(na_kalman(df))}

  n_ts <- nrow(df)
  n_feat <- ncol(df)
  feat_names <- colnames(df)

  deriv <- map_dbl(df, ~ best_deriv(.x))

  max_limit <- floor((n_ts - max(deriv))/(4 + n_windows + 1))
  if(max_limit < 2){stop("not enough data for the validation windows")}

  if(all(c(length(seq_len) == 1, length(method) == n_feat, length(distr) == n_feat, n_sample == 1))){fixed <- TRUE} else {fixed <- FALSE}
  if(!is.null(method)){method <- match.arg(arg = method, choices = c("euclidean", "manhattan", "chebyshev", "gower", "lorentzian", "jaccard", "dice", "squared_euclidean", "divergence", "clark", "avg"), several.ok = TRUE)}
  if(!is.null(distr)){distr <- match.arg(arg = distr, choices = c("norm", "t", "cauchy", "exp", "logis", "empirical"), several.ok = TRUE)}
  if(!is.null(seq_len) && any(seq_len < 2 | seq_len > max_limit)){suppressWarnings(seq_len[seq_len < 2] <- 2); suppressWarnings(seq_len[seq_len > max_limit] <- max_limit); message("setting seq_len within boundaries")}

  if(fixed == FALSE)
  {
    if(is.null(seq_len)){sl_sample_field <- 2:max_limit} else {sl_sample_field <- seq_len}
    if(is.null(method)){m_sample_field <- c("euclidean", "manhattan", "chebyshev", "gower", "lorentzian",
                                            "jaccard", "dice", "squared_euclidean", "divergence", "clark", "avg")} else {m_sample_field <- method} ###

    if(is.null(distr)){dsr_sample_field <- c("norm", "t", "cauchy", "exp", "logis", "empirical")} else {dsr_sample_field <- distr}

    if(length(sl_sample_field) == 1){sl_sample <- rep(sl_sample_field, n_sample)} else {sl_sample <- sample(sl_sample_field, n_sample, replace = TRUE)}
    m_sample <- replicate(n_sample, sample(m_sample_field, n_feat, replace = TRUE), simplify = F)
    dsr_sample <- replicate(n_sample, sample(dsr_sample_field, n_feat, replace = TRUE), simplify = F)
    sl_sample <- map2(sl_sample, replicate(n_sample, deriv, simplify = F), ~ ((.x - max(.y)) == 0) * (.x + 2) + ((.x - max(.y)) == 1) * (.x + 1) + ((.x - max(.y)) > 1) * .x)
    if(length(seq_len)==1 && any(sl_sample != seq_len)){message("fixing seq_len for available data")}

    hyper_params <- list(sl_sample, m_sample, dsr_sample)
  }

  if(fixed == TRUE){hyper_params <- list(sl_sample = list(seq_len), m_sample = list(method), dsr_sample = list(distr))}

  exploration <- lapply(1:n_feat, function(ft) suppressMessages(pmap(hyper_params, ~ windower(df[, ft, drop = T], seq_len = ..1, deriv[ft], ci, method = ..2[ft], distr = ..3[ft], n_windows, dates, error_scale, error_benchmark, feat_names[ft]))))
  exploration <- transpose(exploration)

  if(fixed == FALSE)
  {
    testing_errors <- map(map_depth(exploration, 2, ~ .x$testing_errors), ~ Reduce(rbind, .x))

    aggr_errors <- t(as.data.frame(map(testing_errors, ~ apply(.x, 2, mean, na.omit = T))))
    colnames(aggr_errors) <- paste0("avg_", colnames(aggr_errors))
    history <- data.frame(seq_len = unlist(sl_sample))
    history$method <- m_sample
    history$distr <- dsr_sample
    history <- data.frame(history, aggr_errors)
    rownames(history) <- NULL
    history <- history[order(- unlist(history$avg_pred_scores), abs(unlist(history$avg_me)), unlist(history$avg_mae)),]

    best_index <- as.numeric(rownames(history[1,]))
  }

  if(fixed == TRUE)
  {
    history <- NULL
    best_index <- 1
  }

  testing_errors <- map(map_depth(exploration, 2, ~ .x$testing_errors), ~ Reduce(rbind, .x))
  predictions <- map(exploration[[best_index]], ~ .x$model$predictions)
  names(predictions) <- feat_names
  plots <- map(exploration[[best_index]], ~ .x$model$plot)
  names(plots) <- feat_names

  best <- list(predictions = predictions, testing_errors = testing_errors[[best_index]], plots = plots)

  toc(log = TRUE)
  time_log <- seconds_to_period(round(parse_number(unlist(tic.log())), 0))

  outcome <- list(exploration = exploration, history = history, best = best, time_log = time_log)

  return(outcome)
}
