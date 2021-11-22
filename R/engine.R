#' support functions for tetragon
#'
#' @author Giancarlo Vercellino \email{giancarlo.vercellino@gmail.com}
#'
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
#' @importFrom scales number
#'
#' @param df A data frame with time features on columns
#' @param seq_len Positive integer. Time-step number of the projected sequence
#' @param deriv Integer vector. Number of differentiation operations to perform for each original time feature. 0 = no change; 1: one diff; 2: two diff.
#' @param ci Confidence interval. Default: 0.8
#' @param method Positive integer. Distance method for calculating neighbors. Possibile options are: "euclidean", "manhattan", "chebyshev", "sorensen", "gower", "soergel", "kulczynski_d", "canberra", "lorentzian", "intersection", "non-intersection", "wavehedges", "czekanowski", "motyka", "kulczynski_s", "tanimoto", "ruzicka", "inner_product", "harmonic_mean", "cosine", "hassebrook", "jaccard", "dice",  "fidelity",  "bhattacharyya", "squared_chord", "squared_euclidean", "pearson", "neyman", "squared_chi", "prob_symm", "divergence", "clark", "additive_symm", "taneja", "kumar-johnson", "avg".
#' @param distr String. DIstribution used to calculate kernel densities. Possible options are: "norm", "cauchy", "logis", "t", "exp".
#' @param measure_error Logical. TRUE for measuring validation error. FALSE otherwise.
#' @param n_sim String. Sequencing method: deterministic ("segmented"), or non-deterministic ("sampled").
#' @param dates Date. Vector with dates for time features.


###
engine <- function(df, seq_len, deriv, ci, method, distr, measure_error, n_sim, dates)
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

  recursive_diff <- function(vector, deriv)
  {
    head_value <- vector("numeric", deriv)
    tail_value <- vector("numeric", deriv)
    if(deriv==0){head_value = NULL; tail_value = NULL}
    if(deriv > 0){for(i in 1:deriv){head_value[i] <- head(vector, 1); tail_value[i] <- tail(vector, 1); vector <- diff(vector)}}
    outcome <- list(vector = vector, head_value = head_value, tail_value = tail_value)
    return(outcome)
  }

  invdiff <- function(vector, heads)
  {
    if(is.null(heads)){return(vector)}
    for(d in length(heads):1){vector <- cumsum(c(heads[d], vector))}
    return(vector)
  }

  ###
  sequential_kld <- function(m)
  {
    matrix <- t(as.matrix(m))
    n <- nrow(matrix)
    if(n == 1){return(NA)}
    dens <- apply(matrix, 1, function(x) tryCatch(density(x[is.finite(x)], from = min(matrix[is.finite(matrix)]), to = max(matrix[is.finite(matrix)])), error = function(e) NA))
    backward <- dens[-n]
    forward <- dens[-1]

    finite_index <- map2(forward, backward, ~ is.finite(log(.x$y/.y$y)) & is.finite(.x$y))
    seq_kld <- pmap_dbl(list(forward, backward, finite_index), ~ sum(..1$y[..3] * log(..1$y/..2$y)[..3]))
    avg_seq_kld <- round(mean(seq_kld), 3)

    ratios <- log(dens[[n]]$y/dens[[1]]$y)
    finite_index <- is.finite(ratios)

    end_to_end_kld <- dens[[n]]$y * log(dens[[n]]$y/dens[[1]]$y)
    end_to_end_kld <- tryCatch(round(sum(end_to_end_kld[finite_index]), 3), error = function(e) NA)
    kld_stats <- rbind(avg_seq_kld, end_to_end_kld)

    return(kld_stats)
  }

  ###
  upside_probability <- function(m)
  {
    matrix <- t(as.matrix(m))
    n <- nrow(matrix)
    if(n == 1){return(NA)}
    growths <- matrix[-1,]/matrix[-n,] - 1
    dens <- apply(growths, 1, function(x) tryCatch(density(x[is.finite(x)], from = min(x[is.finite(x)]), to = max(x[is.finite(x)])), error = function(e) NA))
    not_na <- !is.na(dens)
    avg_upp <- tryCatch(round(mean(map_dbl(dens[not_na], ~ sum(.x$y[.x$x>0])/sum(.x$y))), 3), error = function(e) NA)
    end_growth <- matrix[n,]/matrix[1,] - 1
    end_to_end_dens <- tryCatch(density(end_growth[is.finite(end_growth)], from = min(end_growth[is.finite(end_growth)]), to = max(end_growth[is.finite(end_growth)])), error = function(e) NA)
    if(class(end_to_end_dens) == "density"){last_to_first_upp <- round(sum(end_to_end_dens$y[end_to_end_dens$x>0])/sum(end_to_end_dens$y), 3)} else {last_to_first_upp <- NA}
    upp_stats <- rbind(avg_upp, last_to_first_upp)
    return(upp_stats)
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

  actuals <- NULL
  if(measure_error == TRUE){actuals <- tail(df, seq_len)}
  n_feat <- ncol(df)
  n_ts <- nrow(df)
  feat_names <- colnames(df)

  all_positive_check <- map_lgl(df, ~ ifelse(all(.x >= 0), T, F))
  all_negative_check <- map_lgl(df, ~ ifelse(all(.x < 0), T, F))

  diff_models <- map2(df, deriv, ~ recursive_diff(.x, .y))
  diff_vectors <- map(diff_models, ~ .x$vector)

  segmented <- map(diff_vectors, ~ segmenter(.x, seq_len))
  minrow <- min(map_dbl(segmented, ~ nrow(.x)))
  segmented <- map(segmented, ~ tail(as.data.frame(.x), minrow))
  dist <- map(segmented, ~ suppressMessages(distance(.x, method = method)))
  dist <- map(dist, ~ (.x - min(.x))/(max(.x) - min(.x)))
  #dist <- abind::abind(dist, along = 3)
  #dist <- apply(dist, c(1, 2), function(x) sum(x, na.rm = TRUE))

  reduced <- Reduce(cbind, map(dist, ~ as.vector(.x[upper.tri(.x)])))
  #cor_mat <- cor(reduced)

  if(n_feat > 1){cor_mat <- cor(reduced)}
  if(n_feat == 1){cor_mat <- matrix(1, 1, 1)}

  dist <- abind::abind(dist, along = 3)
  dist <- map(narray::split(cor_mat, along = 1), ~ apply(dist, c(1, 2), weighted.mean, w = .x))

  if(distr == "norm")
  {
    mean_weights <- map(dist, ~ apply(.x, 1, mean, na.rm = TRUE))
    sd_weights <- map(dist, ~ apply(.x, 1, sd, na.rm = TRUE))
    weights <- map2(mean_weights, sd_weights, ~ replicate(n_sim, map2_dbl(.x, .y, ~ trunc_norm(1, a = 0, b = Inf, mean = .x, sd = .y)), simplify = FALSE))
  }

  if(distr == "t")
  {
    df_weights <- map(dist, ~ apply(.x, 1, function(x) length(x)))
    ncp_weights <- map(dist, ~ apply(.x, 1, mean, na.rm = TRUE))
    weights <- map2(df_weights, ncp_weights, ~ replicate(n_sim, map2_dbl(.x, .y, ~ trunc_t(1, a = 0, b = Inf, df = .x, ncp = .y)), simplify = FALSE))
  }

  if(distr == "cauchy")
  {
    location_weights <- map(dist, ~ apply(.x, 1, mean, na.rm = TRUE))
    scale_weights <- map(dist, ~ apply(.x, 1, sd, na.rm = TRUE))
    scale_weights <- map(scale_weights, ~ {.x[.x <= 0] <- mean(.x[.x > 0]); return(.x)})
    weights <- map2(location_weights, scale_weights, ~ replicate(n_sim, map2_dbl(.x, .y, ~ trunc_cauchy(1, a = 0, b = Inf, location = .x, scale = .y)), simplify = FALSE))
  }

  if(distr == "exp")
  {
    rate_weights <- map(dist, ~ 1/apply(.x, 1, sd, na.rm = TRUE))
    weights <- map(rate_weights, ~ replicate(n_sim, map_dbl(.x, ~ suppressWarnings(rexp(1, rate = .x))), simplify = FALSE))
  }

  if(distr == "logis")
  {
    location_weights <- map(dist, ~ apply(.x, 1, mean, na.rm = TRUE))
    scale_weights <- map(dist, ~ apply(.x, 1, sd, na.rm = TRUE))
    weights <- map2(location_weights, scale_weights, ~ replicate(n_sim, map2_dbl(.x, .y, ~ trunc_logis(1, a = 0, b = Inf, location = .x, scale = .y)), simplify = FALSE))
  }

  weights <- map_depth(weights, 2, ~ {.x[.x == Inf] <- suppressWarnings(max(.x[is.finite(.x)], na.rm = TRUE)); return(.x)})
  weights <- map_depth(weights, 2, ~ {.x[.x == -Inf] <- suppressWarnings(min(.x[is.finite(.x)], na.rm = TRUE)); return(.x)})
  weights <- map_depth(weights, 2, ~ {.x[is.nan(.x)] <- suppressWarnings(mean(.x[is.finite(.x)], na.rm = TRUE)); return(.x)})
  inv_weights <- map_depth(weights, 2, ~ 1/.x)

  sample_pred <- map2(segmented, inv_weights, ~ mapply(function(wts) apply(.x * .y[[wts]], 2, sum, na.rm = TRUE)/sum(.y[[wts]], na.rm = TRUE), wts = 1:length(.y)))
  heads <- map(diff_models, ~ .x$tail)
  sample_pred <- map2(sample_pred, heads, ~ t(apply(.x, 2, function(x) tail(invdiff(x, .y), seq_len))))

  sample_pred <- map_if(sample_pred, all_positive_check, ~ {.x[.x < 0] <- 0; return(.x)})
  sample_pred <- map_if(sample_pred, all_negative_check, ~ {.x[.x > 0] <- 0; return(.x)})

  errors <- NULL
  if(measure_error == TRUE)
  {
    errors <- pmap(list(actuals, sample_pred, tail(df, - seq_len)), ~ colMeans(t(apply(..2, 1, function(x) greybox::measures(..1, x, ..3, benchmark = "mean")[-c(14, 15)])), na.rm = TRUE))
    if(n_feat > 1){errors <- as.matrix(Reduce(rbind, errors))}
    if(n_feat == 1){errors <- matrix(unlist(errors), 1); colnames(errors) <- c("ME", "MAE", "MSE", "MPE", "MAPE", "sCE", "sMAE", "sMSE", "MASE", "RMSSE", "rMAE", "rRMSE", "rAME")}
    rownames(errors) <- feat_names
  }

  plot <- NULL
  pred_stats <- NULL
  predictions <- NULL
  if(measure_error == FALSE)
  {
    quants <- sort(unique(c((1-ci)/2, 0.25, 0.5, 0.75, ci+(1-ci)/2)))
    p_stats <- function(x){stats <- round(c(min = suppressWarnings(min(x, na.rm = T)), quantile(x, probs = quants, na.rm = T), max = suppressWarnings(max(x, na.rm = T)), mean = mean(x, na.rm = T), sd = sd(x, na.rm = T), mode = tryCatch(suppressWarnings(mlv1(x[is.finite(x)], method = "shorth")), error = function(e) NA), kurtosis = tryCatch(suppressWarnings(kurtosis(x[is.finite(x)], na.rm = T)), error = function(e) NA), skewness = tryCatch(suppressWarnings(skewness(x[is.finite(x)], na.rm = T)), error = function(e) NA)), 3); return(stats)}
    predictions <- map(sample_pred, ~ t(apply(.x, 2, function(x) p_stats(x))))
    names(predictions) <- feat_names
    if(is.Date(dates)){new_dates<- tail(dates, seq_len) + seq_len} else {new_dates <- 1:seq_len}
    predictions <- map(predictions, ~ {.x <- as.data.frame(.x); rownames(.x) <- new_dates; return(.x)})

    avg_iqr_to_range <- round(map_dbl(predictions, ~ mean((.x[,"75%"] - .x[,"25%"])/(.x[,"max"] - .x[,"min"]))), 3)
    last_to_first_iqr <- round(map_dbl(predictions, ~ (.x[seq_len,"75%"] - .x[seq_len,"25%"])/(.x[1,"75%"] - .x[1,"25%"])), 3)
    iqr_stats <- as.data.frame(rbind(avg_iqr_to_range, last_to_first_iqr))
    rownames(iqr_stats) <- c("avg_iqr_to_range", "terminal_iqr_ratio")
    colnames(iqr_stats) <- feat_names

    avg_risk_ratio <- round(map_dbl(predictions, ~ mean((.x[,"max"] - .x[,"50%"])/(.x[,"50%"] - .x[,"min"]))), 3)
    last_to_risk_ratio <- round(map_dbl(predictions, ~ ((.x[seq_len,"max"] - .x[seq_len,"50%"])/(.x[seq_len,"50%"] - .x[seq_len,"min"]))/((.x[1,"max"] - .x[1,"50%"])/(.x[1,"50%"] - .x[1,"min"]))), 3)
    risk_stats <- as.data.frame(rbind(avg_risk_ratio, last_to_risk_ratio))
    rownames(risk_stats) <- c("avg_risk_ratio", "terminal_risk_ratio")
    colnames(risk_stats) <- feat_names

    kld_stats <- map(sample_pred, ~ tryCatch(sequential_kld(.x), error = function(e) c(NA, NA)))
    kld_stats <- as.data.frame(map(kld_stats, ~.x))
    rownames(kld_stats) <- c("avg_kl_divergence", "terminal_kl_divergence")
    colnames(kld_stats) <- feat_names

    upp_stats <- map(sample_pred, ~ tryCatch(upside_probability(.x), error = function(e) c(NA, NA)))
    upp_stats <- as.data.frame(map(upp_stats, ~.x))
    rownames(upp_stats) <- c("avg_upside_prob", "terminal_upside_prob")
    colnames(upp_stats) <- feat_names

    pred_stats <- rbind(iqr_stats, risk_stats, kld_stats, upp_stats)

    if(is.Date(dates)){x_hist <- dates; x_forcat <- new_dates} else {x_hist <- 1:n_ts; x_forcat <- (n_ts + 1):(n_ts + seq_len)}

    plot <- pmap(list(predictions, df, feat_names), ~ ts_graph(x_hist = x_hist, y_hist = ..2, x_forcat = x_forcat, y_forcat = ..1[, 3], lower = ..1[, 1], upper = ..1[, 5], label_y = ..3))
  }

  outcome <- list(predictions = predictions, errors = errors, pred_stats = pred_stats, plot = plot)
  return(outcome)
}
