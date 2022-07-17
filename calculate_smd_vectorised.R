# File  : calculate_smd_vectorised.R
# Date  : 19/10/2021
# Author: Alexander Holloway
library(dplyr)

calculate_s2 <- function(x) {
  x <- data.frame(x)
  n <- x$sample_size
  sd <- x$sd
  t <- nrow(x) - 1
  s2 <- sum((n - 1) * sd^2) / (sum(n) - t - 1)
  x$s2 <- s2
  return(x)
}

calculate_d <- function(x) {
  x <- data.frame(x)
  s <- sqrt(x$s2)
  y <- x$score
  y0 <- y[1]
  d <- (y0 - y) / s
  x$d <- d
  return(x)
}

calculate_df <- function(x) {
  x <- data.frame(x)
  n <- x$sample_size
  nn <- mean(n)
  s <- sqrt(x$s2)
  sd <- x$sd
  t <- nrow(x) - 1
  df <- ((nn - 1) * s^4) * (t + 1)^2 / sum(0.25 * (sd %*% t(sd)))
  x$df <- df
  return(x)
}

calculate_var <- function(x) {
  study.design <- x$study_design[1]
  t <- nrow(x)
  n <- x$sample_size
  nn <- mean(n)
  if (study.design == 1) {
    s2 <- x$s2
    df <- x$df
    sd <- x$sd
    sd0 <- sd[1]
    var <- ((sd0^2 + sd^2) - (sd0 * sd)) / (n * s2)
  } else if (study.design == 2) {
    n0 <- n[1]
    d <- x$d
    var <- (1 / n0) + (1 / n) + (d^2 / (2 * (sum(n) - t - 1)))
  }
  x$var <- var
  return(x)
}


calculate_cov <- function(x,maxnt) {
  study.design <- x$study_design[1]
  t <- nrow(x)
  n <- x$sample_size
  d  <- x$d
  if (study.design == 1) {
    nn <- mean(n)
    sd <- x$sd
    df <- x$df
    s2 <- x$s2
    var <- x$var
    mat <- matrix(sd, nrow=t, ncol=t)
    sd0 <- sd[1]
    cov <- (sd0^2 + (0.5 * sd %*% t(sd)) - (0.5 * sd0 * (mat + t(mat)))) / (nn * s2) + ((d %*% t(d)) / (2 * df))
  } else if (study.design == 2) {
    n0 <- n[1]
    cov <- (1 / n0) + ((d %*% t(d)) / (2 * (sum(n) - t - 1)))
  }
  cov.v <- NULL
  for (i in 2:t) {
    for (j in i:t) {
      if (i != j) {
        cov.v <- c(cov.v, cov[i, j])
      }
    }
  }
  maxcov=dim(combn(maxnt-1,2))[2]
  covnames=paste("cov",apply(combn(2:maxnt,2),2,function(x)paste(x[1],x[2],sep="")),sep="")
  cov.v <- c(cov.v, rep(NA, maxcov - length(cov.v)))
  cov <- do.call("rbind", replicate(t, data.frame(cov.v)))
  colnames(cov) <-covnames
  x <- cbind(x, cov)
  return(x)
}

calculate_smd <- function(x, record_id="record_id", population="population", condition="condition", scale="scale", timepoint="timepoint", score="score", sd="sd", sample_size="sample_size", study_design="study_design") {
  x <- x[, c(record_id, population, condition, scale, timepoint, score, sd, sample_size, study_design)]
  maxnt=max(x %>%group_by(record_id, scale, population) %>% add_count() %>% ungroup %>%  select(n))
  
  colnames(x) <- c("record_id", "population", "condition", "scale", "timepoint", "score", "sd", "sample_size", "study_design")
  out <- x %>%
    group_by(record_id, scale, population) %>%
    filter(n() > 1) %>%
    group_modify(~ calculate_s2(.x)) %>%
    group_modify(~ calculate_d(.x)) %>%
    group_modify(~ calculate_df(.x)) %>%
    group_modify(~ calculate_var(.x)) %>%
    group_modify(~ calculate_cov(.x,maxnt))
  return(out)
}

# Example usage
# # longi <- read.csv("in/longicont.csv")
# test <- longicont %>%
# filter(!is.na(sd)) %>%
# filter(study_design != 0) %>%
# group_by(record_id, scale, population) %>%
# select(record_id, population, condition, scale, timepoint, score, sd, sample_size, study_design) %>%
# filter(n() > 1)
