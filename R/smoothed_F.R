## Smoothed version of empirical distribution -- Adimari (1998)

#' @export
Fs <- function(X, t) {
  dt <- sort(X)
  n <- length(X)
  if (t < dt[1]) res <- 0
  else if (t == dt[1]) res <- 1 / (2 * n)
  else if (t > dt[n]) res <- 1
  else {
    id <- which(dt >= t)[1] # i + 1
    F1 <- (2 * (id - 1) - 1) / (2 * n) # mean(X <= dt[id])
    F2 <- (2 * id - 1) / (2 * n) # mean(X <= dt[id + 1])
    res <- F1 + (F2 - F1) * (t - dt[id - 1]) / (dt[id] - dt[id - 1])
  }
  return(res)
}

## Smoothed version of F in present of ties

#' @export
Fs_ties <- function(X, t) {
  dt <- sort(unique(X))
  n <- length(dt)
  if (t < dt[1]) res <- 0
  else if (t == dt[1]) res <- 1 / (2 * length(X))
  else if (t > dt[n]) res <- 1
  else {
    id <- which(dt >= t)[1] # i + 1
    if (id == 2) {
      F1 <- 0.5 * mean(X <= dt[id - 1])
    } else {
      F1 <- (mean(X <= dt[id - 1]) + mean(X <= dt[id - 2])) * 0.5
    }
    F2 <- (mean(X <= dt[id]) + mean(X <= dt[id - 1])) * 0.5
    res <- F1 + (F2 - F1) * (t - dt[id - 1]) / (dt[id] - dt[id - 1])
  }
  return(res)
}
