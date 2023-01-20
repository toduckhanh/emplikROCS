####========================================================================####
## This file contains functions for estimating the empirical likelihood ratio ##
## for three TCFs at given pair of thresholds                                 ##
####========================================================================####

#' @import stats
#' @import utils

empi_llike_3C <- function(X1, X2, X3, n1, n2, n3, tcf1, tcf2, tcf3, tau,
                          type_F = c("empi", "Adi", "Adi_ties")) {
  ll <- Inf
  r1 <- range(X1)
  r2 <- range(X2)
  r3 <- range(X3)
  ckt11 <- as.numeric(tau[1] >= r1[1] & tau[1] <= r1[2])
  ckt12 <- as.numeric(tau[1] >= r2[1] & tau[1] <= r2[2])
  ckt21 <- as.numeric(tau[2] >= r2[1] & tau[2] <= r2[2])
  ckt22 <- as.numeric(tau[2] >= r3[1] & tau[2] <= r3[2])
  if(ckt11 & (ckt12 | ckt21) & ckt22){
    type_F <- match.arg(type_F)
    F1_tau1 <- switch(type_F,
                      empi = mean(X1 <= tau[1]),
                      Adi = Fs(X1, tau[1]),
                      Adi_ties = Fs_ties(X1, tau[1])
    )
    F2_tau12 <- switch(type_F,
                       empi = mean(X2 <= tau[2]) - mean(X2 <= tau[1]),
                       Adi = Fs(X2, tau[2]) - Fs(X2, tau[1]),
                       Adi_ties = Fs_ties(X2, tau[2]) - Fs_ties(X2, tau[1])
    )
    F3_tau2 <- switch(type_F,
                      empi = mean(X3 <= tau[2]),
                      Adi = Fs(X3, tau[2]),
                      Adi_ties = Fs_ties(X3, tau[2])
    )
    ll1 <- 2 * n1 * (F1_tau1 * log(F1_tau1 / tcf1) +
                       (1 - F1_tau1) * log((1 - F1_tau1) / (1 - tcf1)))
    if (F2_tau12 == 0) {
      ll2 <- Inf
    } else {
      ll2 <- 2 * n2 * (F2_tau12 * log(F2_tau12 / tcf2) +
                         (1 - F2_tau12) * log((1 - F2_tau12) / (1 - tcf2)))
    }
    ll3 <- 2 * n3 * (F3_tau2 * log(F3_tau2 / (1 - tcf3)) +
                         (1 - F3_tau2) * log((1 - F3_tau2) / tcf3))
    ll <- ll1 + ll2 + ll3
   }
  return(ll)
}

## ---- bootstrap procedure to compute w for EL for TCF2 ----
bts_func_3C <- function(X1, X2, X3, n1, n2, n3, tcf1, tcf2, tcf3,
                        enlarged = TRUE, B, type_F) {
  empi_bts <- sapply(1:B, function(i){
    flag <- 0
    while(flag == 0){
      X1.b <- sample(X1, n1, replace = TRUE)
      X2.b <- sample(X2, n2, replace = TRUE)
      X3.b <- sample(X3, n3, replace = TRUE)
      if(enlarged){
        X1.b <- c(X1.b, min(X1), max(X1))
        X2.b <- c(X2.b, min(X2), max(X2))
        X3.b <- c(X3.b, min(X3), max(X3))
      }
      flag <- as.numeric((mean(X1.b) < mean(X2.b)) * (mean(X2.b) < mean(X3.b)))
    }
    tau1_est <- quantile(X1.b, tcf1, names = FALSE)
    tau2_est <- quantile(X3.b, 1 - tcf3, names = FALSE)
    if (tau1_est > tau2_est) {
      temp <- tau2_est
      tau2_est <- tau1_est
      tau1_est <- temp
    }
    res <- empi_llike_3C(X1 = X1.b, X2 = X2.b, X3 = X3.b, n1 = length(X1.b),
                         n2 = length(X2.b), n3 = length(X3.b), tcf1 = tcf1,
                         tcf2 = tcf2, tcf3 = tcf3, tau = c(tau1_est, tau2_est),
                         type_F = type_F)
    return(res)
  })
  empi_bts[is.na(empi_bts)] <- Inf
  r_est <- ((7 / 9)^3) / median(empi_bts)
  return(r_est)
}

## ---- main function ----

#' @export
EL3C <- function(X1, X2, X3, n1, n2, n3, tcf1, tcf2, tcf3, tau, type_F,
                 enlarged = TRUE, B = 200, type_F_B) {
  ll_est <- empi_llike_3C(X1 = X1, X2 = X2, X3 = X3, n1 = n1, n2 = n2, n3 = n3,
                          tcf1 = tcf1, tcf2 = tcf2, tcf3 = tcf3,
                          tau = tau, type_F = type_F)
  if(is.na(ll_est)) ll_est <- Inf
  ll_est_adj <- ll_est
  if(!is.infinite(ll_est)){
    tcf2_est <- Fs(X2, tau[2]) - Fs(X2, tau[1])
    r_est <- bts_func_3C(X1 = X1, X2 = X2, X3 = X3, n1 = n1, n2 = n2, n3 = n3,
                         tcf1 = tcf1, tcf2 = tcf2_est, tcf3 = tcf3,
                         enlarged = enlarged, B = B, type_F = type_F_B)
    ll_est_adj <- r_est * ll_est_adj
  }
  return(ll_est_adj)
}



