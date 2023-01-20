## Influence Function-based Empirical Likelihood -- Hai et al. (2022)

#' @import stats
#' @import utils
#' @import KernSmooth
#' @import rootSolve

#' @export
IF_method <- function(X1, X2, X3, n1, n2, n3, tau, tcf1, tcf2, tcf3) {
  f1hat <- bkde(X1, kernel = "normal")
  f2hat <- bkde(X2, kernel = "normal")
  f3hat <- bkde(X3, kernel = "normal")
  #
  index1 <- which.min((abs(f1hat$x - tau[1])))
  f1hat_phi <- f1hat$y[index1]
  index2 <- which.min(abs(f2hat$x - tau[1]))
  f2hat_phi_1 <- f2hat$y[index2]
  index2 <- which.min(abs(f2hat$x - tau[2]))
  f2hat_phi_3 <- f2hat$y[index2]
  index3 <- which.min(abs(f3hat$x - tau[2]))
  f3hat_phi <- f3hat$y[index3]
  #
  z1 <- ifelse(X1 <= tau[1], 1, 0)
  z2 <- ifelse(tau[1] < X2 & X2 <= tau[2], 1, 0)
  z3 <- ifelse(X3 > tau[2], 1, 0)
  #
  n <- n1 + n2 + n3
  zz <- c((n / n1) * (f2hat_phi_1 / f1hat_phi) * (z1 - tcf1),
          (n / n2) * (z2 - tcf2),
          (n / n3) * (f2hat_phi_3 / f3hat_phi) * (z3 - tcf3))
  ff1 <- function(xx, zz) {
    temp1 <- zz
    temp2 <- 1 + zz * xx
    F1 <- mean(temp1 / temp2)
    return(F1)
  }
  ss <- multiroot(f = ff1, start = 0, maxiter = 100, zz = zz)
  temp <- 1 + ss$root * zz
  if (any(temp <= 0)) {
    ll_IF <- Inf
  } else ll_IF <- 2*sum(log(temp))
  return(ll_IF)
}

