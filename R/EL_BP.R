### Empirical likelihood method of Dong and Tian (2014)

#' @import stats
#' @import utils
#' @import KernSmooth
#' @import BB

Phi <- function(x, y, z){
  Phi <- rep(0, length(y))
  Phi[x < y & y < z] <- 1
  Phi[x == y & y < z | x < y & y == z] <- 1/2
  Phi[x == y & y == z] <- 1/6
  return(Phi)
}

lamb_func <- function(par, W){
  sum(W / (1 + par * W))
}

ll_P2 <- function(Y, tau, theta){
  Wj <- Phi(tau[1], Y, tau[2]) - theta
  lamb_est <- suppressWarnings(BBsolve(par = 0, fn = lamb_func, W = Wj,
                                       method = 2, quiet = TRUE)$par)
  term <- 1 + lamb_est * Wj
  if (any(term <= 0)) res <- Inf
  else res <- 2 * sum(log(term))
  return(res)
}

var_BP <- function(X1, X2, X3, n1, n2, n3, p1, p2, p3, tau_est){
  f1hat <- bkde(X1, kernel = "normal")
  f2hat <- bkde(X2, kernel = "normal")
  f3hat <- bkde(X3, kernel = "normal")
  #
  index1 <- which.min((abs(f1hat$x - tau_est[1])))
  f1hat_phi <- f1hat$y[index1]
  index2 <- which.min(abs(f2hat$x - tau_est[1]))
  f2hat_phi_1 <- f2hat$y[index2]
  index2 <- which.min(abs(f2hat$x - tau_est[2]))
  f2hat_phi_3 <- f2hat$y[index2]
  index3 <- which.min(abs(f3hat$x - tau_est[2]))
  f3hat_phi <- f3hat$y[index3]
  #
  term1 <- p2 * (1 - p2) / n2
  term2 <- p1 * (1 - p1) * f2hat_phi_1^2 / (n1 * f1hat_phi^2)
  term3 <- p3 * (1 - p3) * f2hat_phi_3^2 / (n3 * f3hat_phi^2)
  return(term1 + term2 + term3)
}

boot_BP <- function(X1, X2, X3, n1, n2, n3, p1, p3, B){
  out_bst <- sapply(1:B, function(i){
    x1_b <- sample(X1, size = n1, replace = TRUE)
    x2_b <- sample(X2, size = n2, replace = TRUE)
    x3_b <- sample(X3, size = n3, replace = TRUE)
    tau_empi <- c(quantile(x1_b, p1), quantile(x3_b, 1 - p3))
    P2_est <- mean(x2_b <= tau_empi[2]) - mean(x2_b <= tau_empi[1])
    return(P2_est)
  })
  P2_bst_var <- var(out_bst)
  return(P2_bst_var)
}

#' @export
ELP <- function(X1, X2, X3, n1, n2, n3, tau, tcf1, tcf2, tcf3) {
  ll_est_P <- ll_P2(Y = X2, tau = tau, theta = tcf2)
  if (is.na(ll_est_P)) ll_est_P <- Inf
  ll_est_P_adj <- ll_est_P
  if (!is.infinite(ll_est_P)) {
    tcf2_emp <- mean(X2 <= tau[2]) - mean(X2 <= tau[1])
    var_tcf2_emp <- tcf2_emp * (1 - tcf2_emp)
    var_BP_KDE <- var_BP(X1 = X1, X2 = X2, X3 = X3, n1 = n1, n2 = n2, n3 = n3,
                         p1 = tcf1, p2 = tcf2_emp, p3 = tcf3, tau_est = tau)
    r_est_P <- var_tcf2_emp / (n2 * var_BP_KDE)
    ll_est_P_adj <- ll_est_P_adj * r_est_P
  }
  return(ll_est_P_adj)
}

#' @export
ELB <- function(X1, X2, X3, n1, n2, n3, tau, tcf1, tcf2, tcf3, B = 500) {
  ll_est_P <- ll_P2(Y = X2, tau = tau, theta = tcf2)
  if (is.na(ll_est_P)) ll_est_P <- Inf
  ll_est_BP_adj <- ll_est_P
  if (!is.infinite(ll_est_P)){
    tcf2_emp <- mean(X2 <= tau[2]) - mean(X2 <= tau[1])
    var_tcf2_emp <- tcf2_emp * (1 - tcf2_emp)
    var_BP_bst <- boot_BP(X1 = X1, X2 = X2, X3 = X3, n1 = n1, n2 = n2, n3 = n3,
                          p1 = tcf1, p3 = tcf3, B = B)
    r_est_BP <- var_tcf2_emp / (n2 * var_BP_bst)
    ll_est_BP_adj <- ll_est_BP_adj * r_est_BP
  }
  return(ll_est_BP_adj)
}

#' @export
BTII <- function(X1, X2, X3, n1, n2, n3, p1, p2, p3, B = 500, ci.level){
  z_alp <- qnorm((1 + ci.level)/2)^2
  # bootstrap step
  out_bts <- sapply(1:B, function(i){
    x1_b <- sample(X1, size = n1, replace = TRUE)
    x2_b <- sample(X2, size = n2, replace = TRUE)
    x3_b <- sample(X3, size = n3, replace = TRUE)
    tau_empi <- c(quantile(x1_b, p1), quantile(x3_b, 1 - p3))
    n_P2 <- sum(x2_b >= tau_empi[1] & x2_b <= tau_empi[2])
    P2_bts <- (n_P2 + 0.5*z_alp) / (n2 + z_alp)
    return(P2_bts)
  })
  P2_est <- rowMeans(out_bts)
  sd_P2_est <- apply(out_bts, 1, sd)
  z_P2 <- abs(P2_est - p2)/sd_P2_est
  return(z_P2)
}



