## The following codes are provided by Rahman

#' @import stats
#' @import utils
#' @import emplik

myWdataclean2 <- function (z, wt = rep (1, length (z))) {
  niceorder <- order(z) #get order statistic of z
  sortedz <- z[niceorder] #sorted z, can #sort(z)
  sortedw <- wt[niceorder]
  n <- length (sortedz)
  # y checks for jumps in sortedz using offsets of sortedz
  y <- sortedz[-1] != sortedz[-n]
  #ind stores jump indices ( final index will be n)
  ind <- c(which (y|is.na(y)), n)
  # csum is cumulative sum of the weights
  csumw <- cumsum(sortedw)
  # value contains the ( unique ) obs in sortedz
  # weight has the weights of the obs in sortedz
  list(value = sortedz[ind], weight = diff(c(0, csumw[ind])))
}

myfun5 <- function(x, theta , eps) {
  u <- (x - theta ) * sqrt(5) / eps
  INDE <- (u < sqrt(5)) & (u > -sqrt(5))
  u[u >= sqrt(5)] <- 0
  u[u <= -sqrt(5)] <- 1
  y <- 0.5 - (u - (u)^3 / 15) * 3 / (4 * sqrt(5))
  u[INDE] <- y[INDE]
  return(u)
}

meany <- function(y, c1, c2, eps){
  ny <- length(y)
  vy <- myfun5(y, c2, eps) * (1 - myfun5(y, c1, eps))
  return(vy)
}

eltest <- function (ab = vector("numeric" ,2), x.sample, y.sample, z.sample,
                    mu, eps) {
  # mu=c(P1,P3,P2)
  # ab -vector of quantiles, mu-vector of means, i.e. proportions
  all1 <- el.test(myfun5(x.sample, ab[1], eps), mu[1])
  all3 <- el.test(myfun5(z.sample, ab[2], eps), 1 - mu[2])
  all2 <- el.test(meany(y.sample, ab[1], ab[2], eps), mu[3])
  list('-2LLR' = all1$'-2LLR' + all2$'-2LLR' + all3$'-2LLR')
}

neighb_xyz <- function(sp12 = vector("numeric", 2), x, y, z, true, eps) {
  # true = P2
  # find the threshold
  n1 <- length(x)
  n2 <- length(y)
  n3 <- length(z)
  sortx_w <- myWdataclean2(z = x)
  sortedx <- sortx_w$value
  cx <- cumsum(sortx_w$weight / n1)
  idex1 <- ifelse(sp12[1] == 0, 1, max(which(cx <= sp12[1])))
  ## might have trouble in capturing the first element
  sortz_w <- myWdataclean2(z = z)
  sortedz <- sortz_w$value
  cz <- cumsum(sortz_w$weight / n3)
  idez2 <- ifelse(sp12[2] == 1, 1, max(which(cz <= 1 - sp12[2])))
  if (sortedx[idex1] < sortedz[idez2]) {
    mu <- c(sp12, true) # P1, P3 and P2
    abs1 <- data.frame(x1 = idex1, z2 = idez2)
    res <- eltest(c(sortedx[idex1], sortedz[idez2]), x, y, z, mu, eps)
    best <- res$'-2LLR'
    neg2_llr <- as.vector(best)
    ##############
    ite <- TRUE
    while (ite == TRUE) {
      neib1 <- c(idex1 - 1, idez2)
      if ((idex1 - 1 > 0) & (sortedx[idex1 - 1] < sortedz[idez2])) {
        if ((!any(abs1[, 1] == neib1[1])) ||
            (!any(which(neib1[1] == abs1[, 1]) %in%
                  which(neib1[2] == abs1[, 2])))) {
          abs1 <- rbind(abs1, neib1)
          res <- eltest(ab = c(sortedx[idex1 - 1], sortedz[idez2]), x, y, z,
                        mu, eps)
          neg2_llr <- c(neg2_llr, res$'-2LLR')
        }
      }
      neib2 <- c(idex1, idez2 - 1)
      if ((sortedz[idez2 - 1] > sortedx[idex1]) & (idez2 - 1 > 0)) {
        if ((!any(abs1[, 2] == neib2[2])) ||
            (!any(which(neib2[1] == abs1[, 1]) %in%
                  which(neib2[2] == abs1[, 2])))) {
          abs1 <- rbind(abs1, neib2)
          res <- eltest(ab = c(sortedx[idex1], sortedz[idez2 - 1]), x, y, z,
                        mu, eps)
          neg2_llr <- c(neg2_llr, res$'-2LLR')
        }
      }
      neib3 <- c((idex1 + 1), idez2)
      if ((sortedx[idex1 + 1] < sortedz [idez2]) & (idex1 + 1 <= n1)) {
        if ((!any(abs1[, 1] == neib3[1])) ||
            (!any(which(neib3[1] == abs1[, 1]) %in%
                  which(neib3[2] == abs1[, 2])))) {
          abs1 <- rbind(abs1, neib3)
          res <- eltest(ab = c(sortedx[idex1 + 1], sortedz[idez2]), x, y, z,
                        mu, eps)
          neg2_llr <- c(neg2_llr, res$'-2LLR')
        }
      }
      neib4 <- c(idex1, idez2 + 1)
      if ((idez2 + 1 <= n3) & (sortedx[idex1] < sortedz[idez2 + 1])) {
        if ((!any(abs1[, 2] == neib4[2])) ||
            (!any(which(neib4[1] == abs1[, 1]) %in%
                  which(neib4[2] == abs1[, 2])))) {
          abs1 <- rbind(abs1, neib4)
          res <- eltest(ab = c(sortedx[idex1], sortedz[idez2 + 1]), x, y, z,
                        mu, eps)
          neg2_llr <- c(neg2_llr, res$'-2LLR')
        }
      }
      neib5 <- c(idex1 - 1, idez2 - 1)
      if ((idex1 - 1 > 0) & (idez2 - 1 > 0) &
          (sortedz[idez2 - 1] > sortedx[idex1 - 1])) {
        if (!any(which(neib5[1] == abs1[, 1]) %in%
                 which(neib5[2] == abs1[, 2]))) {
          abs1 <- rbind(abs1, neib5)
          res <- eltest(ab = c(sortedx[idex1 - 1], sortedz[idez2 - 1]), x, y, z,
                        mu, eps)
          neg2_llr <- c(neg2_llr, res$'-2LLR')
        }
      }
      neib6 <- c(idex1 - 1, idez2 + 1)
      if ((idex1 - 1 > 0) & (idez2 + 1 <= n3) &
          (sortedx[idex1 - 1] < sortedz[idez2 + 1])) {
        if (!any(which(neib6[1] == abs1[, 1]) %in%
                 which(neib6[2] == abs1[, 2]))) {
          abs1 <- rbind(abs1, neib6)
          res <- eltest(ab = c(sortedx[idex1 - 1], sortedz[idez2 + 1]), x, y, z,
                        mu, eps)
          neg2_llr <- c(neg2_llr, res$'-2LLR')
        }
      }
      neib7 <- c(idex1 + 1, idez2 - 1)
      if ((sortedx[idex1 + 1] < sortedz[idez2 - 1]) &
          (idex1 + 1 <= n1) & (idez2 - 1 > 0)) {
        if (!any(which(neib7[1] == abs1[, 1]) %in%
                 which(neib7[2] == abs1[, 2]))) {
          abs1 <- rbind(abs1, neib7)
          res <- eltest(ab = c(sortedx[idex1 + 1], sortedz[idez2 - 1]), x, y, z,
                        mu, eps)
          neg2_llr <- c(neg2_llr, res$'-2LLR')
        }
      }
      neib8 <- c(idex1 + 1, idez2 + 1)
      if ((sortedx[idex1 + 1] < sortedz[idez2 + 1]) & (idez2 + 1 <= n3) &
          (idex1 + 1 <= n1)) {
        if (!any(which(neib8[1] == abs1[, 1]) %in%
                 which(neib8[2] == abs1[, 2]))) {
          abs1 <- rbind(abs1, neib8)
          res <- eltest(ab = c(sortedx[idex1 + 1], sortedz[idez2 + 1]), x, y, z,
                        mu, eps)
          neg2_llr <- c(neg2_llr, res$'-2LLR')
        }
      }
      if (min(neg2_llr) == best) {
        return(best)
        ite <- FALSE
      }
      else {
        best <- min(neg2_llr)
        idex1 <- abs1[which(neg2_llr == best), 1]
        idez2 <- abs1[which(neg2_llr == best), 2]
        ite <- TRUE
      }
    }
  } else {
    return(NA)
  }
}

#' @export
PEL <- function(P1, P3, P2, x, y, z, eps = 0.01){
  res <- neighb_xyz(sp12 = c(P1, P3), x = x, y = y, z = z, true = P2, eps = eps)
  return(res)
}

### --- Adjusted EL method ----
myfun5_adj <- function(x, theta, mu, adj, eps){
  # theta = threshold
  u1 <- myfun5(x, theta, eps) - mu
  res <- c(u1, -adj*mean(u1))
  return (res)
}

meany_adj <- function(y, c1, c2, mu, adj, eps){
  vy <- meany(y, c1, c2, eps) - mu
  vy <- c(vy, -adj*mean(vy))
  return(vy)
}

eltest_adj <- function(ab = vector("numeric", 2), x.sample, y.sample, z.sample,
                       mu, adj_vec, eps) {
  # mu=c(P1,P3,P2)
  # ab -vector of quantiles, mu-vector of means, i.e. proportions
  all1 <- el.test(myfun5_adj(x.sample, ab[1], mu[1], adj_vec[1], eps), mu = 0)
  all3 <- el.test(myfun5_adj(z.sample, ab[2], 1 - mu[2], adj_vec[3], eps),
                  mu = 0)
  all2 <- el.test(meany_adj(y.sample, ab[1], ab[2], mu[3], adj_vec[2], eps),
                  mu = 0)
  list ('-2LLR' = all1$'-2LLR' + all2$'-2LLR' + all3$'-2LLR')
}

neighb_xyz_adj <- function(sp12 = vector("numeric", 2), x, y, z, true, eps){
  # true = P2
  # sp12 = c(P1, P3)
  n1 <- length(x)
  n2 <- length(y)
  n3 <- length(z)
  # find the threshold
  sortx_w <- myWdataclean2(z = x)
  sortedx <- sortx_w$value
  cx <- cumsum(sortx_w$weight/n1)
  idex1 <- ifelse(sp12[1] == 0, 1, max(which (cx <= sp12[1])))
  sortz_w <- myWdataclean2(z = z)
  sortedz <- sortz_w$value
  cz <- cumsum(sortz_w$weight/n3)
  idez2 <- ifelse(sp12[2] == 1, 1, max(which(cz <= 1 - sp12[2])))
  t1_est <- sortedx[idex1]
  t2_est <- sortedz[idez2]
  adj_par <- pmax(1, c(0.5*log(n1), 0.5*log(n2), 0.5*log(n3)))
  if(t1_est < t2_est){
    mu <- c(sp12, true) #P1,P3 and P2
    abs1 <- data.frame(x1 = idex1, z2 = idez2)
    res <- eltest_adj(c(t1_est, t2_est), x, y, z, mu, adj_vec = adj_par, eps)
    best <- res$'-2LLR'
    neg2_llr <- as.vector(best)
    ##############
    ite <- TRUE
    while (ite == TRUE) {
      neib1 <- c(idex1 - 1, idez2)
      if ((idex1 - 1 > 0) & (sortedx[idex1 - 1] < sortedz[idez2])) {
        if((!any(abs1[, 1] == neib1[1])) ||
           (!any(which(neib1[1] == abs1[, 1]) %in%
                 which(neib1[2] == abs1[, 2])))) {
          abs1 <- rbind(abs1, neib1)
          res <- eltest_adj(ab = c(sortedx[idex1 - 1], sortedz[idez2]),x, y, z,
                            mu, adj_par, eps)
          neg2_llr <- c(neg2_llr, res$'-2LLR')
        }
      }
      neib2 <- c(idex1, idez2 - 1)
      if ((sortedz[idez2 - 1] > sortedx[idex1]) & (idez2 - 1 > 0)) {
        if ((!any(abs1[, 2] == neib2[2])) ||
            (!any(which(neib2[1] == abs1[, 1]) %in%
                  which(neib2[2] == abs1[, 2])))) {
          abs1 <- rbind(abs1, neib2)
          res <- eltest_adj(ab = c(sortedx[idex1], sortedz[idez2 - 1]), x, y, z,
                            mu, adj_par, eps)
          neg2_llr <- c(neg2_llr, res$'-2LLR')
        }
      }
      neib3 <- c(idex1 + 1, idez2)
      if ((sortedx[idex1 + 1] < sortedz[idez2]) & (idex1 + 1 <= n1)) {
        if ((!any(abs1[, 1] == neib3[1])) ||
            (!any(which(neib3[1] == abs1[, 1]) %in%
                  which(neib3[2] == abs1[, 2])))) {
          abs1 <- rbind(abs1, neib3)
          res <- eltest_adj(ab = c(sortedx[idex1 + 1], sortedz[idez2]), x, y, z,
                            mu, adj_par, eps)
          neg2_llr <- c(neg2_llr, res $'-2LLR')
        }
      }
      neib4 <- c(idex1, idez2 + 1)
      if ((idez2 + 1 <= n3) & (sortedx[idex1] < sortedz[idez2 + 1])) {
        if ((!any(abs1[, 2] == neib4[2])) ||
            (!any(which(neib4[1] == abs1[, 1]) %in%
                  which(neib4[2] == abs1[, 2])))) {
          abs1 <- rbind(abs1, neib4)
          res <- eltest_adj(ab = c(sortedx[idex1], sortedz[idez2 + 1]), x, y, z,
                            mu, adj_par, eps)
          neg2_llr <- c(neg2_llr, res$'-2LLR')
        }
      }
      neib5 <- c(idex1 - 1, idez2 - 1)
      if ((idex1 - 1 > 0) & (idez2 - 1 > 0) &
          (sortedz[idez2 - 1] > sortedx[idex1 - 1])) {
        if (!any(which(neib5[1] == abs1[, 1]) %in%
                 which(neib5[2] == abs1[, 2]))) {
          abs1 <- rbind(abs1, neib5)
          res <- eltest_adj(ab = c(sortedx[idex1 - 1], sortedz[idez2 - 1]), x, y,
                            z, mu, adj_par, eps)
          neg2_llr <- c(neg2_llr, res$'-2LLR')
        }
      }
      neib6 <- c(idex1 - 1, idez2 + 1)
      if ((idex1 - 1 > 0) & (idez2 + 1 <= n3) &
          (sortedx[idex1 - 1] < sortedz[idez2 + 1])) {
        if (!any(which(neib6[1] == abs1[, 1]) %in%
                 which(neib6[2] == abs1[, 2]))) {
          abs1 <- rbind(abs1, neib6)
          res <- eltest_adj(ab = c(sortedx[idex1 - 1], sortedz[idez2 + 1]), x, y,
                            z, mu, adj_par, eps)
          neg2_llr <- c(neg2_llr, res $'-2LLR')
        }
      }
      neib7 <- c(idex1 + 1, idez2 - 1)
      if ((sortedx[idex1 + 1] < sortedz[idez2 - 1]) &
          (idex1 + 1 <= n1) & (idez2 - 1 > 0)) {
        if (!any(which(neib7[1] == abs1[, 1]) %in%
                 which(neib7[2] == abs1[, 2]))) {
          abs1 <- rbind(abs1, neib7)
          res <- eltest_adj(ab = c(sortedx[idex1 + 1], sortedz[idez2 - 1]), x, y,
                            z, mu, adj_par, eps)
          neg2_llr <- c(neg2_llr, res $'-2LLR')
        }
      }
      neib8 <- c(idex1 + 1, idez2 + 1)
      if ((sortedx[idex1 + 1] < sortedz[idez2 + 1]) &
          (idez2 + 1 <= n3) & (idex1 + 1 <= n1)) {
        if (!any(which(neib8[1] == abs1[, 1]) %in%
                 which(neib8[2] == abs1[, 2]))) {
          abs1 <- rbind(abs1, neib8)
          res <- eltest_adj(ab = c(sortedx[idex1 + 1], sortedz[idez2 + 1]), x, y,
                            z, mu, adj_par, eps)
          neg2_llr <- c(neg2_llr, res$'-2LLR')
        }
      }
      if (min(neg2_llr) == best) {
        return(best)
        ite <- FALSE
      }
      else {
        best <- min(neg2_llr)
        idex1 <- abs1[which(neg2_llr == best), 1]
        idez2 <- abs1[which(neg2_llr == best), 2]
        ite <- TRUE
      }
    }
  } else {
    return(NA)
  }
}

#' @export
AEL <- function(P1, P3, P2, x, y, z, eps = 0.01){
  res <- neighb_xyz_adj(sp12 = c(P1, P3), x = x, y = y, z = z, true = P2,
                        eps = eps)
  return(res)
}
