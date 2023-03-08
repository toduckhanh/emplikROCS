## Generate random variable from a mixture-normal model

#' @export
rnorm_mixture <- function(n, mean, sd, prob) {
  u <- runif(n)
  u_prob <- c(0, cumsum(prob))
  x <- numeric(n)
  for (i in 1:n) {
    for (j in 2:length(u_prob)) {
      if (u[i] >= u_prob[j - 1] & u[i] <= u_prob[j]) {
        x[i] <- rnorm(1, mean[j - 1], sd[j - 1])
      }
    }
  }
  return(x)
}
