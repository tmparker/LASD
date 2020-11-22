# Show what the T_3 function should look like in the partially identified normal
# location model - all three distributions are normal.
# It is assumed that the mean of the control and policy A outcomes are zero, and
# the mean of policy B outcomes are varied.  All the scale parameters are set to
# 1.

library(Rcpp)
sourceCpp("../../lasd.cpp") # for the lohi function.

n <- 10000
muB <- 2.9 # the curves start to to above zero at approx. muB = 2.8.
gstep <- 0.1
con <- rnorm(n)
A <- rnorm(n)
B <- rnorm(n, mean = muB)
pool <- c(con, A, B)
big <- ceiling(max(abs(pool)))
grd <- seq(-big, big, by = gstep) # symmetric around 0, odd number of points.
half <- 0.5 * (length(grd) - 1) + 1

t3fun <- function(x, mB) {
  la.part <- 2 * pnorm(x / 2) - 1
  ub.minus <- ifelse(x < -mB, 1, 2 * pnorm((-x - mB) / 2))
  ub.plus <- ifelse(x < mB, 2 * pnorm((x - mB) / 2), 1)
  t3 <- la.part - ub.minus - ub.plus
  return(t3)
}

pdf(file = "pics_partial.pdf", width = 8, height = 4)
par(mfrow = c(1, 3), mar = c(3.1, 3.1, 1.1, 1.1), mgp = 2:0)
curve(t3fun(x, 2.7), from = 0, to = 6, col = "blue", main = 
      expression(mu[B] == 2.7), ylab = expression(T[3](G)(x)))
abline(h = 0, lty = 3, col = "gray")
curve(t3fun(x, 2.8), from = 0, to = 6, col = "blue", main = 
      expression(mu[B] == 2.8), ylab = expression(T[3](G)(x)))
abline(h = 0, lty = 3, col = "gray")
curve(t3fun(x, 2.9), from = 0, to = 6, col = "blue", main = 
      expression(mu[B] == 2.9), ylab = expression(T[3](G)(x)))
abline(h = 0, lty = 3, col = "gray")
dev.off()
