# Explaining the triangular distributions

# Simply using uniform distributions would also work, but would be more boring -
# for those distributions, the dominating alternative would have to dominate in
# a way that FOSD would also order.  It is possible to find triangular
# alternatives that respect LASD but not FOSD because they cross.

dtri <- function(x, lo, mid, hi) {
  d <- ifelse(x < lo, 0, ifelse(lo <= x & x < mid, 
              2 * (x - lo) / ((hi - lo) * (mid - lo)), 
              ifelse(mid <= x & x < hi, 
                      2 * (hi - x) / ((hi - lo) * (hi - mid)), 0)))
  return(d)
}

ptri <- function(q, lo, mid, hi) {
  p <- ifelse(q < lo, 0, ifelse(lo <= q & q < mid, 
              (q - lo)^2 / ((hi - lo) * (mid - lo)),
              ifelse(mid <= q & q < hi, 
                      1 - (hi - q)^2 / ((hi - lo) * (hi - mid)), 1)))
  return(p)
}

rtri <- function(n, lo, mid, hi) {
  Fmid <- (mid - lo) / (hi - lo)
  U <- runif(n)
  sam <- ifelse(U < Fmid, lo + sqrt(U * (hi - lo) * (mid - lo)),
                hi - sqrt((1 - U) * (hi - lo) * (hi - mid)))
  return(sam)
}

#n <- 1000
#L <- -1
#M <- 0
#H <- 1
#sam <- rtri(n, L, M, H)
#
#par(mfrow = c(2, 1), mar = c(3.1, 3.1, 1.1, 1.1), mgp = 2:0)
#hist(sam, freq = FALSE)
#curve(dtri(x, L, M, H), add = TRUE, col = "blue")
#plot(ecdf(sam), do.points = FALSE)
#curve(ptri(x, L, M, H), add = TRUE, col = "blue")

# Here is how alternative distributions look.
loA <- -1
midA <- 0
hiA <- 1

eps <- 0.25
loB <- loA - eps
midB <- midA - eps
hiB <- hiA + eps
# makes the mean of the B distribution -eps / 3.
cvec <- c("dodgerblue3", "firebrick")

pdf(file = "pics_triangle.pdf", width = 6, height = 3)
par(mfrow = c(1, 2), mar = c(3.1, 3.1, 1.1, 1.1), mgp = 2:0)
par(mfrow = c(1, 2), mar = c(3.1, 3.1, 1.1, 1.1), mgp = 2:0)
curve(dtri(x, loA, midA, hiA), from = -2, to = 2, col = cvec[1], lty = 1,
      ylab = expression(f(x)), main = "Densities")
curve(dtri(x, loB, midB, hiB), add = TRUE, col = cvec[2], lty = 2)
abline(v = 0, lty = 3, col = "gray")
legend("topleft", legend = c(expression(f[A]), expression(f[B])), lty = 1:2,
        col = cvec, bty = "n")

curve(ptri(x, loA, midA, hiA), from = -2, to = 2, col = cvec[1], lty = 1,
      ylab = expression(F(x)), main = "CDFs")
curve(ptri(x, loB, midB, hiB), add = TRUE, col = cvec[2], lty = 2)
abline(v = 0, lty = 3, col = "gray")
legend("topleft", legend = c(expression(F[A]), expression(F[B])), lty = 1:2,
        col = cvec, bty = "n")
dev.off()

#t2fun <- function(x, lA, mA, hA, lB, mB, hB) {
#  t21 <- ptri(-x, lA, mA, hA) - ptri(-x, lB, mB, hB)
#  t22 <- ptri(-x, lA, mA, hA) - ptri(-x, lB, mB, hB) + 
#          ptri(x, lA, mA, hA) - ptri(x, lB, mB, hB)
#  return(list(t21 = t21, t22 = t22))
#}
#
#curve(t2fun(x, loA, midA, hiA, loB, midB, hiB)$t21, from = 0, to = 2)
#abline(h = 0, lty = 3, col = "gray")
#curve(t2fun(x, loA, midA, hiA, loB, midB, hiB)$t22, from = 0, to = 2)
#abline(h = 0, lty = 3, col = "gray")


