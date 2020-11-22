# This file demonstrates the use of the statistics derived from in the paper.
# These are loss averse-sensitive dominance test statistics that are either
# supremum norm or L2 norm statistics, and includes one example where the
# distributions are assumed not point identified and bounds are used instead.
# It also has regular first order stochastic dominance test statistics for
# comparison.

library(Rcpp)
sourceCpp("lasd.cpp")
set.seed(123)

n0 <- 100
nA <- 100
nB <- 100
R <- 100
alpha <- 0.05
tiny <- .Machine$double.eps^(1/2)
gstep <- 0.5

# The null hypothesis here is that A dominates B, which means that the test
# statistics used below should be nonpositive.  An unimaginative way to generate
# data that look like this is to let the distributions be identical members of a
# location-scale family, where the location parameter for B is at most that of
# the location of A.  You can generate pairs that violate this by raising the
# location of B above that of A.
# Generating distributions this way his allows us to benchmark the process with
# the FOSD tests used below.

mean0 <- 0
meanA <- -2.75
meanB <- 0
Con <- rnorm(n0, mean = mean0)
A <- rnorm(nA, mean = meanA)
B <- rnorm(nB, mean = meanB)

# Evaluate T1 & T2 processes on R_+.  teval() arranges the samples so that will
# happen the right way, and t1 & t2 generate processes.
tpool <- teval(A, B)
T1 <- t1(tpool, A, B)
T2 <- t2(tpool, A, B)
pool <- deval(A, B)
dpro <- diffedf(pool, A, B)

nAB <- nA + nB
# the next arguments are scaled up by sqrt(n) compared to the paper's notation.
an <- 4 * log(log(nAB)) # for contact set (L2)
bn <- sqrt(log(log(nAB))) # for epsilon-maximization (sup)
cn <- bn # for deciding whether sup m1 > sup m2 or other way (sup)

# Supremum and L2 statistics using each of the three processes.  l2n2
# is different than l2n1 because the T2 process is 2-dimensional.
stats <- c(max(T1), l2n1(tpool, T1), max(T2), l2n2(tpool, T2),
           max(dpro), l2n1(pool, dpro))

# Bootstrap
bT2 <- LASDboot(A, B, R, bn, an, cn)
bFOSD <- FOSDboot(A, B, R, lam)
pvals <- c(bT2$pval, bFOSD$pval)
names(pvals) <- c("sup_T1", "L2_T1", "sup_T2", "L2_T2", "sup_FOSD", "L2_FOSD")

cat("point-identified case: \n")
cat("p-values are: \n")
print(pvals)

### Test with interval-identified criterion
n3 <- n0 + nA + nB
# Also scaled by sqrt(n3) compared to the paper's notation.
an <- 4 * log(log(n3)) # for contact set estimation
bn <- sqrt(log(log(n3))) # for eps-max (inside)
dn <- bn # eps-max (outside)

b3 <- LASDboot_partial(Con, A, B, R, bn, an, dn, gstep)
pval_T3 <- b3$pval
names(pval_T3) <- c("sup_T3", "L2_T3")

cat("interval-identified case: \n")
cat("p-values are: \n")
print(pval_T3)

