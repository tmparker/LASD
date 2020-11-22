# Empirical rejection probabilities for the partially-identified test statistics.

library(Rcpp)
sourceCpp("../lasd.cpp", cacheDir="../lib_cpp")

library(rlecuyer)
nstream <- 100 # number of parallel random number streams
.lec.SetPackageSeed(6:1)
stream.names <- paste0("st", 1:nstream)
.lec.CreateStream(stream.names)
# To see all of them, they are stored in a global variable called 
# .lec.Random.seed.table.
runname <- as.numeric(commandArgs())
runname <- runname[length(runname)] # for some reason it is a vector.
.lec.CurrentStream(paste0("st", runname))

ssize <- c(100, 500, 1000)
R <- c(499, 999, 1999)
alpha <- 0.05
reps <- 10
gstep <- 0.05
mbseq <- seq(-1, 20, by = 1)

mc_partialID <- function(i, ss, m0, mA, mB, R, gstep) {
  con <- rnorm(ss, mean = m0)
  A <- rnorm(ss, mean = mA)
  B <- rnorm(ss, mean = mB)
  rn <- 3 * ss
  an <- 4 * log(log(rn)) # for contact set (L2)
  bn <- sqrt(log(log(rn))) # for epsilon-maximization (sup)
  dn <- bn # for outer epsilon-maximizer set (M^nec, sup-norm stat).
  bT3 <- LASDboot_partial(con, A, B, R, bn, an, dn, gstep)
  pvals <- bT3$pval
  pvals
}

partID <- function(ssize, R, gstep, alpha) {
  meanB <- 2.8 + mbseq / sqrt(ssize) # 2.8 is approx. where to start rejecting
  rej.prob <- matrix(0, 2, length(meanB))
  for (i in 1:length(meanB)) {
    mc <- sapply(1:reps, mc_partialID, ss = ssize, 
                  m0 = 0, mA = 0, mB = meanB[i], R = R, gstep = gstep)
    rej.prob[, i] <- rowMeans((mc < alpha))
  }
  dimnames(rej.prob) <- list(c("sup_T3", "L2_T3"), mbseq)
  rej.prob
}

# coordinates of this list correspond to ssize and R vectors.
results.list <- list(length = length(ssize)) 
for (k in 1:length(ssize)) {
  results.list[[k]] <- partID(ssize[k], R[k], gstep, alpha)
}

.lec.CurrentStreamEnd()
save(results.list, file = paste0("./ind_runs/run", runname, ".rda"))

