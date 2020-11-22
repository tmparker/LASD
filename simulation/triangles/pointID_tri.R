# Empirical rejection probabilities for the point-identified test statistics
# using triangular distributions.

library(Rcpp)
sourceCpp("../lasd.cpp", cacheDir="../lib_cpp")
#set.seed(123)

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
eseq <- seq(-0.5, 0.5, by = 0.1)

rtri <- function(n, lo, mid, hi) {
  Fmid <- (mid - lo) / (hi - lo)
  U <- runif(n)
  sam <- ifelse(U < Fmid, lo + sqrt(U * (hi - lo) * (mid - lo)),
                hi - sqrt((1 - U) * (hi - lo) * (hi - mid)))
  return(sam)
}

mc_pointID_tri <- function(i, ss, eps, R) {
  A <- rtri(ss, -1, 0, 1)
  B <- rtri(ss, -1 - eps, -eps, 1 + eps)
  rn <- 2 * ss
  an <- 4 * log(log(rn)) # for contact set (L2)
  bn <- sqrt(log(log(rn))) # for epsilon-maximization (sup)
  cn <- bn # for deciding sup m1 > sup m2 or opposite
  bT2 <- LASDboot(A, B, R, bn, an, cn)
  bFOSD <- FOSDboot(A, B, R, an)
  bequal <- eqboot(A, B, R)
  pvals <- c(bT2$pval, bFOSD$pval, bequal$pval)
  pvals
}

pointID <- function(ssize, R, alpha) {
  epsB <- eseq #/ sqrt(ssize) # easier not to use local alternatives
  rej.prob <- matrix(0, 8, length(epsB))
  for (i in 1:length(epsB)) {
    mc <- sapply(1:reps, mc_pointID_tri, ss = ssize, eps = epsB[i], R = R)
    rej.prob[, i] <- rowMeans((mc < alpha))
  }
  dimnames(rej.prob) <- list(c("sup_T1", "sup_T2", "L2_T1", "L2_T2", 
                                "sup_FOSD", "L2_FOSD", "sup-eq", "L2-eq"), 
                              eseq)
  rej.prob
}

# coordinates of this list correspond to ssize and R vectors.
results.list <- list(length = length(ssize)) 
for (k in 1:length(ssize)) {
  results.list[[k]] <- pointID(ssize[k], R[k], alpha)
}

.lec.CurrentStreamEnd()
save(results.list, file = paste0("./ind_runs/run", runname, ".rda"))
#save(results.list, file = "triangles.rda")

