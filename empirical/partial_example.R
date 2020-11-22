# This is like tables_and_figures.R but without assuming that identification is
# exact.  All samples, pre-RA and post-RA under different policies, are treated
# as if we did not know the observations were longitudinal, as a demonstration
# of the partially-identified inference procedure.

library(Hmisc)
library(Rcpp)
sourceCpp("../lasd.cpp", cacheDir="../lib_cpp")
bootrep <- 1999
set.seed(321)

# For running locally
prefix <- "./data_copy/"
pre.name <- c("levels/pre_avg_totincome_contr.txt", "levels/avg_tot_pretl_contr.txt")
lev.name <- c("levels/post_avg_total_income_contr.txt", "levels/avg_tot_posttl_contr.txt")

#-----
# LASD tests
# Tests of LASD in changes.
ans <- ans.rev <- vector(mode = "list", length = 2)
for (j in 1:2) {
  ans[[j]] <- ans.rev[[j]] <- matrix(0, 2, 2)
  AFDC.pre <- unlist(read.table(paste0(prefix, pre.name[j]), skip = 1), use.names = FALSE)
  pre.JF.name <- sub("contr", "treat", pre.name[j])
  JF.pre <- unlist(read.table(paste0(prefix, pre.JF.name), skip = 1), use.names = FALSE)
  pre <- c(AFDC.pre, JF.pre)
  AFDC.post <- unlist(read.table(paste0(prefix, lev.name[j]), skip = 1), use.names = FALSE)
  post.JF.name <- sub("contr", "treat", lev.name[j])
  JF.post <- unlist(read.table(paste0(prefix, post.JF.name), skip = 1), use.names = FALSE)
  gstep <- 0.01 # The observations are very closely spaced between 0 and around 10-11.
  n <- length(pre) + length(AFDC.post) + length(JF.post)
  Kap <- sqrt(log(log(n))) # for epsilon-maximization (sup)
  Lam <- 4 * log(log(n)) # for contact set (L2)
  Mu <- Kap # for outer epsilon-maximizer set (M^nec, sup-norm stat).
  tst <- LASDboot_partial(pre, AFDC.post, JF.post, bootrep, Kap, Lam, Mu, gstep)
  tst.rev <- LASDboot_partial(pre, JF.post, AFDC.post, bootrep, Kap, Lam, Mu, gstep)
  ans[[j]][1, ] <- tst$stats
  ans[[j]][2, ] <- tst$pval
  ans.rev[[j]][1, ] <- tst.rev$stats
  ans.rev[[j]][2, ] <- tst.rev$pval
}

mat <- cbind(t(do.call(cbind, ans)), t(do.call(cbind, ans.rev)))
timename <- c("around RA", "around TL")
timepv <- rep(c("statistic", "p-value"), 2)
hypname <- c("$F_{AFDC} \\succeq F_{JF}$", "$F_{JF} \\succeq F_{AFDC}$")
testname <- rep(c("$V_{3n}$", "$W_{3n}$"), 2)
dimnames(mat) <- list(testname, timepv)

results <- latex(mat, file = "partial_table.tex", label = "tab:partialIDtests", 
      rowname = timepv, title = '', ctable = TRUE, math.col.names = TRUE, 
      col.just = rep("c", 4), cgroup = hypname,
      n.cgroup = c(2, 2), caption.loc = "bottom", caption = "Results of dominance tests when no identification of income changes is assumed in the empirical example.  We cannot reject the hypothesis that either policy dominates the other.  1,999 bootstrap repetitions used.")

