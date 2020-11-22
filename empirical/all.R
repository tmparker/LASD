# L2 and supremum norms applied to the T2 process, and L2 or supremum norms
# applied to the difference process used for FOSD tests

library(Hmisc)
library(Rcpp)
sourceCpp("../../lasd.cpp", cacheDir="../../lib_cpp")
bootrep <- 1999
set.seed(321)

# For running locally
prefix <- "../data_copy/"
pre.name <- c("levels/pre_avg_totincome_contr.txt", "levels/avg_tot_pretl_contr.txt")
post.name <- c("levels/post_avg_total_income_contr.txt", "levels/avg_tot_posttl_contr.txt")
chg.name <- c("changes/ch_tot_contr.txt", "changes/ch_tot_tl_contr.txt")

# Tests of LASD in changes, null is that JF dominates AFDC.
chg <- vector(mode = "list", length = 2)
for (j in 1:2) {
  con <- unlist(read.table(paste0(prefix, chg.name[j]), skip = 1), use.names = FALSE)
  tre.name <- sub("contr", "treat", chg.name[j])
  tre <- unlist(read.table(paste0(prefix, tre.name), skip = 1), use.names = FALSE)
  n <- length(con) + length(tre)
  an <- 4 * log(log(n))
  bn <- sqrt(log(log(n)))
  cn <- bn
  chg[[j]] <- LASDboot(tre, con, bootrep, bn, an, cn)
  chg[[j]] <- lapply(chg[[j]], function(x) {x[c(2, 4)]}) # only select V2 and W2
}
chg_vec <- c(sapply(chg, function(x) do.call(rbind, x)))

# Tests of FOSD in levels, null is that JF dominates AFDC.
lev <- vector(mode = "list", length = 2)
for (j in 1:2) {
  con <- unlist(read.table(paste0(prefix, post.name[j]), skip = 1), use.names = FALSE)
  tre.name <- sub("contr", "treat", post.name[j])
  tre <- unlist(read.table(paste0(prefix, tre.name), skip = 1), use.names = FALSE)
  n <- length(con) + length(tre)
  an <- 4 * log(log(n))
  lev[[j]] <- FOSDboot(tre, con, bootrep, an)
}
lev_vec <- c(sapply(lev, function(x) do.call(rbind, x)))

# Create a table.
tabdat <- cbind(chg_vec, lev_vec)
timename <- c("Before JF time limit", "After JF time limit")
testpv <- vector(mode = "character", length = 8)
testpv[1:4*2-1] <- rep(c("supremum norm", "L2 norm"), 2)
testpv[1:2*4] <- "p-value"
hypname <- c("$F_{JF} \\succeq F_{AFDC}$",
              "$G_{JF} \\succeq G_{AFDC}$")

tab <- format.df(tabdat, dec = 4, na.blank = TRUE)
tst <- latex(tab, file = "all_tests.tex", label = "tab:all_pointIDtests",
      rowname = testpv, title = '', ctable = TRUE, math.col.names = TRUE,
      col.just = rep("c", 2), cgroup = c("LASD", "FOSD"),
      n.cgroup = c(1, 1), colheads = hypname, rgroup = timename, 
      n.rgroup = c(4, 4), caption.loc = "bottom",
      caption = "Supremum and $L_2$ tests for inferring whether the Jobs First 
      (JF) program would be preferred to the Aid to Families with Dependent 
      Children (AFDC).  The first column uses changes in income and the second 
      column measures income in levels without regard to pre-policy income.  
      Both types of test statistic agree on rejection decisions.  1999 
      bootstrap repetitions used in each test.")

