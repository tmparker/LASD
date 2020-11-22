# This file produces the tables and figures used in the empirical section.

library(Hmisc)
library(Rcpp)
sourceCpp("../lasd.cpp", cacheDir="../lib_cpp")
bootrep <- 1999
set.seed(321)

# Just get names from control groups, then substitute in treatment labels for
# control labels below.
prefix <- "./data_copy/"
pre.name <- c("levels/pre_avg_totincome_contr.txt", "levels/avg_tot_pretl_contr.txt")
post.name <- c("levels/post_avg_total_income_contr.txt", "levels/avg_tot_posttl_contr.txt")
chg.name <- c("changes/ch_tot_contr.txt", "changes/ch_tot_tl_contr.txt")

#-----
# Tests of LASD in changes, null is that JF dominates AFDC.
chg <- matrix(0, 2, 2)
for (j in 1:2) {
  con <- unlist(read.table(paste0(prefix, chg.name[j]), skip = 1), use.names = FALSE)
  tre.name <- sub("contr", "treat", chg.name[j])
  tre <- unlist(read.table(paste0(prefix, tre.name), skip = 1), use.names = FALSE)
  ev <- teval(con, tre)
  t2proc <- t2(ev, tre, con)
  chg[j, 1] <- l2n2(ev, t2proc)
  n <- length(con) + length(tre)
  an <- 4 * log(log(n))
  bn <- sqrt(log(log(n)))
  cn <- bn
  chg[j, 2] <- LASDboot(tre, con, bootrep, bn, an, cn)$pval[4]
}

#-----
# FOSD tests
# Tests of stochastic dominance in levels.
lev <- matrix(0, 2, 2)
for (j in 1:2) {
  con <- unlist(read.table(paste0(prefix, post.name[j]), skip = 1), use.names = FALSE)
  tre.name <- sub("contr", "treat", post.name[j])
  tre <- unlist(read.table(paste0(prefix, tre.name), skip = 1), use.names = FALSE)
  dev <- deval(con, tre)
  fosd <- diffedf(dev, tre, con)
  n <- length(con) + length(tre)
  an <- 4 * log(log(n))
  lev[j, 1] <- l2n1(dev, fosd)
  lev[j, 2] <- FOSDboot(tre, con, bootrep, an)$pval[2]
}

# Test equality of the levels CDFs using post-TL averages.
con <- unlist(read.table(paste0(prefix, post.name[2]), skip = 1), use.names = FALSE)
tre.name <- sub("contr", "treat", post.name[2])
tre <- unlist(read.table(paste0(prefix, tre.name), skip = 1), use.names = FALSE)
lev.eq.test <- eqboot(con, tre, bootrep)
lev.eq.ans <- c(lev.eq.test$stats[2], lev.eq.test$pval[2])

#-----
# Create a table.
tabdat <- cbind(c(t(chg)), c(lev[1, ], lev[2, ]))#, c(rep(NA, 2), lev.eq.ans))
timename <- c("Before JF time limit", "After JF time limit")
timepv <- vector(mode = "character", length = 4)
timepv[1:2*2-1] <- timename
timepv[1:2*2] <- "p-value"
hypname <- c("$F_{JF} \\succeq F_{AFDC}$",
              "$G_{JF} \\succeq G_{AFDC}$")#, "equality")

#load("emp_pointID.RData")

tab <- format.df(tabdat, dec = 4, na.blank = TRUE)
tst <- latex(tab, file = "table_of_tests.tex", label = "tab:pointIDtests",
      rowname = timepv, title = '', ctable = TRUE, math.col.names = TRUE,
      col.just = rep("c", 3), cgroup = c("LASD", "FOSD"),
      n.cgroup = c(1, 1), colheads = hypname, caption.loc = "bottom",
      caption = "Tests for inferring whether the Jobs First (JF) program would
      be preferred to the Aid to Families with Dependent Children (AFDC).
      Column titles paraphrase the null hypotheses in the tests.  The first 
      column uses changes in income and the second column measures income in 
      levels without regard to pre-policy income.  1999 bootstrap repetitions 
      used in each test.")

# Show this in two pictures.

# LASD picture.
pdf(file = "LASD_pic.pdf", width = 6, height = 2)
par(mfrow = c(1, 3), mar = c(3.1, 3.3, 1.1, 1.1), mgp = 2:0)
for (i in 2:2) {
  con <- unlist(read.table(paste0(prefix, chg.name[i]), skip = 1),
                use.names = FALSE)
  tre.name <- sub("contr", "treat", chg.name[i])
  tre <- unlist(read.table(paste0(prefix, tre.name), skip = 1),
                use.names = FALSE)
  ev <- teval(con, tre)
  t2proc <- t2(ev, tre, con)
  Ftre <- ecdf(tre)
  Fcon <- ecdf(con)
  plot(tre, Ftre(tre), type = "l", lty = 1, main = "empirical DFs (changes)",
       xlab = "gain/loss", xlim = range(c(con, tre)),
       ylab = "empirical DFs (changes)", font.main = 1)
  lines(con, Fcon(con), lty = 2, col = "black")
  legend("topleft", c("JF", "AFDC"), lty = 1:2, bty = "n")
  plot(ev, t2proc[1, ], type = "l", xlab = "absolute gain/loss",
        ylab = expression(paste("scaled ", hat(m)[1])),
        main = expression(paste(m[1], " process")))
  abline(h = 0, lty = 4, col = "gray")
  plot(ev, t2proc[2, ], type = "l", xlab = "absolute gain/loss",
        ylab = expression(paste("scaled ", hat(m)[2])),
        main = expression(paste(m[2], " process")))
  abline(h = 0, lty = 4, col = "gray")
}
dev.off()

# FOSD picture.
i <- 2
pdf(file = "FOSD_pic.pdf", width = 6, height = 3)
par(mfrow = c(1, 2), mar = c(3.1, 3.1, 1.1, 1.1), mgp = 2:0)
con <- unlist(read.table(paste0(prefix, post.name[i]), skip = 1),
              use.names = FALSE)
tre.name <- sub("contr", "treat", post.name[i])
tre <- unlist(read.table(paste0(prefix, tre.name), skip = 1),
              use.names = FALSE)
dev <- deval(con, tre)
fosdTdomC <- diffedf(dev, tre, con)
Ftre <- ecdf(tre)
Fcon <- ecdf(con)
plot(tre, Ftre(tre), type = "l", lty = 1, main = "empirical DFs (levels)",
      xlab = "income", xlim = range(c(con, tre)),
      ylab = "empirical DFs", font.main = 1)
lines(con, Fcon(con), lty = 2, col = "black")
legend("topleft", c("JF", "AFDC"), lty = 1:2, bty = "n")
plot(dev, fosdTdomC, type = "l", main = "FOSD process", xlab = "income",
      ylab = "scaled difference of DFs", font.main = 1)
abline(h = 0, lty = 4, col = "gray")
dev.off()

