# R code for the tests proposed in the paper.
# This translates the code in pdom.cpp into pure R code.

tevalR <- function(A, B) {
  tiny = .Machine$double.eps^(1/2)
  ev <- c(0, A, B)
  ev[ev < 0] <- -ev[ev < 0] + tiny
  ev <- sort(unique(ev))
  ev
}

devalR <- function(A, B) {
  ev <- sort(unique(c(A, B)))
  ev
}

diffecdfR <- function(ev, A, B) {
  scl <- sqrt(length(A) + length(B))
  Fa <- ecdf(A)
  Fb <- ecdf(B)
  diffproc <- scl * (Fa(ev) - Fb(ev))
  diffproc
}

l2n1R <- function(ev, fun) {
  ne <- length(ev)
  f2 <- pmax(0, fun)^2
  ans <- sqrt(sum(f2[-ne] * diff(ev)))
  ans
}

# For this function, fun should be a 2xne matrix, where ne is the length of ev.
l2n2R <- function(ev, fun) {
  ne <- length(ev)
  f2 <- matrix(pmax(0, fun)^2, nrow = 2)
  ans <- sqrt(sum(t(f2[, -ne]) * diff(ev)))
  ans
}

t1R <- function(ev, A, B) {
  scl <- sqrt(length(A) + length(B))
  Fa <- ecdf(A)
  Fb <- ecdf(B)
  na <- Fa(-ev)
  nb <- Fb(-ev)
  pa <- Fa(ev)
  pb <- Fb(ev)
  pos <- pmax(0, pa - pb)
  f <- scl * (pos + na - nb)
  f
}

t2R <- function(ev, A, B) {
  scl <- sqrt(length(A) + length(B))
  Fa <- ecdf(A)
  Fb <- ecdf(B)
  na <- Fa(-ev)
  nb <- Fb(-ev)
  pa <- Fa(ev)
  pb <- Fb(ev)
  f <- matrix(0, 2, length(ev))
  f[1, ] <- scl * (na - nb)
  f[2, ] <- scl * (pa - pb + na - nb)
  f
}

lasdbootR <- function(A, B, R, kap, lam, mu) {
  tiny = .Machine$double.eps^(1/2)
  ev <- tevalR(A, B)
  nev <- length(ev)
  T1 <- t1R(ev, A, B)
  T2 <- t2R(ev, A, B)
  stats <- c(max(T1), max(T2), l2n1R(ev, T1), l2n2R(ev, T2))
  maxm1 <- max(T2[1, ])
  maxm2 <- max(T2[2, ])
  mind <- 0
  if (maxm1 > maxm2 + mu) mind = 1
  if (maxm1 < maxm2 - mu) mind = 2
  epsmax21 <- (T2[1, ] >= maxm1 - kap)
  epsmax22 <- (T2[2, ] >= maxm2 - kap)
  epsmax2 <- rbind(epsmax21, epsmax22)
  con21 <- (abs(T2[1, ]) < lam)
  con22 <- (abs(T2[2, ]) < lam)
  if (!any(con21)) con21 <- rep(TRUE, ne)
  if (!any(con22)) con22 <- rep(TRUE, ne)
  conT2 <- rbind(con21, con22)
  t2B <- matrix(0, 2, R)
  for (i in 1:R) {
    sA <- sample(A, replace = TRUE)
    sB <- sample(B, replace = TRUE)
    f2b <- t2R(ev, sA, sB) - T2
    v21b <- max(0, max(f2b[1, ] * epsmax21))
    v22b <- max(0, max(f2b[2, ] * epsmax22))
    v2b <- max(v21b, v22b) * (mind == 0) + v21b * (mind == 1) + 
            v22b * (mind == 2)
    w2b <- l2n2(ev, f2b * conT2)
    t2B[, i] <- c(v2b, w2b)
  }
  rejv1 <- stats[1] < t2B[1,] + tiny
  rejv2 <- stats[2] < t2B[1,] + tiny
  rejw1 <- stats[3] < t2B[2,] + tiny
  rejw2 <- stats[4] < t2B[2,] + tiny
  pval <- c(mean(rejv1), mean(rejv2), mean(rejw1), mean(rejw2))
  ans <- list(stats = stats, pval = pval)
  ans
}

FOSDbootR <- function(A, B, R, lam) {
  ev <- devalR(A, B)
  dF <- diffecdfR(ev, A, B)
  conF <- (abs(dF) < lam)
  if (!any(conF)) conF <- rep(TRUE, length(ev))
  FSDB <- matrix(0, 2, R)
  for (i in 1:R) {
    sA <- sample(A, replace = TRUE)
    sB <- sample(B, replace = TRUE)
    bootF <- diffecdfR(ev, sA, sB)
    bootF <- (bootF - dF) * conF
    FSDB[1, i] <- max(bootF)
    FSDB[2, i] <- l2n1R(ev, bootF)
  }
  vn <- max(dF * conF)
  wn <- l2n1R(ev, dF * conF)
  rej <- c(vn, wn) < FSDB + tiny
  pval <- rowMeans(rej)
  ans <- list(dF=dF, conF=conF, FSDB=FSDB, stats=c(vn, wn), rej=rej, 
              pval=pval)
  ans
}

# Estimate lower and upper bounds for a pair of marginal distribution functions.
lohiR <- function(ev, tr, co) {
  ne <- length(ev)
  Ft <- ecdf(tr)
  Ftr <- Ft(ev)
  bounds <- matrix(0, 2, ne)
  for (i in 1:ne) {
    coshift <- co + ev[i]
    Fc <- ecdf(coshift)
    Fco <- Fc(ev)
    dpro <- c(0, Ftr - Fco)
    bounds[1, i] <- max(dpro)
    bounds[2, i] <- 1 + min(dpro)
  }
  bounds
}

gevalR <- function(con, A, B, gstep) {
  lims <- c(range(A) - rev(range(con)), range(B) - rev(range(con)))
  glen <- ceiling(max(abs(lims)) / gstep)
  gtop <- glen * gstep
  ulen <- 2 * glen + 1
  grd <- -gtop + gstep * 1:ulen
  grd
}

# Taking the lohiR function and wrapping it up into a T3 form.  This produces a
# T3 process, and later you can use it to calculate statistics.
t3R <- function(ev, con, A, B) {
  scl <- sqrt(length(con) + length(A) + length(B))
  epos <- ev[ev >= 0]
  nLA <- lohiR(-epos, A, con)[1, ]
  nUB <- lohiR(-epos, B, con)[2, ]
  pLA <- lohiR(epos, A, con)[1, ]
  pUB <- lohiR(epos, B, con)[2, ]
  f <- scl * (pLA - pUB + nLA - nUB)
  f
}

# Make estimates that are needed for bootstrapping the bounding functions
# consistently.  In particular, sets of marginal epsilon-maximizers and
# minimizers are needed, and sets that denote where maps of the A- and
# B-empirical processes (mapped through the marginal pointwise maximizer or
# minimizer) are above zero, below zero or close to zero.

bndestR <- function(ev, con, A, B, kap, lam, mu) {
  scl <- sqrt(length(con) + length(A) + length(B))
  epos <- ev[ev >= 0]
  ne <- length(ev)
  nx <- (ne - 1) / 2 + 1
  gA <- ecdf(A)
  gB <- ecdf(B)
  GA <- gA(ev)
  GB <- gB(ev)
  m1 <- m2 <- m3 <- m4 <- matrix(0, ne, nx)
  mm1 <- mm2 <- mm3 <- mm4 <- double(nx)
  M1 <- M2 <- M3 <- M4 <- matrix(0, ne, nx)
  for (i in 1:nx) {
    coshift1 <- con - epos[i]
    coshift2 <- con + epos[i]
    Gco1 <- ecdf(coshift1)
    Gco2 <- ecdf(coshift2)
    m1[, i] <- GA - Gco1(ev)
    m2[, i] <- GA - Gco2(ev)
    m3[, i] <- Gco1(ev) - GB
    m4[, i] <- Gco2(ev) - GB
    mm1[i] <- max(m1[, i])
    mm2[i] <- max(m2[, i])
    mm3[i] <- max(m3[, i])
    mm4[i] <- max(m4[, i])
    M1[, i] <- (m1[, i] >= mm1[i] - kap)
    M2[, i] <- (m2[, i] >= mm2[i] - kap)
    M3[, i] <- (m3[, i] >= mm3[i] - kap)
    M4[, i] <- (m4[, i] >= mm4[i] - kap)
  }
  t3pro <- mm1 + mm2 + mm3 + mm4 - 2
  csest <- (abs(t3pro) <= lam)
  if (!any(csest)) csest <- rep(TRUE, ne)
  Mnec <- (t3pro >= max(t3pro) - mu)
  ans <- list(m1 = m1, M1 = M1, M2 = M2, M3 = M3, M4 = M4, csest = csest, Mnec = Mnec)
  ans
}

LASDboot_partialR <- function(con, A, B, bounds, kap, lam, mu, gstep) {
  con <- sort(con)
  A <- sort(A)
  B <- sort(B)
  scl <- sqrt(length(con) + length(A) + length(B))
  ev <- gevalR(con, A, B, gstep)
  ests <- bndestR(ev, con, A, B, kap / scl, lam / scl, mu / scl)
  T3 <- t3R(ev, con, A, B)
  ne <- length(ev)
  nx <- (ne - 1) / 2 + 1
  epos <- ev[ev >= 0]
  tiny = .Machine$double.eps^(1/2)
  v3n <- max(T3)
  w3n <- l2n1R(epos, T3)
  GnA <- ecdf(A)
  GnB <- ecdf(B)
  GA <- GnA(ev)
  GB <- GnB(ev)
  vstar <- wstar <- double(R)
  mstar1 <- mstar2 <- mstar3 <- mstar4 <- matrix(0, ne, nx)
  mmax1 <- mmax2 <- mmax3 <- mmax4 <- double(nx)
  statboot <- matrix(0, 2, R)
  for (r in 1:R) {
    sC <- sample(con, replace = TRUE)
    sA <- sample(A, replace = TRUE)
    sB <- sample(B, replace = TRUE)
    GsA <- ecdf(sA)
    GsB <- ecdf(sB)
    GstarA <- GsA(ev)
    GstarB <- GsB(ev)
    for (j in 1:nx) {
      coshift1 <- con - epos[j]
      coshift2 <- con + epos[j]
      scoshift1 <- sC - epos[j]
      scoshift2 <- sC + epos[j]
      Gc1 <- ecdf(coshift1)
      Gc2 <- ecdf(coshift2)
      Gsc1 <- ecdf(scoshift1)
      Gsc2 <- ecdf(scoshift2)
      mstar1[, j] <- GstarA - Gsc1(ev) - GA + Gc1(ev)
      mstar2[, j] <- GstarA - Gsc2(ev) - GA + Gc2(ev)
      mstar3[, j] <- Gsc1(ev) - GstarB - Gc1(ev) + GB
      mstar4[, j] <- Gsc2(ev) - GstarB - Gc2(ev) + GB
      mmax1[j] <- max(mstar1[, j] * ests$M1[, j])
      mmax2[j] <- max(mstar2[, j] * ests$M2[, j])
      mmax3[j] <- max(mstar3[, j] * ests$M3[, j])
      mmax4[j] <- max(mstar4[, j] * ests$M4[, j])
    }
    gfun <- scl * (mmax1 + mmax2 + mmax3 + mmax4)
    gfun <- pmax(0, gfun)
    vstar[r] <- max(gfun * ests$Mnec)
    wstar[r] <- l2n1R(epos, gfun * ests$csest)
  }
  rej <- c(v3n, w3n) < rbind(vstar, wstar) + tiny
  pval <- rowMeans(rej)
  ans <- list(T3 = T3, stats = c(v3n, w3n), vstar = vstar, wstar = wstar,
              rej = rej, pval = pval)
  ans
}

# For comparison
n <- 100
A <- sort(rnorm(n))
B <- sort(rnorm(n, mean = 2.8 + 8 / sqrt(3 * n)))
con <- sort(rnorm(n))
kap <- sqrt(log(log(3 * n)))
lam <- 4 * log(log((3 * n)))
mu <- kap
gstep <- 0.1
R <- 199

set.seed(1)
rtst <- LASDboot_partialR(con, A, B, R, kap, lam, mu, gstep)
cat("Result from R:", rtst$stats, "\n")

# To check that the R-only code matches the R with C++ code, uncomment these
# lines.
#library(Rcpp)
#sourceCpp("./lasd.cpp")
#set.seed(1)
#ctst <- LASDboot_partial(con, A, B, R, kap, lam, mu, gstep)
#cat("Result from R with C++:", ctst$stats)

