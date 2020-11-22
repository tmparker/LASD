// This file contains functions used to to conduct inference for 
// loss averse-sensitive dominance, whether distributions are point- 
// or partially-identified.

// The C++11 library is needed for an inline function call in ord().
// [[Rcpp::plugins(cpp11)]] 
#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// ----- First part contains utilities used by other functions -----

// Construct evaluation vectors for the T1 and T2 processes.
// A and B are samples, 
// tiny is something like .Machine$double.eps^(1/2), for making sure that
// the process levels are correctly computed when the process is caglad
// instead of cadlag (just evaluate right before a change).
// [[Rcpp::export]]
NumericVector teval(NumericVector A, NumericVector B) {
  double tiny = pow(std::numeric_limits<double>::epsilon(), 0.5);
  std::vector<double> eval;
  eval.insert(eval.end(), A.begin(), A.end());
  eval.insert(eval.end(), B.begin(), B.end());
  eval.push_back(0.0);
  for (int i = 0; i < eval.size(); ++i)
    if (eval[i] < 0.0) eval[i] = -eval[i] + tiny;
  std::sort(eval.begin(), eval.end());
  std::vector<double>::iterator it;
  it = std::unique(eval.begin(), eval.end());
  eval.resize(std::distance(eval.begin(), it));
  return wrap(eval);
}

// Collect a pooled sample for the difference empirical process, used
// for FOSD process.
// A and B are samples.
// [[Rcpp::export]]
NumericVector deval(NumericVector A, NumericVector B) {
  std::vector<double> eval;
  eval.insert(eval.end(), A.begin(), A.end());
  eval.insert(eval.end(), B.begin(), B.end());
  std::sort(eval.begin(), eval.end());
  std::vector<double>::iterator it;
  it = std::unique(eval.begin(), eval.end());
  eval.resize(std::distance(eval.begin(), it));
  return wrap(eval);
}

// An ecdf() function
// eval is where to evaluate the empirical CDF 
// sobs is sorted observations.
// [[Rcpp::export]]
NumericVector edf(NumericVector eval, NumericVector sobs) {
  int neval = eval.size();
  IntegerVector ans(neval);
  for (int i = 0; i < neval; ++i)
    ans[i] = std::upper_bound(sobs.begin(), sobs.end(), eval[i]) - sobs.begin();
  return as<NumericVector>(ans) / sobs.size();
}

// Compute the difference empirical process; used for standard first-order 
// stochastic dominance tests.  The process is scaled.
// eval is where to evaluate it, A and B are samples
// [[Rcpp::export]]
NumericVector diffedf(NumericVector eval, NumericVector A, NumericVector B) {
  double nA = double(A.size());
  double nB = double(B.size());
  double scl = pow(nA + nB, 0.5);
  std::sort(A.begin(), A.end());
  std::sort(B.begin(), B.end());
  NumericVector Fa = edf(eval, A);
  NumericVector Fb = edf(eval, B);
  NumericVector diffproc = scl * (Fa - Fb);
  return diffproc;
}

// Compute the one-sided L2 norm of a one-dimensional process.
// eval is where the changes in the function occur (the function is assumed to be
// a step function) and fun is the function values at the eval points.
// [[Rcpp::export]]
double l2n1(NumericVector eval, NumericVector fun) {
  int pm1 = eval.size() - 1;
  NumericVector inside = no_init(pm1);
  for (int i = 0; i < pm1; ++i) {
    inside(i) = pow(std::max(0.0, fun(i)), 2.0);
    inside(i) *= (eval(i + 1) - eval(i));
  }
  return pow(sum(inside), 0.5);
}

// Compute the one-sided L2 norm of a two-dimensional process.
// eval is where the changes in the function occur (the function is assumed to be
// a step function) and fun is the function values at the eval points.
// [[Rcpp::export]]
double l2n2(NumericVector eval, NumericMatrix fun) {
  int pm1 = eval.size() - 1;
  NumericMatrix inside = no_init(2, pm1);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < pm1; ++j) {
      inside(i, j) = pow(std::max(0.0, fun(i, j)), 2.0);
      inside(i, j) *= (eval(j + 1) - eval(j));
    }
  }
  return pow(sum(inside), 0.5);
}

// ----- Second part is for point-identified processes -----

// Compute the T_1 process on a set of points.
// eval is where to evaluate it, A and B are samples
// [[Rcpp::export]]
NumericVector t1(NumericVector eval, NumericVector A, NumericVector B) {
  double nA = double(A.size());
  double nB = double(B.size());
  double scl = pow(nA + nB, 0.5);
  std::sort(A.begin(), A.end());
  std::sort(B.begin(), B.end());
  NumericVector na = edf(-eval, A);
  NumericVector nb = edf(-eval, B);
  NumericVector pa = edf(eval, A);
  NumericVector pb = edf(eval, B);
  NumericVector pos = pmax(0.0, pa - pb);
  NumericVector f = scl * (pos + na - nb);
  return f;
}

// Compute the T_2 process on a set of points.
// eval is where to evaluate it, A and B are samples
// [[Rcpp::export]]
NumericMatrix t2(NumericVector eval, NumericVector A, NumericVector B) {
  double nA = double(A.size());
  double nB = double(B.size());
  double scl = pow(nA + nB, 0.5);
  std::sort(A.begin(), A.end());
  std::sort(B.begin(), B.end());
  NumericVector na = edf(-eval, A);
  NumericVector nb = edf(-eval, B);
  NumericVector pa = edf(eval, A);
  NumericVector pb = edf(eval, B);
  NumericMatrix f = no_init(2, eval.size());
  f(0, _) = scl * (na - nb);
  f(1, _) = scl * (pa - pb + na - nb);
  return f;
}

// Resampling for the T_1 or T_2 processes (asymptotic behavior is same).
// A, B are the observed samples
// R is the number of bootstrap repetitions
// lam is the contact set estimator cutoff for the process when it
// is evaluated with sup-norm or L2 functionals.
// Returns a 2xR matrix of sup-T2 and L2-T2 statistics.

// [[Rcpp::export]]
List LASDboot(NumericVector A, NumericVector B, int R, 
              double kap, double lam, double mu) {
  double tiny = pow(std::numeric_limits<double>::epsilon(), 0.5);
  int nA = A.size();
  int nB = B.size();
  NumericVector eval = teval(A, B);
  int ne = eval.size();
  NumericVector T1 = t1(eval, A, B); // sorts the observations
  NumericMatrix T2 = t2(eval, A, B);
  NumericVector stats(4);
  stats[0] = std::max(0.0, double(max(T1)));
  stats[1] = std::max(0.0, double(max(T2)));
  stats[2] = l2n1(eval, T1);
  stats[3] = l2n2(eval, T2);
  // Setup for bootstrap - estimated pieces.
  double maxm1 = double(max(T2(0, _)));
  double maxm2 = double(max(T2(1, _)));
  int mind = 0;
  if (maxm1 > maxm2 + mu) mind = 1;
  if (maxm1 < maxm2 - mu) mind = 2;
  NumericVector epsmax1 = ifelse(T2(0, _) >= maxm1 - kap, 1.0, 0.0);
  NumericVector epsmax2 = ifelse(T2(1, _) >= maxm2 - kap, 1.0, 0.0);
  NumericVector con1 = ifelse(abs(T2(0, _)) < lam, 1.0, 0.0);
  NumericVector con2 = ifelse(abs(T2(1, _)) < lam, 1.0, 0.0);
  if (sum(con1) == 0.0) std::fill(con1.begin(), con1.end(), 1.0);
  if (sum(con2) == 0.0) std::fill(con2.begin(), con2.end(), 1.0);
  // Bootstrap repetitions
  double v1b;
  double v2b;
  double vb;
  double wb;
  NumericVector sA(nA);
  NumericVector sB(nB);
  NumericMatrix f2b(2, ne);
  NumericMatrix pcount(4, R);
  NumericVector pval(4);
  for (int i = 0; i < R; ++i) {
    sA = sample(A, nA, true);
    sB = sample(B, nB, true);
    f2b = t2(eval, sA, sB); // sorts the bootstrap observations
    for (int j = 0; j < 2 * ne; ++j)
      f2b[j] = f2b[j] - T2[j];
    v1b = std::max(0.0, double(max(f2b(0, _) * epsmax1)));
    v2b = std::max(0.0, double(max(f2b(1, _) * epsmax2)));
    if (mind == 0) vb = std::max(v1b, v2b);
    if (mind == 1) vb = v1b;
    if (mind == 2) vb = v2b;
    for (int j = 0; j < ne; ++j) {
      f2b(0, j) = f2b(0, j) * con1(j);
      f2b(1, j) = f2b(1, j) * con2(j);
    }
    wb = l2n2(eval, f2b);
    // Calculate p-values
    if (stats[0] < vb + tiny) pcount(0, i) = 1.0;
    if (stats[1] < vb + tiny) pcount(1, i) = 1.0;
    if (stats[2] < wb + tiny) pcount(2, i) = 1.0;
    if (stats[3] < wb + tiny) pcount(3, i) = 1.0;
  }
  for (int j = 0; j < 4; ++j)
    pval[j] = mean(pcount(j, _));
  return List::create(_["stats"] = stats, _["pval"] = pval);
}

// Resampling for the FOSD process.
// A, B are the observed samples
// R is the number of bootstrap repetitions
// lam is the contact set estimator cutoff for the process when it
// is evaluated with sup-norm or L2 functionals.
// Returns a 2xR matrix of sup-norm and L2-norm statistics.

// [[Rcpp::export]]
List FOSDboot(NumericVector A, NumericVector B, int R, double lam) {
  double tiny = pow(std::numeric_limits<double>::epsilon(), 0.5);
  int nA = A.size();
  int nB = B.size();
  std::sort(A.begin(), A.end());
  std::sort(B.begin(), B.end());
  NumericVector eval = deval(A, B); // for FOSD process
  int ne = eval.size();
  NumericVector dF = diffedf(eval, A, B);
  NumericVector conF = ifelse(abs(dF) < lam, 1.0, 0.0);
  if (is_true(all(conF == 0.0)))
    std::fill(conF.begin(), conF.end(), 1.0);
  for (int i = 0; i < ne; ++i)
    dF[i] *= conF[i];
  double vn = std::max(0.0, double(max(dF)));
  double wn = l2n1(eval, dF);
  double vb;
  double wb;
  // Bootstrap repetitions
  NumericVector sA(nA);
  NumericVector sB(nB);
  NumericVector bootF(ne);
  NumericMatrix rej(2, R);
  NumericVector pval(2);
  for (int i = 0; i < R; ++i) {
    sA = sample(A, nA, true);
    sB = sample(B, nB, true);
    bootF = diffedf(eval, sA, sB);
    for (int j = 0; j < ne; ++j)
      bootF[j] = (bootF[j] - dF[j]) * conF[j];
    vb = std::max(0.0, double(max(bootF)));
    wb = l2n1(eval, bootF);
    if (vn < vb + tiny) rej(0, i) = 1.0;
    if (wn < wb + tiny) rej(1, i) = 1.0;
  }
  pval[0] = mean(rej.row(0));
  pval[1] = mean(rej.row(1));
  NumericVector stats(2);
  stats[0] = vn;
  stats[1] = wn;
  return List::create(_["stats"] = stats, _["pval"] = pval);
}

// A bootstrap test for the equality of two distribution functions.
// [[Rcpp::export]]
List eqboot(NumericVector A, NumericVector B, int R) {
  int nA = A.size();
  int nB = B.size();
  std::sort(A.begin(), A.end());
  std::sort(B.begin(), B.end());
  NumericVector eval = deval(A, B);
  int ne = eval.size();
  NumericVector dF = diffedf(eval, A, B);
  double vn = std::max(0.0, double(max(abs(dF))));
  double wn = l2n1(eval, abs(dF));
  double vb;
  double wb;
  // Bootstrap repetitions
  NumericVector sA(nA);
  NumericVector sB(nB);
  NumericVector bootF(ne);
  NumericMatrix rej(2, R);
  NumericVector pval(2);
  for (int i = 0; i < R; ++i) {
    sA = sample(A, nA, true);
    sB = sample(B, nB, true);
    bootF = diffedf(eval, sA, sB);
    for (int j = 0; j < ne; ++j)
      bootF[j] = bootF[j] - dF[j];
    vb = std::max(0.0, double(max(abs(bootF))));
    wb = l2n1(eval, abs(bootF));
    if (vn < vb) rej(0, i) = 1.0;
    if (wn < wb) rej(1, i) = 1.0;
  }
  pval[0] = mean(rej.row(0));
  pval[1] = mean(rej.row(1));
  NumericVector stats(2);
  stats[0] = vn;
  stats[1] = wn;
  return List::create(_["stats"] = stats, _["pval"] = pval);
}

// ----- This part is for interval-identified inference processes -----

// Compute lower and upper bounds for the CDF of the difference 
// between treatment and control using samples tr and co.
// Avoids use of diffedf because of the repeated sorting and evaluation of the 
// treatment ecdf (which only needs one evaluation) that involves.
// First row of bounds is lower bound, second row is upper bound.
// [[Rcpp::export]]
NumericMatrix lohi(NumericVector eval, NumericVector tr, NumericVector co) {
  int ne = eval.size();
  std::sort(tr.begin(), tr.end());
  std::sort(co.begin(), co.end());
  NumericVector Gtr = edf(eval, tr);
  NumericMatrix bounds = no_init(2, ne);
  NumericVector coshift(co.size());
  NumericVector Gco(ne);
  NumericVector dpro(ne);
  for (int i = 0; i < ne; ++i) {
    coshift = co + eval[i];
    Gco = edf(eval, coshift);
    dpro = Gtr - Gco;
    dpro.push_back(0.0);
    bounds(0, i) = max(dpro);
    bounds(1, i) = 1.0 + min(dpro);
  }
  return bounds;
}

// Construct an evaluation grid for the T3 process
// A B and C are samples (for treatments A&B and C=Control)
// gstep is the grid diameter.
// [[Rcpp::export]]
NumericVector geval(NumericVector C, NumericVector A, NumericVector B, double gstep) {
  NumericVector lims = NumericVector::create(min(A) - max(C),
    max(A) - min(C), min(B) - max(C), max(B) - min(C));
  int glen = ceil(max(abs(lims)) / gstep);
  double gtop = double(glen) * gstep;
  int ulen = 2 * glen + 1;
  NumericVector grd(ulen);
  for (int i = 0; i < ulen; ++i)
    grd[i] = -gtop + double(i) * gstep;
  return grd;
}

// Compute the T_3 process
// eval is where the functions are evaluated - in later functions, assumed that
// geval() is used to construct eval but not here.
// C, A and B correspond to X_0, X_A and X_B realizations.
// [[Rcpp::export]]
NumericVector t3(NumericVector eval, NumericVector C, NumericVector A, 
                  NumericVector B) {
  double nA = double(A.size());
  double nB = double(B.size());
  double n0 = double(C.size());
  double scl = pow(nA + nB + n0, 0.5);
  NumericVector epos = eval[eval >= 0];
  NumericVector nLA = lohi(-epos, A, C)(0, _);
  NumericVector nUB = lohi(-epos, B, C)(1, _);
  NumericVector pLA = lohi(epos, A, C)(0, _);
  NumericVector pUB = lohi(epos, B, C)(1, _);
  NumericVector f = scl * (pLA - pUB + nLA - nUB);
  return f;
}

// Make all the contact set & epsilon-maximizer estimates that are needed to
// bootstrap the sup-T3 or L2-T3 statistics.
// Assumed that eval is a geval() object.
// kap = epsilon-maximizer constant for the four subfunctions
// lam = contact set constant
// mu = epsilon-maximizer constant for the outer optimization problem.

// [[Rcpp::export]]
List estT3(NumericVector eval, NumericVector C, NumericVector A, 
          NumericVector B, double kap, double lam, double mu) {
  int ne = eval.size();
  int nx = (ne - 1) / 2 + 1;
  NumericVector epos = eval[eval >= 0];
  std::sort(A.begin(), A.end());
  std::sort(B.begin(), B.end());
  std::sort(C.begin(), C.end());
  NumericVector coshift1(C.size());
  NumericVector coshift2(C.size());
  NumericVector GA = edf(eval, A);
  NumericVector GB = edf(eval, B);
  NumericMatrix m1(ne, nx);
  NumericMatrix m2(ne, nx);
  NumericMatrix m3(ne, nx);
  NumericMatrix m4(ne, nx);
  NumericVector mm1(nx);
  NumericVector mm2(nx);
  NumericVector mm3(nx);
  NumericVector mm4(nx);
  LogicalMatrix M1(ne, nx);
  LogicalMatrix M2(ne, nx);
  LogicalMatrix M3(ne, nx);
  LogicalMatrix M4(ne, nx);
  LogicalVector csest(nx);
  LogicalVector Mnec(nx);
  for (int j = 0; j < nx; ++j) {
    coshift1 = C - epos[j];
    coshift2 = C + epos[j];
    m1(_, j) = GA - edf(eval, coshift1);
    m2(_, j) = GA - edf(eval, coshift2);
    m3(_, j) = edf(eval, coshift1) - GB;
    m4(_, j) = edf(eval, coshift2) - GB;
    mm1[j] = max(m1(_, j));
    mm2[j] = max(m2(_, j));
    mm3[j] = max(m3(_, j));
    mm4[j] = max(m4(_, j));
    for (int i = 0; i < ne; ++i) {
      if (m1(i, j) >= mm1[j] - kap) M1(i, j) = true;
      if (m2(i, j) >= mm2[j] - kap) M2(i, j) = true;
      if (m3(i, j) >= mm3[j] - kap) M3(i, j) = true;
      if (m4(i, j) >= mm4[j] - kap) M4(i, j) = true;
    }
  }
  NumericVector t3pro = mm1 + mm2 + mm3 + mm4 - 2.0;
  for (int j = 0; j < nx; ++j) {
    if (std::abs(t3pro[j]) < lam) csest[j] = true;
    if (t3pro[j] >= max(t3pro) - mu) Mnec[j] = true;
  }
  if (sum(csest) == 0) csest = rep(true, nx);
  return List::create(_["M1"] = M1, _["M2"] = M2, 
      _["M3"] = M3, _["M4"] = M4, //_["t3pro"] = t3pro, 
      _["csest"] = csest, _["Mnec"] = Mnec);
}

// Bootstrap test for sup-T3 and L2-T3 statistics.
// Assumes that eval is a geval() object.
// [[Rcpp::export]]
List LASDboot_partial(NumericVector C, NumericVector A, NumericVector B, 
                      int R, double kap, double lam, double mu, 
                      double gstep) {
  double nA = double(A.size());
  double nB = double(B.size());
  double n0 = double(C.size());
  double scl = pow(nA + nB + n0, 0.5);
  std::sort(A.begin(), A.end());
  std::sort(B.begin(), B.end());
  std::sort(C.begin(), C.end());
  NumericVector eval = geval(C, A, B, gstep);
  // parameters scaled in the next part to match the other routines' args.
  List ests = estT3(eval, C, A, B, kap / scl, lam / scl, mu / scl);
  NumericVector T3 = t3(eval, C, A, B);
  double tiny = pow(std::numeric_limits<double>::epsilon(), 0.5);
  int ne = eval.size();
  int nx = (ne - 1) / 2 + 1;
  NumericVector epos(nx);
  std::copy(eval.begin() + nx - 1, eval.end(), epos.begin());
  double v3n = std::max(0.0, double(max(T3)));
  double w3n = l2n1(epos, T3);
  NumericVector sC(C.size());
  NumericVector sA(A.size());
  NumericVector sB(B.size());
  NumericVector coshift1(C.size());
  NumericVector coshift2(C.size());
  NumericVector scoshift1(C.size());
  NumericVector scoshift2(C.size());
  NumericVector GA = edf(eval, A);
  NumericVector GB = edf(eval, B);
  NumericVector GstarA(ne);
  NumericVector GstarB(ne);
  NumericVector Gc1(ne);
  NumericVector Gc2(ne);
  NumericVector Gsc1(ne);
  NumericVector Gsc2(ne);
  NumericMatrix mstar1(ne, nx);
  NumericMatrix mstar2(ne, nx);
  NumericMatrix mstar3(ne, nx);
  NumericMatrix mstar4(ne, nx);
  NumericVector mmax1(nx);
  NumericVector mmax2(nx);
  NumericVector mmax3(nx);
  NumericVector mmax4(nx);
  NumericVector gfun(nx);
  NumericVector vstar(R);
  NumericVector wstar(R);
  NumericMatrix rej(2, R);
  NumericMatrix M1 = ests["M1"];
  NumericMatrix M2 = ests["M2"];
  NumericMatrix M3 = ests["M3"];
  NumericMatrix M4 = ests["M4"];
  NumericVector Mnec = ests["Mnec"];
  NumericVector csest = ests["csest"];
  for (int r = 0; r < R; ++r) {
    sC = sample(C, n0, true);
    sA = sample(A, nA, true);
    sB = sample(B, nB, true);
    std::sort(sA.begin(), sA.end());
    std::sort(sB.begin(), sB.end());
    std::sort(sC.begin(), sC.end());
    GstarA = edf(eval, sA);
    GstarB = edf(eval, sB);
    for (int j = 0; j < nx; ++j) {
      coshift1 = C - epos[j];
      scoshift1 = sC - epos[j];
      coshift2 = C + epos[j];
      scoshift2 = sC + epos[j];
      Gc1 = edf(eval, coshift1);
      Gc2 = edf(eval, coshift2);
      Gsc1 = edf(eval, scoshift1);
      Gsc2 = edf(eval, scoshift2);
      mstar1(_, j) = GstarA - Gsc1 - GA + Gc1;
      mstar2(_, j) = GstarA - Gsc2 - GA + Gc2;
      mstar3(_, j) = Gsc1 - GstarB - Gc1 + GB;
      mstar4(_, j) = Gsc2 - GstarB - Gc2 + GB;
      mmax1[j] = max(mstar1(_, j) * M1(_, j));
      mmax2[j] = max(mstar2(_, j) * M2(_, j));
      mmax3[j] = max(mstar3(_, j) * M3(_, j));
      mmax4[j] = max(mstar4(_, j) * M4(_, j));
    }
    gfun = scl * (mmax1 + mmax2 + mmax3 + mmax4);
    gfun = pmax(0.0, gfun);
    vstar[r] = max(gfun * Mnec);
    wstar[r] = l2n1(epos, gfun * csest);
    if (v3n < vstar[r] + tiny) rej(0, r) = 1.0;
    if (w3n < wstar[r] + tiny) rej(1, r) = 1.0;
  }
  NumericVector pval = NumericVector::create(mean(rej.row(0)), 
                        mean(rej.row(1)));
  NumericVector stats = NumericVector::create(v3n, w3n);
  return List::create(_["stats"] = stats, _["pval"] = pval);
//  return List::create(_["T3"] = T3, _["stats"] = stats, _["vstar"] = vstar, 
//                      _["wstar"] = wstar, _["rej"] = rej, _["pval"] = pval);
}

