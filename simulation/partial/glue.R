# Put things back together.

ind <- 1:100
#ind <- ind[-c(12)]

load(paste0("./ind_runs/run1.rda"))
rej <- 10 * simplify2array(results.list) # 10 repetitions per file

for (run in ind[-1]) {
  load(paste0("./ind_runs/run", run, ".rda"))
  rpiece <- 10 * simplify2array(results.list)
  rej <- rej + rpiece
}
rej <- rej / 1000

save(rej, file = "partialID_normal.rda")

