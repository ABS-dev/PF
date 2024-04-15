################################################################################
############ Utilities to check the code for style and other issues ############
################################################################################

library(data.table)
library(stringr)

# Lint the package to check for formatting errors
# Path should be the root of the package directory.
path <- getwd()
linter <- N <- NULL

# Variables that don't follow stylistic guidelines, which we are not changing.
okay_vars <- c("A", "B", "C", "L", "N", "Q", "R", "V", "Y", "L1", "L2",
               "Lk", "MN", "Q1", "Q2", "U1", "U2", "Xx", "Y1", "Y2", "Phi",
               "Tau", "Y$C", "d.i", "m.i", "n.i", "p.i", "u.p", "v.i", "y.i",
               "z.s", "RRmh", "RRor", "RRsc", "ci.d", "g.ph", "logR", "r.ij",
               "varR", "wTau", "z.ph", "IDRsc", "RRlsi", "RRstr", "Terms",
               "V.eta", "ci.dl", "gradR", "p.set", "phiWt", "q.set", "r.hat",
               "r.max", "r.min", "rsbWt", "tauWt", "var.b", "y.dot", "z.ah2",
               "z.al2", "z.new", "z.old", "z.phi", "IDRlsi", "Phi.ML",
               "V.beta", "ll.max", "mu.hat", "p.test", "rr.est", "rr.opt",
               "scst.y", "zi.phi", "RRotsst", "RRtosst", "df.test", "idr.hat",
               "num.dig", "old.low", "phi.mle", "phi.new", "phi.old",
               "rho.mle", "tau.hat", "tau.new", "varlogR", "zis.phi", "zsc.phi",
               "zsk.phi", "RRmpWald", "beta.hat", "bird.fit", "ci.asymp",
               "fit.only", "gradlogR", "hom.test", "iter.max", "old.high",
               "rr.ci.hi", "rr.ci.lo", "tauOptim", "trace.it", "Table6$tx",
               "birdm.fit", "this.call", "use.alpha", "var.log.r", "full.track",
               "grad.log.r", "p.set$n1y1", "p.set$n2y2", "q.set$n1y1",
               "q.set$n2y2", "show.warns", "var.log.rr", "drop.levels",
               "score.start", "options.warn", "var.beta.hat", "show.warnings",
               "newfamily.name", ".rr.score.asymp", "nuisance.points")

okay_vars <- c()

## Change This Value ##
# Determine whose files should be linted.  Change the number in the [.]

excluded_files <- included_files <- NULL

idx <- 1
excluded_linters <- c("indentation_linter",
                      "cyclocomp_linter",
                      "commented_code_linter",
                      "object_usage_linter"
                      )[idx]

# Run lintr on the selected files.
tmp <- codeDiagnostics::lint_package_extended(
  path = path,
  okay_vars = okay_vars,
  excluded_files = excluded_files,
  included_files = included_files,
  excluded_linters = excluded_linters)

# Send results to Markers screen
if (length(tmp) == 0) {
  message("Sucess!  No lintr issues!")
} else {
  tmp
}

## Run to this point to see lint results.
# CTRL-SHIFT-B

# Format results as a data.table for further exploration
dt <- data.table::as.data.table(tmp)
dt[, .N, linter][order(N)]
dt[, .N, filename]
dt[, .N, keyby = .(filename, linter)]


tmp[dt[, which(grepl("name", linter))]]



lf <- paste0("(", paste0(lindsay_files, collapse = "|"), ")")
dt[!str_detect(filename, lf), .N, filename]
tmp[!dt[, str_detect(filename, lf)]]


?assign


tmp <- dt[str_detect(linter, "name"),
          str_extract(trimws(line), "^[^ ]*"), line][, .N, keyby = V1]
cat(tmp[, paste0("\"", V1, "\"", collapse = ", ")])
