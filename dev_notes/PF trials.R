library(data.table)
library(PF)

n = 20

dt = data.table(pair = rep(1:n, each = 2),
                trt  = rep(c("v", "c"), n),
                diagnosis = sample(c(0, 1), 2 * n, replace = TRUE))

temp = RRmpWald(diagnosis ~ trt + cluster(pair), dt, compare = c("v", "c"), rnd = 4)


estimate = data.table(PF = numeric(1330),
                      LL = numeric(1330),
                      UL = numeric(1330),
                      PPi = numeric(1330),
                      NPi = numeric(1330),
                      PNi = numeric(1330))
k = 1
for (PP in 0:n) {
  nPP = n - PP
  if (nPP >= 0) {
    for (NP in 0:nPP) {
      nPP_NP = nPP - NP
      if (nPP_NP >= 0) {
        for (PN in 0:nPP_NP) {
          NN = nPP_NP - PN
          if (NP > 0 & PN > 0) { # There are strange errors if both or one is 0
            dt[, diagnosis := c(rep(1, 2 * PP),
                                rep(c(1, 0), PN),
                                rep(c(0, 1), NP),
                                rep(c(0), 2 * NN))]
            temp = RRmpWald(diagnosis ~ trt + cluster(pair), dt, compare = c("v", "c"), rnd = 4)$estimate
            estimate[k, PF := temp[1]]
            estimate[k, LL := temp[2]]
            estimate[k, UL := temp[3]]
            estimate[k, PPi := PP]
            estimate[k, PNi := PN]
            estimate[k, NPi := NP]
            k = k + 1
          }
        }
      }
    }
  }
}
estimate

estimate[, vacPos := PPi + PNi]
estimate[, conPos := PPi + NPi]
setorder(estimate, vacPos, conPos)
estimate[1:30]
pf = estimate[, .(PF = mean(PF),
                  min_LL = min(LL),
                  max_LL = max(LL),
                  min_UL = min(UL),
                  max_UL = max(UL)),
               by = .(vacPos, conPos)]
pf
pf[PF > 0 & min_LL != max_LL]
estimate[PF == 0]

pf[PF > 0, max_LL - min_LL]
pf[PF > 0, max_UL - min_UL]
