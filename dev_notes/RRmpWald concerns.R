
save(.wide, file = "C:/Users/tfkent/Documents/CVB Packages/wide.rdata")
load("C:/Users/tfkent/Documents/CVB Packages/wide.rdata")
.wide

# Observed that the NA's for group g2 in chcage 3 and 4 are actually treated as 0's for calcuating
# PF, but shouldn't this effect the confidence interval?  That doesn't seem right since there is 
# no matched pair.

t1 = dcast(.wide, chcage ~ group, value.var = "positive")[1:36]
colnames(t1) = c("chcage", "g1", "g2", "g3", "g4", "g5")
t1[, c("g4", "g5") := NULL]
t1[, g2. := g2]
t1[is.na(g2.), g2. := 0]
t1


t1[, .N, by = .(g1, g3)]
t1[, .N, by = .(g2., g3)]
t1[, .N, by = .(g2, g3)]
