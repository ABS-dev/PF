library(data.table)
library(PF)
rm(list = ls())
cat("\014")

trial = data.table(animalid = 1:60,
                   group = rep(c("control", "vaccinate"), 30),
                   cage  = rep(c(1, 2, 3), each = 20),
                   sick = FALSE)

set.seed(6022)


trial[group == "control"   & cage == 1, sick := sample(c(TRUE, FALSE), 10, prob = c(.2, .8), replace = TRUE)]
trial[group == "control"   & cage == 2, sick := sample(c(TRUE, FALSE), 10, prob = c(.1, .9), replace = TRUE)]
trial[group == "control"   & cage == 3, sick := sample(c(TRUE, FALSE), 10, prob = c(.3, .7), replace = TRUE)]
trial[group == "vaccinate" & cage == 1, sick := sample(c(TRUE, FALSE), 10, prob = c(.8, .2), replace = TRUE)]
trial[group == "vaccinate" & cage == 2, sick := sample(c(TRUE, FALSE), 10, prob = c(.9, .1), replace = TRUE)]
trial[group == "vaccinate" & cage == 3, sick := sample(c(TRUE, FALSE), 10, prob = c(.7, .3), replace = TRUE)]


tbl = trial[, table(group, sick)]

M = matrix(c(12, 18, 22, 8), ncol = 2, nrow = 2, byrow = 2)
rownames(M) = c("vaccine", "control")
colnames(M) = c("sick", "well")

tbl
M


RRtosst(c(12, 30, 22, 30), rnd = 4, pf = TRUE) # This appears to be faster. Avg 1.51 & 1.68 over 100 trials
RRtosst(M, rnd = 4, pf = TRUE)                 # This apperas to be slower. Avg 1.62 & 1.72 over 100 trials

# Why isn't this working?
RRtosst(formula = diagnosis ~ group, data = trial, compare = c("vaccinate", "control"))


require(dplyr)
data1 <- data.frame(group = rep(c("treated", "control"), each = 2),
                    y = c(1, 3, 7, 5),
                    n = c(12, 12, 14, 14),
                    cage = rep(paste('cage', 1:2), 2))
data1

data2 <- data1 |>
  group_by(group) |>
  summarize(sum_y = sum(y),
            sum_n = sum(n))

# I really don't like how this formula things works.  And for some reason, I'm crashing the function...
RRtosst(data = data2, formula =  cbind(sum_y, sum_n) ~ group,
        compare = c("vaccinate", "control"))

# A way that would be more inline with how lm, glm, etc works is something like
# RRtosst(formula = sick ~ group, data = trial, compare = c("vaccinate", "control"))
#Then the function does the counting for me.
#What if we have something like:....
# RRtost(formula = sick ~ group + block(cage), data = trail, compare = c("vaccinate", "control"))
# Sick would need to be 0 / 1 or FALSE / TRUE
#
#




trial
Table6

summary = trial[,  .N, by = .(group, cage, sick)]
wide = dcast(summary, group + cage ~ sick, value.var = "N", drop = FALSE)
setnames(wide, c("FALSE", "TRUE"), c("well", "sick"))
wide[is.na(well), well := 0]
wide[is.na(sick), sick := 0]
wide[, n := well + sick]
wide

Table6
RRstr(cbind(y, n)    ~ tx    + cluster(clus), Table6, compare = c('a', 'b'))


wide
wide[, sum(sick), by = group]

RRtosst(c(8, 30, 18, 30), rnd = 4, pf = TRUE) # This appears to be faster. Avg 1.51 & 1.68 over 100 trials

RRstr(cbind(sick, n) ~ group + cluster(cage), data.frame(wide),   compare = c("vaccinate", "control")) # Crash because one value is 0
RRstr(cbind(sick, n) ~ group + cluster(cage), data.frame(wide),   compare = c("control", "vaccinate"))

# in title: "165A 12D5.30 232785 Reovirus Eff", I had to do "vaccinte" vs "control" to get it to work.  Something is strange here


wide$sick[5] = 9
wide$well[5] = 1
wide

RRstr(cbind(sick, n) ~ group + cluster(cage), data.frame(wide),   compare = c("vaccinate", "control"))
RRstr(cbind(sick, n) ~ group + cluster(cage), data.frame(wide),   compare = c("control", "vaccinate"))

