## ------------------------------------------------------------------------
require(PF)
RRsc(c(4,24,12,28))

## ------------------------------------------------------------------------
RRotsst(c(4,24,12,28))
RRtosst(c(4,24,12,28))

## ------------------------------------------------------------------------
RRstr(cbind(y,n) ~ tx + cluster(clus), Table6 , pf = F)
# Data matrix input:
# RRstr(Y = table6, pf = F)

## ------------------------------------------------------------------------
RRmh(cbind(y,n) ~ tx + cluster(clus), Table6 , pf = F)
# Data matrix input:
# RRmh(Y = table6, pf = F)

## ------------------------------------------------------------------------
bird.fit <- glm(cbind(y,n-y) ~ tx - 1, binomial, bird)
RRor(bird.fit)

## ------------------------------------------------------------------------
phiWt(bird.fit, fit.only = F)$phi
summary(update(bird.fit, family = quasibinomial))$disp

## ------------------------------------------------------------------------
# model weighted by phi hat
RRor(phiWt(bird.fit))

## ------------------------------------------------------------------------
# model with separate phi for vaccinates and controls
RRor(phiWt(bird.fit, subset.factor = bird$tx))

## ------------------------------------------------------------------------
# subtract 2 degrees of freedom
RRor(phiWt(bird.fit, subset.factor = bird$tx), degf = 2)

## ------------------------------------------------------------------------
# model weighted using tau hat
RRor(tauWt(bird.fit, subset.factor = bird$tx))

## ------------------------------------------------------------------------
# different cluster sizes, same cluster fractions
birdm.fit <- glm(cbind(y,n-y) ~ tx - 1, binomial, birdm)
RRor(tauWt(birdm.fit, subset.factor = birdm$tx))

## ------------------------------------------------------------------------
# Compare phi and tau weights
#  
phi.wts <-phiWt(birdm.fit,fit.only = F, subset.factor = birdm$tx)$weights
tau.wts <- tauWt(birdm.fit,fit.only = F, subset.factor = birdm$tx)$weights
w <- cbind(w.phi=phi.wts,w.tau=tau.wts,nw.phi=phi.wts*birdm$n,
		nw.tau=tau.wts*birdm$n)
print(cbind(birdm[,c(3,1,2)],round(w, 2)), row.names=F)

## ------------------------------------------------------------------------
# model weighted with Rao-Scott weights
RRor(rsbWt(birdm.fit, subset.factor = birdm$tx))
# just the design effect estimates
rsb(birdm$y, birdm$n, birdm$tx)$d

## ------------------------------------------------------------------------
RRlsi(c(4,24,12,28))
RRlsi(c(4,24,12,28), use.alpha = T)

## ------------------------------------------------------------------------
IDRsc(c(26,204,10,205), pf = F)

## ------------------------------------------------------------------------
IDRlsi(c(26,204,10,205), pf = F)
IDRlsi(c(26,204,10,205), pf = F, use.alpha = T)

