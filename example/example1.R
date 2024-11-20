# Example 1 ####################################################################
## Compare Tropfish with gillnetTMB
library(RTMB)

## Data --------------------------------------------------------------------
library(TropFishR)
data(gillnet)

rtypes <- c("norm.sca","norm.loc","lognorm",'gamma')
dists <- c("poisson","nbinom")

## Tropfish --------------------------------------------------------------------
dat0 <- matrix(c(gillnet$midLengths, gillnet$CatchPerNet_mat),byrow = FALSE, ncol=(dim(gillnet$CatchPerNet_mat)[2]+1))

tropfits <- lapply(rtypes, function(x){
    g <- gillnetfit(data = dat0, meshsizes = gillnet$meshSizes,rtype=x,details=TRUE) 
    g$gear.pars
})

names(tropfits) <- rtypes

## GillnetTMB --------------------------------------------------------------------
dimnames(gillnet$CatchPerNet_mat) <- list(gillnet$midLengths,gillnet$meshSizes)
dat <- reshape2::melt(gillnet$CatchPerNet_mat,varnames = c("length","mesh"),value.name = "cpn")

p0 <- ggplot(dat,aes(x=length,y=cpn,col=as.factor(mesh)))+geom_line()

# data format (for every different year, region or period, will estimate different N with same selecti)
x <- list(
    year = rep(1,nrow(dat)),
    region = rep(1,nrow(dat)),
    period = rep(1,nrow(dat)),
    length = dat$length,
    mesh = dat$mesh,
    cpn = dat$cpn,
    rtype="norm.sca",
    distr = "poisson"         # poisson or nbinom
)

# parameters
par <- defpar(x)

# fit model
m1 <- gillnetfitTMB(x,par) # warnings if poor starting values
m1

partable(m1)[1:2,2:3] # estimated paras
tropfits$norm.sca[1:2,] # how do they compare to tropfits?

# Example 1b ####################################################################
# Demo of other options

## Other selectivity curves ######################
x$rtype <- "norm.loc"
par <- defpar(x)
m2 <- gillnetfitTMB(x,par)
m2

partable(m2)[1:2,2:3]
tropfits$norm.loc[1:2,]

x$rtype <- "lognorm"
par <- defpar(x)
m3 <- gillnetfitTMB(x,par)
m3

partable(m3)[1:2,2:3]
tropfits$lognorm[1:2,]

x$rtype <- "gamma"
par <- defpar(x)
m4 <- gillnetfitTMB(x,par)
m4

partable(m4)[1:2,2:3]
tropfits$gamma[1:2,]

## Negative binomial  ######################
x$rtype <- "norm.sca"
x$distr <- "nbinom"
par <- defpar(x)
m5 <- gillnetfitTMB(x,par) # warning is ok if all paras have sd, and max gradient <0.001
m5

partable(m5)[1:2,2:3]

# Example 1c ####################################################################
# Demo of functions

## Plots and tables ------------------------------------------------------

### 1 model
fittable(m1)
partable(m1)
seltable(m1)
Ntable(m1)
selmaxtable(m1)

plotSel(m1)
plotN(m1)
plotRes(m1)
plotOP(m1)

### multiple models
x$distr <- "poisson"
ms <- lapply(rtypes, function(i){
    x$rtype <- i
    par <- defpar(x)
    gillnetfitTMB(x,par)
})
names(ms) <- rtypes
ms <- do.call('c',ms) # equivalent to c(m1,m2,m3,m4)
ms

fittable(ms)
AIC(ms)
partable(ms)
Ntable(ms)
seltable(ms)
selmaxtable(ms)
plotRes(ms)

plotSel(ms)
plotN(ms)
plotRes(ms)
plotOP(ms)


