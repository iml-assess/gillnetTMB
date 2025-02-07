# Example 1: Compare Tropfish with gillnetTMB ##################################
#library(gillnetTMB)
library(TropFishR)
library(RTMB)

## A) Baseline ######################

data(gillnet)
rtypes <- c("norm.sca","norm.loc","lognorm",'gamma')

### Tropfish --------------------------------------------------------------------
dat0 <- matrix(c(gillnet$midLengths, gillnet$CatchPerNet_mat),byrow = FALSE, ncol=(dim(gillnet$CatchPerNet_mat)[2]+1))

tropfits <- lapply(rtypes, function(x){
    g <- gillnetfit(data = dat0, meshsizes = gillnet$meshSizes,rtype=x,details=TRUE) 
    g$gear.pars
})

names(tropfits) <- rtypes

### GillnetTMB --------------------------------------------------------------------
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

### Comparison -----------------------------------------------------------------

partable(m1)[1:2,2:3] # estimated paras
tropfits$norm.sca[1:2,] # how do they compare to tropfits?

if(!identical(round(unlist(partable(m1)[1:2,2:3],use.names=F),4),
              round(as.vector(tropfits$norm.sca[1:2,]),4))) warning("model output not identical")

## B) Selectivity curves #######################################################

# fit all selectivity curves

ms <- lapply(rtypes, function(i){
    x$rtype <- i
    par <- defpar(x)
    gillnetfitTMB(x,par)
})
names(ms) <- rtypes
ms <- do.call('c',ms) # equivalent to c(m1,m2,m3,m4)
ms

# % error in parameters

aa <- round(sapply(ms,function(x)unlist(partable(x)[1:2,2:3])),4)
bb <- round(sapply(tropfits,function(x)as.vector(x[1:2,])),4)
(aa-bb)*100/aa

## C) Distributions ####################################################################
# Tropfish only has poisson distribution, gillnetTMB gives the choice with the negative binomial
x$rtype <- "norm.sca"
x$distr <- "nbinom"
par <- defpar(x)
m2 <- gillnetfitTMB(x,par) # warning is ok if all paras have sd, and max gradient <0.001
m2

partable(m2)[1:2,2:3]

## D) DEMO of functions ########################################################

### 1 model --------------------------------------------------------------------
fittable(m1)
partable(m1)
seltable(m1)
Ntable(m1)
selmaxtable(m1)

plotSel(m1)
plotN(m1)
plotRes(m1)
plotOP(m1)

### multiple models  -----------------------------------------------------------
fittable(ms)
AIC(ms)
partable(ms)
Ntable(ms)
seltable(ms)
selmaxtable(ms)

plotSel(ms)
plotN(ms)
plotRes(ms)
plotOP(ms)



