# Example 1: Compare Tropfish with gillnetTMB ##################################
library(gillnetTMB)
library(ggplot2)

## A) Baseline ######################

data(gillnet)
rtypes <- c("norm.sca","norm.loc","lognorm",'gamma')

### Tropfish --------------------------------------------------------------------
dat0 <- matrix(c(gillnet$midLengths, gillnet$CatchPerNet_mat),byrow = FALSE, ncol=(dim(gillnet$CatchPerNet_mat)[2]+1))
add <- cbind(seq(74.5,80.5,by=2),matrix(0,nrow=4,ncol=8))
dat0 <- rbind(dat0,add)

tropfits <- lapply(rtypes, function(x){
    g <- TropFishR::gillnetfit(data = dat0, meshsizes = gillnet$meshSizes,rtype=x,details=TRUE) 
    g$gear.pars
})

names(tropfits) <- rtypes

### GillnetTMB --------------------------------------------------------------------
d <- dat0[,-1]
dimnames(d) <- list("length"= dat0[,1],"mesh" = gillnet$meshSizes)
dat <- reshape2::melt(d,value.name = "cpn")

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
par <- gillnetTMB::defpar(x)

# fit model
m1 <- gillnetTMB::gillnetfit(x,par)
m1 


### Comparison -----------------------------------------------------------------
partable(m1)[1:2,2:3]     # estimated paras
tropfits[[x$rtype]][1:2,] # how do they compare to tropfits?

## B) Selectivity curves #######################################################

# fit all selectivity curves

ms <- lapply(rtypes, function(i){
    x$rtype <- i
    par <- defpar(x)
    gillnetfit(x,par)
})
names(ms) <- rtypes
ms <- do.call('c',ms)
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
m2 <- gillnetfit(x,par) # warning is ok if all paras have sd, and max gradient <0.001
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
plotOP(m1,facet = T)

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

## E) Predictions ########################################################

new <- predict(ms,length=seq(40,90,0.01))
ggplot(new,aes(x=length,y=y,col=fit))+
    geom_line()+
    facet_wrap(~mesh)


