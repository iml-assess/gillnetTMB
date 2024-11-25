# Example 3 ####################################################################
# 4R herring data (not as perfect)
library(RTMB)

## 1) data ---------------------------------------------------------------------
load("Rdata/input.Rdata",verbose = T)

# basic info
x$mesh <- x$mesh*2.54     # make sure mesh sizes are in cm, numeric vector 
x$length <- x$length/10   # lengths are equally in cm
x$cpn <- round(x$cpn*100,0)   # catches are treated as "counts". FAILS if all values <1!!!

# pad with zeros (which are ALSO observations)
l <- data.frame(x)
dat <- dcast(l,length~mesh,value.var = "cpn",fill = 0)
dat <- melt(dat,id.vars = "length",variable.name = "mesh",value.name = "cpn")

# get some idea about how the data looks
p0 <- ggplot(dat,aes(x=length,y=cpn))+
    geom_bar(stat='identity')+
    facet_grid(mesh~.)

# reformat for model
x <- list(
    year = rep(1,nrow(dat)),
    region = rep(1,nrow(dat)),
    period = rep(1,nrow(dat)),
    length = as.numeric(as.character(dat$length)),
    mesh = as.numeric(as.character(dat$mesh)),
    cpn =as.numeric(as.character(dat$cpn))*100,
    rtype="norm.sca",
    distr = "poisson"         # poisson or nbinom
)

## 2) fit model ---------------------------------------------------------------------

# get parameters
par <- defpar(x)
par

# run model
m1 <- gillnetfitTMB(x,par) # 12 and 10 converges, but flat selectivity
m1

## 3) check output ---------------------------------------------------------------------

# tables
fittable(m1)
partable(m1)
seltable(m1)
Ntable(m1)
selmaxtable(m1)

# plots

meshconv <- data.frame(inch=c(2, 2.25, 2.5, 2.63, 2.75, 3),
                       mesh.lab=c("2 in","2 1/4 in","2 1/2 in","2 5/8 in","2 3/4 in","3 in"),
                       mesh=c(2, 2.25, 2.5, 2.63, 2.75, 3)*2.54)


plotSel(m1,meshlabs = setNames(as.character(meshconv$mesh.lab), meshconv$mesh))
plotN(m1)
plotRes(m1)
plotOP(m1)

# other check
mypred <- seltable(m1)
mypred <- merge(mypred,meshconv)
mypred$mesh.lab <- factor(mypred$mesh.lab,levels=unique(meshconv$mesh.lab))
p0 <- ggplot(mypred,aes(x=length,y=cpn))+
    geom_bar(stat='identity')+
    geom_line(aes(y=sel*max(cpn)),col="darkred")+ #selectivity rescaled to max cpn
    facet_grid(mesh.lab~.)
    
## 4) fit (scattergun) ---------------------------------------------------------------------
rtypes <- c("norm.sca","norm.loc","lognorm",'gamma')
dists <- c("poisson","nbinom")

co <- expand.grid(rtypes,dists) # 7 and 8 don't converge yet (-Inf in nll; but lognorm nbinon has great starters and bounded params??)

ms <- apply(co[1:6,],1, function(i){
    x$rtype <- i[1]
    x$distr <- i[2]
    par <- defpar(x)
    gillnetfitTMB(x,par)
})
names(ms) <- do.call(paste, c(co[1:6,], sep=" - "))

ms <- do.call('c',ms) # equivalent to c(m1,m2,m3,m4)
ms

AIC(ms)
partable(ms)

plotSel(ms,meshlabs = setNames(as.character(meshconv$mesh.lab), meshconv$mesh))
plotN(ms)
plotRes(ms)
plotOP(ms)

mypred <- seltable(ms)
mypred <- merge(mypred,meshconv)
mypred$mesh.lab <- factor(mypred$mesh.lab,levels=unique(meshconv$mesh.lab))
p0 <- ggplot(mypred[mypred$fit==names(AIC(ms))[1],],aes(x=length,y=cpn))+
    geom_bar(stat='identity')+
    geom_line(data=mypred,aes(y=sel*max(cpn),col=fit))+ #selectivity rescaled to max cpn
    facet_grid(mesh.lab~.)
