# Example 2 ####################################################################
## averaging npue across a factor (year/region/period) VERSUS using factor-specific model
library(RTMB)
library(ggplot2)
theme_set(theme_classic())



## Data ------------------------------------------------------------------------
# 1) dummy data from TropFish
library(TropFishR)
data(gillnet)
dimnames(gillnet$CatchPerNet_mat) <- list(gillnet$midLengths,gillnet$meshSizes)
dat0 <- reshape2::melt(gillnet$CatchPerNet_mat,varnames = c("length","mesh"),value.name = "cpn")

# 2) generate more years (copy-paste data)
ny <- 3
dat <- data.frame(sapply(dat0, rep.int,ny),year=rep(1:ny,each=nrow(dat0))) 

# 3) create differences between years (add noise) 
set.seed(123)
dat$cpn <- pmax(round(rnorm(nrow(dat),dat$cpn,75),0),0) # add some noise to make years different

p0 <- ggplot(dat,aes(x=length,y=cpn,col=as.factor(mesh)))+geom_line()+facet_wrap(~year,ncol=1)



## Option 1: take average across factor (here year) and then fit-----------------
library(plyr)
dat1 <- ddply(dat,c("length","mesh"),summarise,cpn=mean(cpn))

x <- list(
    year = rep(1,nrow(dat1)),
    region = rep(1,nrow(dat1)),
    period = rep(1,nrow(dat1)),
    length = dat1$length,
    mesh = dat1$mesh,
    cpn = dat1$cpn,
    rtype="norm.sca",
    distr = "poisson"         # poisson or nbinom
)

# parameters
par <- defpar(x)

# fit model
m1 <- gillnetfitTMB(x,par)
m1

sel <- seltable(m1)

back <- merge(dat,sel[,-1])
back$est <- with(back,cpn/sel)

p1S <- plotSel(m1)
p1N <- ggplot(back,aes(x=length,y=est))+
    geom_bar(stat="identity",position=position_dodge())+
    # labs(y=ylab,x=xlab,fill=filllab)+
    facet_wrap(year+region~period)

## Option 2: estimate in model -----------------------------------
x <- list(
    year = dat$year,
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
m2 <- gillnetfitTMB(x,par) # warnings if poor starting values
m2

p2N <- plotN(m2)
p2S <- plotSel(m2)

## Compare -----------------------------------
pars <- partable(c(m1,m2))
pars[pars$par %in% c("k1","k2"),]

library(gridExtra)
grid.arrange(p1N,p2N)
grid.arrange(p1S,p2S)
