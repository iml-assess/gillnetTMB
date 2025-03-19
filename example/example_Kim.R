#############################################################################
### The way I would analyse the data based on the scripts produced by EVB ###
#############################################################################

# note: some packages already loaded in .Rprofile

### Packages
library(RTMB)
library(dplyr)
library(ggplot2)
library(FSA)
library(lubridate)
library(reshape2)
library(tidyr)

### Functions (not necessary: if .Rprofile is in Rproject these load automatically when opening project)
source("./EVB/gillnetTMB-main/R/fit.R")
source("./EVB/gillnetTMB-main/R/methods.R")
source("./EVB/gillnetTMB-main/R/plots.R")
source("./EVB/gillnetTMB-main/R/predict.R")
source("./EVB/gillnetTMB-main/R/tables.R")

### Import dataset
gillnet <- read.csv("data/gillnet_data_2021_2023_kim.csv",sep=",")

### Data preparation

# Length class
# eli: gillnet$LCat <- 0.5*round(gillnet$longt/10/0.5) 
gillnet$long.cm <- gillnet$longt/10
minl <- floor(min(gillnet$long.cm, na.rm=T)) # replaced with round with floor for future (minl should be lower than observed)
gillnet <- lencat(~long.cm, data=gillnet, startcat=minl, w=0.5) # I don't get why 27.336 should be 27 and not 27.5 (row 1). L29 is an easier oneliner?


# Panel surface area
gillnet$surf.net= ifelse(gillnet$an == "2021", 5.4, 37.9)

# Add region
# I consider only two regions (4Sw and 4R)
gillnet$region <- ifelse(gillnet$peche=="4SW", "4SW", "4R")

# Add period
gillnet$period <- ifelse(month(gillnet$datetime_set)<=7, "spring", "summer")

# Mesh size in cm
gillnet$mesh.cm <- gillnet$maille*2.54 
# Remove S_ECH 71810 and 71811 (no mesh sizes) and keep only spring and fall spawners (remove unknown)
gillnet <- gillnet %>% filter(!S_ECH %in% c(71810,71811)) %>% filter(choixc %in% c("A","P"))

# Improve estimates of catch weight (in cases where the entire catch was sampled)
gillnet <- gillnet %>%
  group_by(S_ECH,maille) %>%
  mutate(catch2 = case_when(an==2021 & no_db<30 ~ sum(poidt)/1000,
                            an>2021 & no_db<50 ~ sum(poidt)/1000,
                            .default=catch))

# Remove observations with no length
cpn <- gillnet[!is.na(gillnet$longt),]

# Standardize length frequency distributions for soak time, net size and sample size
cpn <- cpn %>%
  group_by(S_ECH, an, choixc, LCat, region, period, mesh.cm) %>%
  mutate(cpn=n(), cpue=cpn*catch2/soak_time*surf.net*(poidt/1000))

# Average over replicates (and period) for each region and pad with zeros
cpn.avg <- cpn %>%
  group_by(an, choixc, region, mesh.cm, LCat) %>%
  summarize(mean.cpue=mean(cpue)) %>%
  ungroup() %>%
  complete(an, choixc, region, mesh.cm, LCat, fill=list(mean.cpue=0))

# What does the data look like?
ggplot(cpn.avg,aes(x=LCat,y=mean.cpue))+
  geom_bar(stat="identity")+
  facet_grid(mesh.cm~region+choixc,scale="free")

# Catches are treated as "counts"; FAILS if all values < 1
cpn.avg$mean.cpue <- round(cpn.avg$mean.cpue*1000,0) 

### Fit model for fall spawners in 4Sw (like Mélanie)
dat.fall <- cpn.avg %>%
  filter(choixc=="A" & region=="4SW")

# What does the data look like? ELI: removed scale free and added year 
ggplot(dat.fall,aes(x=LCat,y=mean.cpue))+
  geom_bar(stat='identity')+
  facet_grid(mesh.cm~an)

# eli: temproary work without poor year (2021) and mesh (6.68)
dat.fall <- dat.fall[dat.fall$mesh.cm!=6.6802&dat.fall$an!=2021,]

# ELI: fill gaps
ggplot(dat.fall,aes(x=LCat,y=mesh.cm))+ # eli: bubbles to easier spot gaps -> still some left
    geom_point(aes(size=mean.cpue))+
    facet_grid(an~.)

lookup <- expand.grid(an=unique(dat.fall$an),
                      region=unique(dat.fall$region),
                      mesh.cm=unique(dat.fall$mesh.cm),
                      LCat=seq(min(dat.fall$LCat),max(dat.fall$LCat),by=0.5))
dat.fall <- merge(lookup,dat.fall,all=T)
dat.fall[is.na(dat.fall$mean.cpue),"mean.cpue"] <- 0

# Reformat for model
x <- list(
  year = dat.fall$an,
  #year = rep(1,nrow(dat.fall)),
  #region = dat.spr$region,
  region = rep(1,nrow(dat.fall)),
  #period = dat.fall$period,
  period = rep(1,nrow(dat.fall)),
  length = dat.fall$LCat,
  mesh = dat.fall$mesh.cm,
  cpn = dat.fall$mean.cpue,
  rtype="norm.sca",
  distr = "poisson"         # poisson or nbinom
)

# Get parameters
(par <- defpar(x))

# Run model 
(m1 <- gillnetfitTMB(x,par)) # this should say "convergence" at the end

# no convergence -> can I get better par estimates?
dat.falla <- ddply(dat.fall,c("LCat","mesh.cm"),summarise,mean.cpue=mean(mean.cpue))
x <- list(
    #year = dat.fall$an,
    year = rep(1,nrow(dat.falla)),
    #region = dat.spr$region,
    region = rep(1,nrow(dat.falla)),
    #period = dat.fall$period,
    period = rep(1,nrow(dat.falla)),
    length = dat.falla$LCat,
    mesh = dat.falla$mesh.cm,
    cpn = dat.falla$mean.cpue,
    rtype="norm.sca",
    distr = "poisson"         # poisson or nbinom
)
(par <- defpar(x))

# Run model 
(m1a <- gillnetfitTMB(x,par)) # this should say "convergence" at the end


# Get parameters
(par <- defpar(x))

# Run model 
(m1 <- gillnetfitTMB(x,par)) # this should say "convergence" at the end



### Check output 

# Tables
fittable(m1) # max.grad SHOULD BE roughly BELOW 0.0001
partable(m1) # there should be NO NaN in column sd
seltable(m1)
Ntable(m1)
selmaxtable(m1)

# Plots
meshconv <- data.frame(inch=c(2, 2.25, 2.5, 2.63, 2.75, 3),
                       mesh.lab=c("2 in","2 1/4 in","2 1/2 in","2 5/8 in","2 3/4 in","3 in"),
                       mesh.cm=c(2, 2.25, 2.5, 2.63, 2.75, 3)*2.54)

plotSel(m1,meshlabs = setNames(as.character(meshconv$mesh.lab), meshconv$mesh.cm))
plotN(m1)
plotRes(m1)
plotOP(m1)

# Fit all possible selectivity curves

rtypes <- c("norm.sca","norm.loc","lognorm",'gamma')
dists <- c("poisson","nbinom")

co <- expand.grid(rtypes,dists) 

ms <- apply(co[1:6,],1, function(i){
  x$rtype <- i[1]
  x$distr <- i[2]
  par <- defpar(x)
  gillnetfitTMB(x,par)
})
names(ms) <- do.call(paste, c(co[1:6,], sep=" - "))

ms <- do.call('c',ms) # equivalent to c(m1,m2,m3,m4)
ms
# Only norm.sca - nbinom and norm.loc - nbinom -> Convergence OK

AIC(ms)
partable(ms)

plotSel(ms,meshlabs = setNames(as.character(meshconv$mesh.lab), meshconv$mesh))
plotN(ms)
plotRes(ms)
plotOP(ms)

mypred <- seltable(ms)
mypred <- merge(mypred,meshconv)
mypred$mesh.lab <- factor(mypred$mesh.lab,levels=unique(meshconv$mesh.lab))
ggplot(mypred[mypred$fit==names(AIC(ms))[1],],aes(x=length,y=cpn))+
  geom_bar(stat='identity')+
  geom_line(data=mypred,aes(y=sel*max(cpn),col=fit))+ #selectivity rescaled to max cpn
  facet_grid(mesh.lab~.)
