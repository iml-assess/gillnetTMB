# GillnetTMB

Get selectivity and population index estimates, by year/period/region. Selectivity is assumed constant over space and time.

logObs = logN + logSel 

1) 4 selectivity curves
2) poisson or negative binomial distribution

Based on RTMB.

# Installation

Currently not packaged (just Rproject)
devtools::install_github("iml-assess/gillnetTMB")

# Example

See example file. In short:

x <- as.list(mydata)

par <- defpar(x)

m <- gillnetfitTMB(x,par) 

m

partable(m)

fittable(m)

seltable(m)

Ntable(m)

plotSel(m)

plotN(m)

plotRes(m)
