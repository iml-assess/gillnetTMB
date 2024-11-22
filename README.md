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


partable(m) # parameters

fittable(m) # model output

seltable(m) # selectivity estimates

selmaxtable(m) # lengths at peak selectivity

Ntable(m) # abundance estimates


plotSel(m) # plot selectivity curves

plotN(m) # plot abundance

plotRes(m) # plot residuals (bubbles)

plotOP(m) # plot observed versus predicted

# model validation

1) model convergence (print m, or m$opt$convergence==0)
2) maximum gradient <0.001 (see fittable)
3) all parametesr have a standard error (see partable)
