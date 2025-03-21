# GillnetTMB

Model to estimate gear selectivity and population abundance by factor (year/period/region). Selectivity is assumed constant over space and time.

The model follows;

logPred(mlf) = logN(lf) + logSel(l) (based on Surette et al. 2016),

where logPred(mlg) is the predicted mean for mesh size m and fish length l, for factor f. logN is the estimated population abundance corrected for selectivity, on a log scale. logSel is time time- and space invariant gear selectivity at length. The package has four types of selectivity curves (see Millar 1997 and 1999). Two count distributions can be used; the Poisson and the negative binomial. Note that each mesh size is assumed to have equal power (i.e., their relative fishing power  is one).

The model is based on RTMB. 

# Installation

```
devtools::install_github("iml-assess/gillnetTMB")
```

The package requires installation of Rtools.

# Example

See example files. In short:

```
x <- as.list(mydata)
par <- defpar(x)
m <- gillnetfitTMB(x,par) 
m
```

Tables and plots (including generic methods such as AIC, logLik, coef):
```
fittable(m)
partable(m)
AIC(m)
seltable(m)
Ntable(m)
selmaxtable(m)

plotSel(m)
plotN(m)
plotRes(m)
plotOP(m)
```

# Model validation

1) model convergence (print m, or m$opt$convergence==0)
2) maximum gradient <0.001 (see fittable)
3) all parametesr have a standard error (see partable)

# Notes

Advantages over the TropFishR package, using Millar's original gillnet selectivity fitting function (function "gillnetfit"):
1) Estimation of abundance for various factors (combinations of year-region-period) assuming constant selectivity.
2) Confidence intervals around the abundance and selectivity estimates are readily provided.
3) Two error distributions (posibility to explore the negative binomial).

# References

1. Millar, R., 1997. Estimation of gillnet and hook selectivity using log-linear models. ICES J. Mar. Sci. 54, 471–477. https://doi.org/10.1006/jmsc.1996.0196
2. Millar, R.B., Freyer, R.J., 1999. Estimating the size-selection curves of towed gears, traps, nets and hooks. Rev. Fish Biol. Fish. 9, 89–116. https://doi.org/https://doi.org/10.1023/A:1008838220001
3. Surette, T.., LeBlanc, C.., Mallet, A., 2016. Abundance indices and selectivity curves from experimental multi-panel gillnets for the southern Gulf of St. Lawrence fall herring fishery. Can. Sci. Advis. Secr. Res. Doc. 067, vi + 23 p.