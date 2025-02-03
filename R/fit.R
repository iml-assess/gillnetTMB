##' Default parameters
##' @param x input data (list)
##' @details Defaut parameters for gillnet selectivity function. Based on code Surette 2016.
##' mu initiation: In logsel = -(length-k1 x mesh)^2/(2 x k2^2), at peak selectivity (logSel=0), the equation can be solved so k1 = length/mesh.
##' @rdname defpar
##' @examples 
##' x <- defpar(data)
##' @export
defpar <-function(x){

    mus <- aggregate(rep(x$length, round(x$cpn)), by = list(rep(x$mesh, round(x$cpn))), mean)
    sds <- aggregate(rep(x$length, round(x$cpn)), by = list(rep(x$mesh, round(x$cpn))), sd)
    
    switch(x$rtype, 
           norm.loc = {
               ret <- list(
                   k1=unname(coef(lm(mus[,2] ~ mus[,1] - 1))), # very similar to   mu <- mean(x$length/x$mesh) ,         
                   k2=log(mean(sds[, 2]))
               )
           }, # sd = k2 
           norm.sca = {
               ret <- list(
                   k1=unname(coef(lm(mus[,2] ~ mus[,1] - 1))),        
                   k2=unname(coef(lm(sds[,2] ~ sds[,1] - 1)))
               )
           }, 
           gamma = {
               xx <- x$length/ x$mesh
               m <- glm(x$cpn ~ xx + log(xx), family = "poisson")
               ret <- list(
                   k1=unname(coef(m)["log(xx)"] + 1),  # Gamma 'alpha' shape parameter.       
                   k2=unname(-1 / coef(m)["xx"])
               )
           }, # k1 = alpha
           lognorm = {
               mul <- mus[,2]+sds[,2]^2/2
               sdl <- sqrt((sds[,2]^2-1)+(2*mus[,2]+sds[,2]^2))
               ret <- list(
                   k1=unname(coef(lm(mul ~ mus[,1] - 1))),        
                   k2=unname(coef(lm(sdl ~ sds[,1] - 1))) # if log here:0.4
               )
           } 
    )
    
    di <- match.arg(x$distr,c("poisson","nbinom"))
    if(di=="nbinom") ret$logtheta <- 1
    
    temp <- aggregate(list(Ninit = x$cpn), by = x[1:4], mean)
    ret$N <- temp[,"Ninit"]
    
    check <- selTMB(ret,x)
    if (is.na(check) | !is.finite(check)) warning("Default initial parameters will result in undefined likelihood.")

    return(ret)
}

##' Fit gillnet selectivity model
##' @param x input data (list). mesh needs to be in increasing order.
##' @param par start parameters (list)
##' @param rtype type
##' @param distr distribution (poisson or dnbinom)
##' @details Fits gillnet selecitivty model (obs = N * sel)
##' @rdname defpar
##' @import RTMB
##' @export
gillnetfitTMB <- function(x,par){
    
    x$rtype <- match.arg(x$rtype,c("norm.sca","norm.loc","gamma","lognorm"))
    x$distr <- match.arg(x$distr,c("poisson","nbinom"))
    
    xx <- unique(do.call("cbind",x[1:4]))
    
    cmb <- function(f, d) function(p) f(p, d)
    obj <- MakeADFun(cmb(selTMB, x), par, silent=T)
    
    lower <- rep(0,length(obj$par))
    upper <- rep(Inf,length(obj$par))
    names(lower) <- names(upper) <- names(obj$par)
    upper[which(names(upper)=='N')] <- max(x$cpn)*1.01
    
    opt <- nlminb(obj$par, obj$fn, obj$gr,control=list(trace=1, eval.max=2000, iter.max=1000),lower=lower,upper=upper)
    
    sdr <- sdreport(obj)

    sel <- sdr$value[names(sdr$value)=="sel"]
    
    idN <- which(colnames(sdr$cov.fixed)=="N")
    Ncov <- sdr$cov.fixed[idN,idN]
    Nest <- data.frame(xx, est=sdr$par.fixed[idN],sd=sqrt(diag(Ncov)))

    ret <- list(sel = sel, 
                N = Nest, 
                conv = opt$message, 
                nll = opt$objective, 
                df = length(opt$par), 
                AIC = 2*length(opt$par)-2*opt$objective, 
                sdrep=sdr,  
                data = x, 
                para = as.list(sdr,"Est"),
                parasd = as.list(sdr,"Std"),
                opt=opt, 
                obj=obj)
    #attr(ret, "RemoteSha") <- substr(packageDescription("gillnetTMB")$RemoteSha, 1, 12)
    #attr(ret, "Version") <- packageDescription("gillnetTMB")$Version
    class(ret)<-"gillnet"
    return(ret)
}

##' Fit gillnet selectivity model
##' @param par parms
##' @param x data
##' @details RTMB model code
##' @rdname selTMB
selTMB <- function(par, x) {
    getAll(x, par)
    
    # Selectivity curve (Millar 1997&1999)
    switch(rtype, 
           norm.loc = {selfun <- function(length,mesh,k1,k2,...){-(length-k1*mesh)^2/(2*k2^2)}}, # mu=k1, sd = k2 
           norm.sca = {selfun <- function(length,mesh,k1,k2,...){-(length-k1*mesh)^2/(2*k2*mesh^2)}}, # mu=k1, sd = k2*mesh, spread is scaled to mesh
           gamma = {selfun <- function(length,mesh,k1,k2,...){(k1-1)*(log(length)-log((k1-1)*k2*mesh))+k1-1-length/(k2*mesh)}}, # k1 = alpha, k2=k (idem surette)
           lognorm = {selfun <- function(length,mesh,k1,k2,m1,...){log(1/length)+k1+log(mesh/m1)-k2^2/2-((log(length)-k1-log(mesh/m1))^2)/(2*k2^2)}} #  k1 = mu, k2 = sd
    )
    logsel <- selfun(length,mesh,k1,k2,m1=min(mesh))

    # Prediction
    n <- length(cpn)
    logpred <- numeric(n)
    
    d <- do.call("cbind",x[1:6])
    dn <- unique(do.call("cbind",x[1:4]))
    for(i in 1:n){
        r <- d[i,]
        id <- which(year == r[1] & region==r[2] & period==r[3] & mesh==r[5])
        Nid <- which(dn[,1]==r[1] & dn[,2]==r[2] & dn[,3]==r[3])
        logpred[id] <- log(N[Nid])+logsel[id]       #+log(p[Nid])
    }

    # likelihood
    pred <- exp(logpred)
    switch(distr, 
           poisson = {nll <- dpois(cpn,pred,TRUE)}, 
           nbinom = {nll <- dnbinom(x=cpn,prob=exp(logtheta)/(exp(logtheta)+pred),size=exp(logtheta),log=TRUE)}
    )
    
    # report
    sel <- exp(logsel)
    res <- cpn - pred
    
    ADREPORT(res)
    ADREPORT(sel)
    ADREPORT(pred)
    
    # return
    -sum(nll)
}

