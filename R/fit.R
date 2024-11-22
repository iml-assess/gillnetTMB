##' Default parameters
##' @param x input data (list)
##' @details Defaut parameters for gillnet selectivity function
##' @rdname defpar
##' @examples 
##' x <- defpar(data)
##' @export
defpar <-function(x){
    Ns <- nrow(unique(do.call("cbind",x[1:4])))
    mu <- mean(x$length/x$mesh)
    sd <- diff(range(x$length))/(length(unique(x$mesh))-1)
    switch(x$rtype, 
           norm.loc = {
               ret <- list(
                   k1=mu,         
                   k2=sd
               )
           }, # sd = k2 
           norm.sca = {
               ret <- list(
                   k1=mu,        
                   k2=(sd^2)/(mean(x$mesh)^2)
               )
           }, 
           gamma = {
               ret <- list(
                   k1=150,      
                   k2=0.01
               )
           }, # k1 = alpha
           lognorm = {
               ret <- list(
                   k1=4,          
                   k2=0.1
               )
           } 
    )
    
    di <- match.arg(x$distr,c("poisson","nbinom"))
    if(di=="nbinom") ret$logtheta <- log(3)
    
    ret$N <- rep(max(x$cpn),Ns) 
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
           norm.loc = {selfun <- function(length,mesh,k1,k2,...){exp(-(length-k1*mesh)^2/(2*k2^2))}}, # sd = k2 
           norm.sca = {selfun <- function(length,mesh,k1,k2,...){exp(-(length-k1*mesh)^2/(2*k2*mesh^2))}}, # sd ~ k2*mesh, spread is scaled to mesh
           gamma = {selfun <- function(length,mesh,k1,k2,...){((length/((k1-1)*k2*mesh))^(k1-1))*exp(k1-1-length/(k2*mesh))}}, # k1 = alpha
           lognorm = {selfun <- function(length,mesh,k1,k2,m1,...){(1/length)*exp(k1+log(mesh/m1)-k2^2/2-((log(length)-k1-log(mesh/m1))^2)/(2*k2^2))}} #  k1 = mu, k2 = sd
    )
    sel <- selfun(length,mesh,k1,k2,m1=min(mesh))
    
    # Prediction (can I get to logpred without all the for loops? tried with merge Nm -> didn't work)
    n <- length(cpn)
    logpred <- numeric(n)
    Nm <- cbind(unique(do.call("cbind",x[1:4])))

    for(y in unique(year)){
        for(r in unique(region)){
            for(p in unique(period)){
                for(m in unique(mesh)){
                    id <- which(mesh==m & year == y & region==r & period==p)
                    Nid <- which(Nm[,1]==y & Nm[,2]==r & Nm[,3]==p)
                    logpred[id] <- log(N[Nid])+log(sel[id])#+log(p[Nid]) 
                }
            }
        }
    }

    # likelihood
    pred <- exp(logpred)
    switch(distr, 
           poisson = {nll <- dpois(cpn,pred,TRUE)}, 
           nbinom = {nll <- dnbinom(x=cpn,prob=exp(logtheta)/(exp(logtheta)+pred),size=exp(logtheta),log=TRUE)}
    )
    
    # report
    res <- cpn - pred
    
    ADREPORT(res)
    ADREPORT(sel)
    ADREPORT(pred)
    
    # return
    -sum(nll)
}

