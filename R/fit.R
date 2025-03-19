##' Default parameters
##' @param x input data (list)
##' @param plot plot selectivity curves at inital parameter estimates.
##' @details Defaut parameters for gillnet selectivity function. Based on code Surette 2016.
##' mu initiation: In logsel = -(length-k1 x mesh)^2/(2 x k2^2), at peak selectivity (logSel=0), the equation can be solved so k1 = length/mesh.
##' @rdname defpar
##' @import ggplot2
##' @examples 
##' x <- defpar(data)
##' @export
defpar <-function(x,plot=FALSE){
    
    xx <- defdat(x)
    mus <- aggregate(rep(xx$length, round(xx$cpn)), by = list(rep(xx$mesh, round(xx$cpn))), mean)
    sds <- aggregate(rep(xx$length, round(xx$cpn)), by = list(rep(xx$mesh, round(xx$cpn))), sd)
    
    switch(x$rtype, 
           norm.loc = {
               ret <- list(par=c(
                   k1=unname(coef(lm(mus[,2] ~ mus[,1] - 1))), # very similar to   mu <- mean(x$length/x$mesh) ,         
                   k2=mean(sds[, 2]) # surette: log(mean)
               ))
           }, # sd = k2 
           norm.sca = {
               ret <- list(par=c(
                   k1=unname(coef(lm(mus[,2] ~ mus[,1] - 1))),        
                   k2=unname(coef(lm(sds[,2] ~ sds[,1] - 1)))
               ))
           }, 
           gamma = {
               r <- xx$length/ xx$mesh
               m <- glm(xx$cpn ~ r + log(r), family = "poisson")
               ret <- list(par=c(
                   k1=unname(coef(m)["log(r)"] + 1),  # Gamma 'alpha' shape parameter.       
                   k2=unname(-1 / coef(m)["r"])
               ))
           }, # k1 = alpha
           lognorm = {
               mul <- mus[,2]+sds[,2]^2/2
               sdl <- sqrt((sds[,2]^2-1)+(2*mus[,2]+sds[,2]^2))
               ret <- list(par=c(
                   k1=unname(coef(lm(mul ~ mus[,1] - 1))),        
                   k2=unname(coef(lm(sdl ~ sds[,1] - 1))) # if log here:0.4
               ))
           } 
    )
    
    if(xx$distr==2) ret$par <- c(ret$par,logtheta=1)
    
    temp <- aggregate(list(Ninit = xx$cpn), by = xx[c(1,6)], mean)
    temp$Ninit[temp$Ninit==0] <- min(temp$Ninit[temp$Ninit!=0])
    ret$logN <- matrix(log(temp[,"Ninit"]),ncol=length(unique(temp$length))) # check if multiple groups if byrow or not
    dimnames(ret$logN) <- list(group=1:nrow(ret$logN),length=unique(temp$length))

    check <- gillnetR(ret,xx)
    if (is.na(check) | !is.finite(check)) warning("Default initial parameters will result in undefined likelihood.")
    
    if(plot){
        df <- data.frame(do.call(cbind,x[1:5]))
        df$selectivity <- do.call(predSel,c(x,ret$par))
        p <- ggplot(df,aes(x=length,y=selectivity,col=as.factor(mesh)))+
            geom_line()+
            labs(title="Initial parameters values")+
            scale_y_continuous(expand=c(0,0),limits=c(0,1))
        print(p)
    }
    
    return(ret)
}

##' Data validation
##' @param x input data (list)
##' @param warning logical
##' @details Transforms data input, for easy integration in TMB. Lengths are removed for which only zeros were observed across all mesh sizes, years, regions and periods (see Surette).
##' @rdname defdat
##' @examples 
##' x <- defdat(data)
defdat <-function(x,warning=FALSE){
    
    if(missing(x)) stop("argument 'x' is missing, with no default")
    rtypes <- c("norm.loc","norm.sca","gamma","lognorm")
    distrs <- c("poisson","nbinom")
    x$rtype <- match.arg(x$rtype,rtypes)
    x$distr <- match.arg(x$distr,distrs)
    
    # Remove lengths for which always zero (idem Surette)
    checkzero <- apply(do.call("cbind",x[1:4]),1,paste,collapse=".")
    checkzero <- as.numeric(as.factor(checkzero))
    temp <- aggregate(x$cpn, by = list(checkzero), sum)
    ix <- temp[,2][checkzero]==0
    if(any(ix)){
        if(warning) print(paste("Removed", length(ix[ix]), "observations with zero counts for all mesh sizes.") )
        x[1:6] <- lapply(x[1:6],function(i)i[!ix])
    }
    
    # sort
    s <- order(x[[1]], x[[2]],x[[3]],x[[5]],x[[4]])
    x[c(1:5)] <- lapply(x[c(1:5)],function(x) x[s])

    # Groups
    group <- apply(do.call("cbind",x[1:3]),1,paste,collapse=".")
    group <- as.numeric(as.factor(group))
    
    # Transform input data for cpp
    ret <- list(length=x$length,
                 mesh=x$mesh,
                 cpn=x$cpn,
                 rtype=which(rtypes==x$rtype),
                 distr=which(distrs==x$distr),
                 group=group,
                 lengthid=as.numeric(as.factor(x$length)))
    
    return(ret)
}

##' Fit gillnet selectivity model
##' @param x input data (list). mesh needs to be in increasing order.
##' @param par start parameters (list)
##' @details Fits gillnet selecitivty model (obs = N * sel)
##' @rdname gillnetfit
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats nlminb
##' @useDynLib gillnetTMB
##' @examples
##' fit <- gillnetfit(data, defpar(par))
##' @export
gillnetfit <- function(x,par){
    
    # Checks
    if(missing(x)) stop("argument 'x' is missing, with no default")
    if(missing(par))  par <- defpar(x)
    xx <- defdat(x,TRUE)
    
    obj <- MakeADFun(xx,par,DLL = "gillnetTMB") #compile('src/gillnetTMB.cpp');dyn.load(dynlib("src/gillnetTMB"))

    lower <- rep(0,length(obj$par))
    upper <- rep(Inf,length(obj$par))
    names(lower) <- names(upper) <- names(obj$par)
    upper[grep("N",names(upper))] <- max(x$cpn)*1.01
    
    opt <- nlminb(obj$par,obj$fn,obj$gr,upper=upper,lower=lower,
                  control = list(trace=10,eval.max=2000,iter.max=1000))
    
    sdr <- sdreport(obj)
    
    ret <- list(conv = opt$message, 
                nll = opt$objective, 
                df = length(opt$par), 
                AIC = 2*length(opt$par)-2*opt$objective, 
                sdrep = sdr,  
                data = x, 
                para = as.list(sdr,"Est"),
                parasd = as.list(sdr,"Std"),
                opt = opt, 
                obj = obj)
    
    attr(ret, "Version") <- packageDescription("gillnetTMB")$Version
    class(ret)<-"gillnet"
    return(ret)
}

##' Fit gillnet selectivity model
##' @param par parms
##' @param x data
##' @details Gillnet model code (R version). Used to validate initial parameter estimation (defpar). Note: does not compile with RTMB (complex number conversion).
##' @rdname gillnetR
gillnetR <- function(par, x) {
    list2env(x, env = environment())     # getAll(x, par) for RTMB
    list2env(par, env = environment())
    
    # Selectivity curve (Millar 1997&1999)
    switch(rtype, 
           norm.loc = {selfun <- function(length,mesh,k1,k2,...){-(length-k1*mesh)^2/(2*k2^2)}}, # mu=k1, sd = k2 
           norm.sca = {selfun <- function(length,mesh,k1,k2,...){-(length-k1*mesh)^2/(2*k2*mesh^2)}}, # mu=k1, sd = k2*mesh, spread is scaled to mesh
           gamma = {selfun <- function(length,mesh,k1,k2,...){(k1-1)*(log(length)-log((k1-1)*k2*mesh))+k1-1-length/(k2*mesh)}}, # k1 = alpha, k2=k (idem surette)
           lognorm = {selfun <- function(length,mesh,k1,k2,m1,...){log(1/length)+k1+log(mesh/m1)-k2^2/2-((log(length)-k1-log(mesh/m1))^2)/(2*k2^2)}} #  k1 = mu, k2 = sd
    )
    logsel <- selfun(length,mesh,par[1],par[2],m1=min(mesh))
   
    # Prediction
    n <- length(cpn)
    logpred <- numeric(n)
    
    for(i in 1:n){
        g <- group[i]
        l <- lengthid[i]
        logpred[i] <- logN[g,l]+logsel[i]
    }
    pred <- exp(logpred)
    
    # likelihood
     switch(distr, 
           poisson = {nll <- dpois(cpn,pred,TRUE)}, 
           nbinom = {nll <- dnbinom(x=cpn,prob=exp(par[3])/(exp(par[3])+pred),size=exp(par[3]),log=TRUE)}
    )
    
    # return
    -sum(nll)
}

