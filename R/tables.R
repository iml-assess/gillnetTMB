##' combine.df
##' @param  x list of data.frames to which rbind can be applied
##' @details ...
combine.df <- function(x){
    n <- names(x)
    if(is.null(n)) n <- 1:length(x)
    do.call(rbind, lapply(n, function(name) {cbind(fit = name, x[[name]])}))
} 

##' fittable
##' @param  x...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
fittable <-function(x, ...){
    UseMethod("fittable")
}

##' @rdname fittable
##' @method fittable gillnet
##' @export
fittable.gillnet <- function(x){
    AIC <- AIC(x)
    nll <- logLik(x)
    df <- x$df
    max.grad <- max(-x$obj$gr(x$opt$par))
    coef <- coef(x)[1:2]
    ret <- c(AIC=AIC,df=df,nll=nll,max.grad=max.grad,coef)
    return(ret)
}

##' @rdname fittable
##' @method fittable gillnetset
##' @export
fittable.gillnetset <- function(x){
    ret <- sapply(x,fittable)
    nm <- names(x)
    if(!is.null(nm)) colnames(ret) <- nm
    return(ret)
}


##' partable
##' @param  x...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
partable <-function(x, ...){
    UseMethod("partable")
}

##' @rdname partable
##' @method partable gillnet
##' @export
partable.gillnet <- function(x, ...){
    ret <- data.frame(
        par=names(x$sdrep$par.fixed),
        est=x$sdrep$par.fixed,
        sd=sqrt(diag(x$sdrep$cov.fixed))
    )
    return(ret)
}

##' @rdname partable
##' @method partable gillnetset
##' @export
partable.gillnetset <- function(x, ...){
    ret <- lapply(x,partable)
    ret <- combine.df(ret)
    return(ret)
}


##' seltable
##' @param  x...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
seltable <-function(x, ...){
    UseMethod("seltable")
}

##' @rdname seltable
##' @method seltable gillnet
##' @export
seltable.gillnet <- function(x){
    sdr <- x$sdrep
    id <- names(sdr$value)=="logsel"
    logp <- sdr$value[id]
    logpsd <- sdr$sd[id]
    logp <- cbind(logp,logp+logpsd%o%c(-1.96,1.96))
    colnames(logp)<-c("estimate","low","high")
    ret <- data.frame(attr(x$data,"group"),length=x$data$length,mesh=x$data$mesh,exp(logp))
    rownames(ret) <- 1:nrow(ret)
    return(ret)
}

##' @rdname seltable
##' @method seltable gillnetset
##' @export
seltable.gillnetset <- function(x){
    ret <- lapply(x,seltable)
    ret <- combine.df(ret)
    return(ret)
}

##' predtable
##' @param  x...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
predtable <-function(x, ...){
    UseMethod("predtable")
}

##' @rdname predtable
##' @method predtable gillnet
##' @export
predtable.gillnet <- function(x){
    sdr <- x$sdrep
    id <- names(sdr$value)=="logpred"
    logp <- sdr$value[id]
    logpsd <- sdr$sd[id]
    logp <- cbind(logp,logp+logpsd%o%c(-1.96,1.96))
    colnames(logp)<-c("estimate","low","high")
    ret <- data.frame(attr(x$data,"group"),length=x$data$length,mesh=x$data$mesh,cpn=x$data$cpn,exp(logp))
    rownames(ret) <- 1:nrow(ret)
    return(ret)
}

##' @rdname predtable
##' @method predtable gillnetset
##' @export
predtable.gillnetset <- function(x){
    ret <- lapply(x,predtable)
    ret <- combine.df(ret)
    return(ret)
}

##' Ntable
##' @param  x...
##' @param ... extra arguments not currently used
##' @details ...
##' @importFrom reshape2 melt
##' @export
Ntable <-function(x, ...){
    UseMethod("Ntable")
}

##' @rdname Ntable
##' @method Ntable gillnet
##' @export
Ntable.gillnet <- function(x){
    sdr <- x$sdrep
    idN <- which(colnames(sdr$cov.fixed)=="logN")
    logN <- sdr$par.fixed[idN]
    logNsd <- sqrt(diag(sdr$cov.fixed[idN,idN]))
    logN <- cbind(logN,logN+logNsd%o%c(-1.96,1.96))
    colnames(logN)<-c("estimate","low","high")
    ret <- data.frame(unique(cbind(length=x$data$length,attr(x$data,"group"))),exp(logN))
    ret$distr <- x$data$distr
    ret$rtype <- x$data$rtype
    return(ret)
}

##' @rdname Ntable
##' @method Ntable gillnetset
##' @export
Ntable.gillnetset <- function(x){
    ret <- lapply(x,Ntable)
    ret <- combine.df(ret)
    return(ret)
}

##' selmaxtable
##' @param  x... gillnet or gillnetset object
##' @param ... extra arguments not currently used
##' @details ...
##' @export
selmaxtable <-function(x, ...){
    UseMethod("selmaxtable")
}

##' @rdname selmaxtable
##' @method selmaxtable gillnet
##' @export
selmaxtable.gillnet <- function(x){
    ret <- seltable(x)
    ret <- merge(aggregate(estimate ~ mesh, max, data = ret), ret[,c("length","mesh","estimate")])
    return(ret)
}

##' @rdname selmaxtable
##' @method selmaxtable gillnetset
##' @export
selmaxtable.gillnetset <- function(x){
    ret <- lapply(x,selmaxtable)
    ret <- combine.df(ret)
    return(ret)
}






