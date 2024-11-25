##' Print gillnet object
##' @method print gillnet
##' @param  x ...
##' @param  ... extra arguments
##' @details ...
##' @export
print.gillnet<-function(x, ...){
    cat("Gillnet selectivity model: log likelihood is", x$nll,"Convergence", ifelse(0==x$opt$convergence, "OK\n", "failed\n"))
}

##' Log likelihood of gillnet object
##' @method logLik gillnet
##' @param  object gillnet fitted object (result from gillnetfitTMB)
##' @param  ... extra arguments
##' @details ...
##' @export
logLik.gillnet<-function(object, ...){
    ret<- object$nll
    attr(ret,"df")<-object$df
    class(ret)<-"logLik"
    ret
}

##' AIC of gillnet object
##' @method AIC gillnet
##' @param  object gillnet fitted object (result from gillnetfitTMB)
##' @param  ... extra arguments
##' @details ...
##' @export
AIC.gillnet<-function(object, ...){
    ret<- object$AIC
    class(ret)<-"AIC"
    ret
}

##' AIC of gillnetset object
##' @method AIC gillnetset
##' @param  object gillnet fitted object (result from gillnetfitTMB)
##' @param  ... extra arguments
##' @details ...
##' @export
AIC.gillnetset<-function(object, ...){
    ret<- sort(sapply(object,AIC))
    class(ret)<-"AIC"
    ret
}

##' Extract fixed selectivity coefficients of gillnet object
##' @method coef gillnet
##' @param  object gillnet fitted object (result from illnetfitTMB)
##' @param  ... extra arguments
##' @details ...
##' @importFrom stats coef
##' @export
coef.gillnet <- function(object, ...){
    ret <- object$sdrep$par.fixed
    attr(ret,"cov") <- object$sdrep$cov.fixed
    attr(ret,"sd") <- sqrt(diag(object$sdrep$cov.fixed))
    class(ret)<-"gillnetcoef"
    ret
}

##' Print gillnetcoef object
##' @method print gillnetcoef
##' @param  x ...
##' @param  ... extra arguments
##' @details ...
##' @export
print.gillnetcoef<-function(x, ...){
    y<-as.vector(x)
    names(y)<-names(x)
    print(y)
}

##' Extract residuals from gillnet object
##' @method residuals gillnet
##' @param object gillnet fitted object (result from gillnet.fit)
##' @importFrom stats residuals
##' @details ...
##' @export
residuals.gillnet<-function(object){
    ret <- do.call('cbind',object$data[1:5])
    ret <- cbind(ret,res=object$sdrep$value[names(object$sdrep$value)=="res"])
    ret
}

##' Extract residuals from gillnetset object
##' @method residuals gillnetset
##' @param object gillnetset fitted objects (results from gillnet.fit)
##' @importFrom stats residuals
##' @details ...
##' @export
residuals.gillnetset<-function(object){
    ret <- lapply(object,function(x) data.frame(residuals(x)))
    ret <- combine.df(ret)
    return(ret)
}

##' Collect gillnet objects
##' @method c gillnet
##' @param  ... gillnet fits to be combined
##' @details ...
##' @export
c.gillnet<-function(...){
    ret<-list(...)
    class(ret)<-"gillnetset"
    ret
}
