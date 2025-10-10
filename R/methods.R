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
    ret <- predtable(object)
    ret$res <- with(ret,cpn-estimate)
    ret
}

##' Extract residuals from gillnetset object
##' @method residuals gillnetset
##' @param object gillnetset fitted objects (results from gillnet.fit)
##' @importFrom stats residuals
##' @details ...
##' @export
residuals.gillnetset<-function(object){
    ret <- lapply(object,function(x) residuals(x))
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
    return(ret)
}

##' Predict selectivity
##' @method predict gillnet
##' @param object gillnet fit object
##' @param  ... 
##' @details ...
##' @export
predict.gillnet<-function(object,length=seq(min(object$data$length), max(object$data$length), by = 0.01),...){
    rtypes <- c("norm.loc","norm.sca","gamma","lognorm")
    k1 <- object$para$par[1]
    k2 <- object$para$par[2]
    r <-  rtypes[object$data$rtype]   
    
    mesh <- unique(object$data$mesh)
    ret <- expand.grid(length=length,mesh=mesh)
    ret$y <- predSel(ret$length,ret$mesh,k1,k2,r)
    return(ret)
}

##' Predict selectivity
##' @method predict gillnetset
##' @param object gillnetset object
##' @param  ... 
##' @details ...
##' @export
predict.gillnetset <- function(object,...){
    ret <- lapply(object,function(x) predict(x,...))
    ret <- combine.df(ret)
    return(ret)
}

