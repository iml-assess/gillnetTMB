##' plotN
##' @param  x...
##' @param ... extra arguments not currently used
##' @details ...
##' @import ggplot2
##' @export
plotN <-function(x, ...){
    UseMethod("plotN")
}

##' @rdname plotN
##' @param ylab ylab
##' @param xlab xlab
##' @param ci logical
##' @method plotN gillnet
##' @import ggplot2
##' @export
plotN.gillnet <- function(x,ylab="Relative abundance index",xlab="Length",ci=TRUE){
    d <- Ntable(x)
    p <- ggplot(d,aes(x=length,y=estimate))+
        geom_bar(stat="identity")+
        labs(y=ylab,x=xlab)+
        facet_wrap(year+region~period)
    if(ci) p <- p+geom_errorbar(aes(ymin=low,ymax=high),width=0.2)
    return(p)
}

##' @rdname plotN
##' @param ylab ylab
##' @param xlab xlab
##' @method plotN gillnetset
##' @import ggplot2
##' @export
plotN.gillnetset <- function(x,ylab="Relative abundance index",xlab="Length",ci=TRUE,collab="Fit",dodge=1.8){
    d <- Ntable(x)
    p <- ggplot(d,aes(x=length,y=estimate,fill=as.factor(fit),ymin=low,ymax=high))+
        geom_bar(stat="identity",aes(fill=as.factor(fit)),position=position_dodge())+
        labs(y=ylab,x=xlab,fill=collab)+
        facet_wrap(year+region~period)
    if(ci) p <- p+geom_errorbar(width=0.2,position=position_dodge(dodge))
    return(p)
}

##' plotSel
##' @param  x...
##' @param ... extra arguments not currently used
##' @details ...
##' @import ggplot2
##' @export
plotSel <-function(x, ...){
    UseMethod("plotSel")
}

##' @rdname plotSel
##' @param x gillnet object
##' @param ylab ylab
##' @param xlab xlab
##' @param collab color lab
##' @param meshlabs color labels
##' @method plotSel gillnet
##' @import ggplot2 viridis
##' @export
plotSel.gillnet <- function(x,ylab="Selectivity",xlab="Length",collab="Mesh",meshlabs=NULL){
    d <- seltable(x)
    m <- selmaxtable(x)
    if(is.null(meshlabs)) meshlabs <- setNames(as.character(unique(d$mesh)), unique(d$mesh)) 
    ggplot(d,aes(x=length,y=estimate))+
        #geom_vline(data=m,aes(xintercept=peak),col="grey",linetype="dashed")+
        geom_ribbon(aes(ymin=low,ymax=high,fill=as.factor(mesh)),alpha=0.5)+
        geom_line(aes(col=as.factor(mesh)))+
        labs(y=ylab,x=xlab,col=collab,fill=collab)+
        scale_color_viridis_d(labels=meshlabs)+
        scale_fill_viridis_d(labels=meshlabs)+
        scale_y_continuous(limits=c(0,1.05),expand=c(0,0))
}

##' @rdname plotSel
##' @param x gillnetset object
##' @param ylab ylab
##' @param xlab xlab
##' @param collab color lab
##' @param meshlabs color labels
##' @method plotSel gillnetset
##' @import ggplot2 viridis
##' @export
plotSel.gillnetset <- function(x,ylab="Selectivity",xlab="Length",collab="Mesh",meshlabs=NULL){
    d <- seltable(x)
    m <- selmaxtable(x)
    if(is.null(meshlabs)) meshlabs <- setNames(as.character(unique(d$mesh)), unique(d$mesh)) 
    ggplot(d,aes(x=length,y=estimate))+
        geom_vline(data=m,aes(xintercept=peak),col="grey",linetype="dashed")+
        geom_ribbon(aes(ymin=low,ymax=high,fill=as.factor(mesh)),alpha=0.5)+
        geom_line(aes(col=as.factor(mesh)))+
        labs(y=ylab,x=xlab,col=collab,fill=collab)+
        scale_color_viridis_d(labels=meshlabs)+
        scale_fill_viridis_d(labels=meshlabs)+
        facet_wrap(~fit,ncol=1)+
        scale_y_continuous(limits=c(0,1.05),expand=c(0,0))+
        scale_x_continuous(limits=range(d$length))
}

##' plotRes
##' @param  x...
##' @param ... extra arguments not currently used
##' @details ...
##' @import ggplot2
##' @export
plotRes <-function(x, ...){
    UseMethod("plotRes")
}

##' @rdname plotRes
##' @param ylab ylab
##' @param xlab xlab
##' @method plotRes gillnet
##' @import ggplot2
##' @export
plotRes.gillnet <- function(x,xlab="Length",ylab="Mesh size"){
    d <- residuals(x)
    d$sign <- ifelse(d$res<0,"-","+")
    ggplot(d,aes(x=length,y=as.factor(mesh),size=res^2,col=sign))+
        geom_point()+
        #theme(legend.position = "none")+
        labs(x=xlab,y=ylab)+
        facet_wrap(year+region~period)
}

##' @rdname plotRes
##' @param ylab ylab
##' @param xlab xlab
##' @method plotRes gillnetset
##' @import ggplot2
##' @export
plotRes.gillnetset <- function(x,ylab="Selectivity",xlab="Length",collab="Mesh"){
    d <- residuals(x)
    d$sign <- ifelse(d$res<0,"-","+")
    ggplot(d,aes(x=length,y=as.factor(mesh),size=res^2,col=sign))+
        geom_point()+
        theme(legend.position = "none")+
        labs(x=xlab,y=ylab)+
        facet_grid(year+region+period~fit)
}

##' plot observed vs predicted
##' @param  x...
##' @param ... extra arguments not currently used
##' @details ...
##' @import ggplot2
##' @export
plotOP <-function(x, ...){
    UseMethod("plotOP")
}

##' @rdname plotOP
##' @param  x gillnet object
##' @param log logical
##' @param ylab ylab
##' @param xlab xlab
##' @param collabl color lab
##' @param facet plot mesh as facet
##' @param meshlabs labels for mesh facets (if facet=TRUE)
##' @method plotOP gillnet
##' @import ggplot2 viridis
##' @export
plotOP.gillnet <- function(x, log=FALSE, xlab="Observed", ylab="Predicted", collab="Mesh", facet=FALSE, meshlabs=NULL, ...){
    d <- predtable(x)
    if(log){
        d$cpn <- log(d$cpn)
        d$estimate <- log(d$estimate)
    }
    p <- ggplot(d,aes(x=cpn,y=estimate,col=as.factor(mesh)))+
        geom_point()+
        geom_abline(intercept = 0,slope=1,col="darkgrey")+
        #theme(legend.position = "none")+
        labs(x=xlab,y=ylab,col=collab)+
        scale_color_viridis_d()
     if(facet){
         meshlabs <- ifelse(is.null(meshlabs),unique(d$mesh),meshlabs)
         p <- p + facet_wrap(~mesh,labeller=as_labeller(meshlabs))
     }  
    p
}

##' @rdname plotOP
##' @param  x gillnetset object
##' @param log logical
##' @param ylab ylab
##' @param xlab xlab
##' @param collabl color lab
##' @method plotOP gillnetset
##' @import ggplot2 viridis
##' @export
plotOP.gillnetset <- function(x, log=FALSE, xlab="Observed", ylab="Predicted", collab="Mesh", ...){
    d <- predtable(x)
    if(log){
        d$cpn <- log(d$cpn)
        d$pred <- log(d$pred)
    }
    ggplot(d,aes(x=cpn,y=estimate,col=as.factor(mesh)))+
        geom_point()+
        geom_abline(intercept = 0,slope=1,col="darkgrey")+
        #theme(legend.position = "none")+
        labs(x=xlab,y=ylab,col=collab)+
        scale_color_viridis_d()+
        facet_wrap(~fit)
}

