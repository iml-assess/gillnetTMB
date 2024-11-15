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
##' @method plotN gillnet
##' @import ggplot2
##' @export
plotN.gillnet <- function(x,ylab="Relative abundance index",xlab="Length"){
    d <- Ntable(x)
    ggplot(d,aes(x=length,y=est))+
        geom_bar(stat="identity")+
        geom_errorbar(aes(ymin=est-2*sd,ymax=est+2*sd),width=0.2)+
        labs(y=ylab,x=xlab)+
        facet_wrap(year+region~period)
}

##' @rdname plotN
##' @param ylab ylab
##' @param xlab xlab
##' @method plotN gillnet
##' @import ggplot2
##' @export
plotN.gillnetset <- function(x,ylab="Relative abundance index",xlab="Length",filllab="Fit",dodge=1.8){
    d <- Ntable(x)
    ggplot(d,aes(x=length,y=est,fill=as.factor(fit),ymin=est-1.96*sd,ymax=est+1.96*sd))+
        geom_bar(stat="identity",aes(fill=as.factor(fit)),position=position_dodge())+
        geom_errorbar(width=0.2,position=position_dodge(dodge))+
        labs(y=ylab,x=xlab,fill=filllab)+
        facet_wrap(year+region~period)
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
##' @param ylab ylab
##' @param xlab xlab
##' @method plotSel gillnet
##' @import ggplot2 viridis
##' @export
plotSel.gillnet <- function(x,ylab="Selectivity",xlab="Length",collab="Mesh"){
    d <- seltable(x)
    ggplot(d,aes(x=length,y=sel))+
        #geom_ribbon(aes(ymin=predsel-2*predsel.sd,ymax=predsel+2*predsel,fill=as.factor(mesh)),alpha=0.5)+
        geom_line(aes(col=as.factor(mesh)))+
        labs(y=ylab,x=xlab,col=collab)+
        scale_color_viridis_d()
}

##' @rdname plotSel
##' @param ylab ylab
##' @param xlab xlab
##' @method plotSel gillnet
##' @import ggplot2
##' @export
plotSel.gillnetset <- function(x,ylab="Selectivity",xlab="Length",collab="Mesh"){
    d <- seltable(x)
    ggplot(d,aes(x=length,y=sel))+
        #geom_ribbon(aes(ymin=predsel-2*predsel.sd,ymax=predsel+2*predsel,fill=as.factor(mesh)),alpha=0.5)+
        geom_line(aes(col=as.factor(mesh)))+
        labs(y=ylab,x=xlab,col=collab)+
        scale_color_viridis_d()+
        facet_wrap(~fit,ncol=1)
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
##' @import ggplot2 viridis
##' @export
plotRes.gillnet <- function(x,xlab="Length",ylab="Mesh size"){
    d <- data.frame(residuals(x))
    d$sign <- ifelse(d$res<0,"-","+")
    ggplot(d,aes(x=length,y=as.factor(mesh),size=res^2,col=sign))+
        geom_point()+
        #theme(legend.position = "none")+
        labs(x=xlab,y=ylab)
}

##' @rdname plotRes
##' @param ylab ylab
##' @param xlab xlab
##' @method plotRes gillnet
##' @import ggplot2
##' @export
plotRes.gillnetset <- function(x,ylab="Selectivity",xlab="Length",collab="Mesh"){
    d <- data.frame(residuals(x))
    d$sign <- ifelse(d$res<0,"-","+")
    ggplot(d,aes(x=length,y=as.factor(mesh),size=res^2,col=sign))+
        geom_point()+
        theme(legend.position = "none")+
        labs(x=xlab,y=ylab)+
        facet_wrap(~fit)
}
