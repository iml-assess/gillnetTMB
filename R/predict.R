##' Predict selectivity
##' @param length vector
##' @param mesh vector
##' @param k1 param
##' @param k2 param
##' @param rtype type of selectivity curve
##' @details Predicts selectivity based on knowm parameter values
##' @rdname predSel
##' @export
predSel <- function(length,mesh,k1,k2,rtype="norm.loc"){
    rtype <- match.arg(rtype,c("norm.sca","norm.loc","gamma","lognorm"))
    
    switch(rtype, 
           norm.loc = {selfun <- function(length,mesh,k1,k2,...){-(length-k1*mesh)^2/(2*k2^2)}}, # mu=k1, sd = k2 
           norm.sca = {selfun <- function(length,mesh,k1,k2,...){-(length-k1*mesh)^2/(2*k2*mesh^2)}}, # mu=k1, sd = k2*mesh, spread is scaled to mesh
           gamma = {selfun <- function(length,mesh,k1,k2,...){(k1-1)*(log(length)-log((k1-1)*k2*mesh))+k1-1-length/(k2*mesh)}}, # k1 = alpha, k2=k (idem surette)
           lognorm = {selfun <- function(length,mesh,k1,k2,m1,...){log(1/length)+k1+log(mesh/m1)-k2^2/2-((log(length)-k1-log(mesh/m1))^2)/(2*k2^2)}} #  k1 = mu, k2 = sd
    )
    logsel <- selfun(length,mesh,k1,k2,m1=min(mesh))
    
    return(exp(logsel))
}
