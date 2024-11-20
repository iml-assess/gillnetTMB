.First <- function(){
  
  ### load packages
  cran.packages <- c('reshape2','ggplot2','gridExtra','viridis','plyr')
  new.packages <- cran.packages[!(cran.packages %in% utils::installed.packages()[,"Package"])]
  if(length(new.packages)>0) install.packages(new.packages)
  
  invisible(lapply(c(cran.packages), function(x) require(x, character.only = TRUE)))
  
  ### source R directory
  srcf <- unlist(sapply("./R",list.files,full.names=T))
  invisible(sapply(srcf,source))
  
  ### ggplot layout
  theme_set(theme_classic())
}
