ggplotlimdil <- function(x, col.group=NULL, cex=1, lwd=1, legend.pos="bottomleft", ...)
  #	Plot method for limdil objects
  #	Yifang Hu and Gordon Smyth
  #	20 February 2009.  Last revised 6 February 2013.
{
  x$group <- factor(x$group)
  num.group <- nlevels(x$group)
  if(is.null(col.group)) 
    col.group <- 1:num.group
  else
    col.group <- rep(col.group,num.group)
  
  col <- x$group
  levels(col) <- col.group
  col <- as.character(col)
  dose <- x$dose
  maxx <- max(dose)	
  
  i <- x$response==x$tested
  x$response[i] <- x$response[i]-0.5
  
  nonres <- log(1-x$response/x$tested)
  if(num.group>1 && any(i)) nonres <- pmin(0,jitter(nonres))
  
  miny <- min(nonres)
  # par(mar = c(3, 4, 1, 1))
  # plot(x=1,y=1,xlim=c(0,maxx),ylim=c(min(miny,-0.5),0),xlab="dose (number of cells)",
  #      ylab="log fraction nonresponding",type="n",las = 1,...)
  # points(dose[!i],nonres[!i],pch=1,col=col[!i],cex=cex)
  # points(dose[i],nonres[i],pch=6,col=col[i],cex=cex)
  # 
  p <- ggplot() +
        geom_point(aes(dose[!i],nonres[!i]), col = as.factor(col[!i]), shape = as.factor(isTRUE(i))) +
        geom_point(aes(dose[i],nonres[i]), col = col[i]) + xlim(0,maxx) + ylim(min(miny,-0.5),0) +
        geom_abline(intercept = 0, slope=-1/x$CI[1,2],col=col.group[1]) +
        geom_abline(intercept = 0, slope=-1/x$CI[1,1],col=col.group[1],linetype = "dashed") +
        geom_abline(intercept = 0, slope=-1/x$CI[1,3],col=col.group[1],linetype = "dashed") +
        geom_abline(intercept = 0, slope=-1/x$CI[2,2],col=col.group[2]) +
        geom_abline(intercept = 0, slope=-1/x$CI[2,1],col=col.group[2],linetype = "dashed") +
        geom_abline(intercept = 0, slope=-1/x$CI[2,3],col=col.group[2],linetype = "dashed") +
        ylab("log fraction without spheres") + xlab("dose (number of cells)") +
        theme(axis.title.y = element_text(size = 12)) + theme_cowplot()

  # for(g in 1:num.group) {
  #   abline(a=0,b=-1/x$CI[g,2],col=col.group[g],lty=1,lwd=lwd)
  #   abline(a=0,b=-1/x$CI[g,1],col=col.group[g],lty=2,lwd=lwd)
  #   abline(a=0,b=-1/x$CI[g,3],col=col.group[g],lty=2,lwd=lwd)
  # }
  # if(num.group>1) legend(legend.pos,legend=paste("Group",levels(x$group)),text.col=col.group,cex=0.6*cex)
  # invisible(list(x=dose,y=nonres,group=x$group))
  return(p)
}

elda_bar <- function (df, color.group){
  df <- as.data.frame(100/df$CI)
  df$Group <- c("-Dox","+Dox")
  p <- df %>% 
    ggplot(aes(Group, Estimate, fill = Group)) +
    geom_col(position = "dodge") +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), position = "dodge", width = 0.25) + 
    scale_y_continuous(expand = c(0, 0))+
    theme_cowplot() + font_size + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), legend.position = "none") +
    scale_fill_manual(values = color.group) 
  print(p)
}
