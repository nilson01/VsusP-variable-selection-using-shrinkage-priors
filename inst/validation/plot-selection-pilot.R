# Make the pilot figure from the complete summary table; base R only.
plot_selection_pilot <- function(summary.file, output.file) {
  tab <- utils::read.csv(summary.file,check.names=FALSE)
  methods <- c("SAVS","S2M_legacy_b0.5","S2M_exact_refine_b0.5")
  colors <- c("#0072B2","#D55E00","#555555"); symbols <- c(16,15,17)
  grDevices::pdf(output.file,width=6.7,height=3.7,pointsize=9,useDingbats=FALSE)
  on.exit(grDevices::dev.off(),add=TRUE)
  graphics::layout(matrix(c(1,1,2,3),2,2,byrow=TRUE),heights=c(.16,1))
  graphics::par(mar=rep(0,4));graphics::plot.new()
  graphics::legend("center",c("SAVS","S2M, b = 0.5","S2M exact/refined, b = 0.5"),
                   col=colors,pch=symbols,horiz=TRUE,bty="n",cex=.85)
  for(k in seq_along(c("empty","exact"))) {
    metric <- c("empty","exact")[k]
    graphics::par(mar=c(4.5,3.6,2,0.8),mgp=c(2.2,.65,0),tcl=-.22)
    graphics::plot(NA,xlim=c(.5,6.5),ylim=c(-.025,1.025),xaxt="n",
                   xlab="",ylab="Proportion of datasets",bty="l")
    graphics::abline(h=c(0,.5,1),col="#DDDDDD",lwd=.6)
    labels <- lapply(seq_len(6),function(i) {
      r <- tab[tab$cell==i,][1,]
      bquote(atop(.(r$effect),.(r$rho)))
    })
    graphics::axis(1,at=1:6,labels=as.expression(labels),cex.axis=.85)
    graphics::mtext(expression(paste("Top: ",beta[1]==-beta[2],"; bottom: ",rho)),1,line=3,cex=.85)
    graphics::title(main=c("a  Empty selection","b  Exact support recovery")[k],adj=0,cex.main=.95,font.main=2)
    for(j in seq_along(methods)) {
      rows <- tab[tab$method==methods[j],];rows <- rows[order(rows$cell),]
      stopifnot(nrow(rows)==6L,all(rows$n_success==10))
      x <- (1:6)+(j-2)*.18
      graphics::segments(x,rows[[paste0(metric,".lower")]],x,rows[[paste0(metric,".upper")]],col=colors[j],lwd=1)
      graphics::points(x,rows[[paste0(metric,".mean")]],pch=symbols[j],col=colors[j],cex=.8)
    }
  }
}
