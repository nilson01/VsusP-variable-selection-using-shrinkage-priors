# Rebuild the documentation figures from the supplied numerical inputs.
main<-function(root=commandArgs(trailingOnly=TRUE)[1]) {
  if(is.na(root))stop("Supply the package source directory.")
  library(VsusP)
  root<-normalizePath(root,mustWork=TRUE)
  data<-file.path(root,"inst","doc-data");out<-file.path(root,"man","figures")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  read<-function(name)utils::read.csv(file.path(data,name),check.names=FALSE)
  colors<-c("#0072B2","#D55E00","#009E73","#555555")
  render<-function(name,draw,width=5.5,height=3.1) {
    for(type in c("pdf","png")) {
      file<-file.path(out,paste0(name,".",type))
      if(type=="pdf")grDevices::pdf(file,width=width,height=height,
                                    pointsize=9,useDingbats=FALSE)
      else grDevices::png(file,width=width*300,height=height*300,res=300,pointsize=9)
      tryCatch({graphics::par(mar=c(3.2,3.5,.8,.7),mgp=c(2,.55,0),tcl=-.2,bty="l");draw()},
               finally=grDevices::dev.off())
    }
  }
  r<-read("ridge-draws.csv")
  ridge<-structure(list(beta=array(as.matrix(r[,c("x1","x2")]),c(1000,2,2)),
    intercept=matrix(r$intercept,1000,2),sigma2=matrix(r$sigma2,1000,2),
    prior="ridge",feature.names=c("x1","x2")),class="vsusp_fit")
  B<-rbind(c(.01,.02,.03,.3,.31,10,11),c(-.01,-.02,-.03,-.3,-.31,-10,-11))
  tuning<-S2MSensitivity(B,c(.1,.5),clustering="exact")
  render("tuning",function() {
    values<-plot(tuning,pch=19,col=colors[1],lwd=1.3)
    stopifnot(identical(values$count,c(4,2)))
  })
  render("coefficients",function() {
    values<-PlotCoefficients(ridge)
    stopifnot(isTRUE(all.equal(values$mean,unname(colMeans(r[,c("x1","x2")])))))
  })
  spike.fixed<-rbind(c(1,0,0),c(2,0,0),c(0,-1,0),c(0,0,0))
  colnames(spike.fixed)<-c("x1","x2","x3")
  inclusion<-SpikeSlabSelection(spike.fixed)
  render("inclusion",function() {
    values<-PlotInclusion(inclusion,border=NA)
    stopifnot(isTRUE(all.equal(unname(values),c(.5,.25,0))))
  })
  hs<-read("horseshoe-draws.csv");ss<-read("spike-draws.csv")
  hs.array<-array(hs$beta,c(3000,4,1),list(NULL,NULL,"x1"))
  trace<-hs.array[1:300,,,drop=FALSE]
  render("trace",function() {
    PlotTrace(trace,col=colors,lty=1:4,lwd=.75,ylab=expression(beta[1]))
    graphics::legend("topright",paste("Chain",1:4),col=colors,lty=1:4,bty="n",cex=.75,ncol=2)
  })
  render("autocorrelation",function() {
    values<-PlotAutocorrelation(hs.array,max.lag=50,col=colors,lty=1:4,lwd=1)
    graphics::abline(h=0,col="grey70",lty=3)
    graphics::legend("topright",paste("Chain",1:4),col=colors,lty=1:4,bty="n",cex=.75,ncol=2)
    stopifnot(all(values[1,]==1))
    utils::write.csv(values,file.path(data,"autocorrelation-values.csv"),row.names=FALSE)
  })
  target<-read("posterior-targets.csv")
  render("posterior-agreement",function() {
    graphics::par(mar=c(3.2,7.5,.8,.7))
    ratio<-target$difference/target$tolerance;rows<-rev(seq_len(nrow(target)))
    graphics::plot(NA,xlim=c(-1.1,1.1),ylim=c(.5,5.5),yaxt="n",ylab="",
      xlab="Difference / audit tolerance",bty="n")
    graphics::rect(-1,.5,1,5.5,col="#F2F2F2",border=NA)
    graphics::abline(v=c(-1,0,1),col=c("grey65","grey35","grey65"),lty=c(2,1,2))
    labels<-expression(paste("SS: ",Pr(beta[1]!=0~"|"~Y)),
      paste("SS: ",E(beta[1]~"|"~Y)),paste("SS: ",E(sigma^2~"|"~Y)),
      paste("HS: ",E(beta[1]~"|"~Y)),paste("HS: ",E(sigma^2~"|"~Y)))
    graphics::axis(2,at=rows,labels=labels,las=1,tick=FALSE)
    graphics::points(ratio,rows,pch=ifelse(target$method=="horseshoe",17,16),
      col=ifelse(target$method=="horseshoe",colors[1],colors[2]),cex=1)
    stopifnot(all(abs(ratio)<=1))
  },height=3.5)
  render("joint-models",function() {
    graphics::par(mfrow=c(1,2),mar=c(3.3,3.5,1.4,.6),mgp=c(2,.5,0))
    p<-target$estimate[target$quantity=="inclusion"]
    exact<-target$quadrature[target$quantity=="inclusion"]
    error<-1.96*target$MCSE[target$quantity=="inclusion"]
    graphics::plot(NA,xlim=c(.6,2.4),ylim=c(0,1),xaxt="n",xlab="",ylab="Posterior probability")
    graphics::axis(1,at=1:2,labels=expression(beta[1]==0,beta[1]!=0))
    graphics::points((1:2)-.06,c(1-exact,exact),pch=18,col="grey30")
    graphics::segments((1:2)+.06,c(1-p,p)-error,(1:2)+.06,c(1-p,p)+error,col=colors[1])
    graphics::points((1:2)+.06,c(1-p,p),pch=16,col=colors[1])
    graphics::legend("topright",c("Integration","Sampling"),pch=c(18,16),col=c("grey30",colors[1]),bty="n",cex=.7)
    graphics::title(main="a  Inclusion models",adj=0,font.main=1,cex.main=.9)
    PlotTrace(array(ss$beta,c(3000,4,1))[1:300,,,drop=FALSE],col=colors,lty=1:4,lwd=.65,
      ylab=expression(beta[1]))
    graphics::title(main="b  Coefficient traces",adj=0,font.main=1,cex.main=.9)
    stopifnot(abs(mean(ss$beta!=0)-p)<1e-12)
  },width=6.8,height=3.2)
  source(file.path(root,"inst","validation","plot-selection-pilot.R"),local=TRUE)
  plot_selection_pilot(file.path(root,"inst","validation","summary.csv"),
                       file.path(out,"selection-pilot.pdf"))
  # PDF is the vector source; PNG is the HTML alternative.
  command<-Sys.which("pdftoppm")
  if(!nzchar(command))stop("pdftoppm is required to export the pilot PNG.")
  status<-system2(command,c("-singlefile","-r","300","-png",
                 shQuote(file.path(out,"selection-pilot.pdf")),
                 shQuote(file.path(out,"selection-pilot"))),stdout=FALSE,stderr=FALSE)
  if(status!=0L)stop("Pilot PNG conversion failed.")
  cat("Built eight figure pairs; data and generation source are retained.\n")
}
main()
