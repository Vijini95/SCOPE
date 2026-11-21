# =============================================================================
# FIGURES -- BASE R PDF and PNG output; called automatically by experiments.
# Produces per-condition and aggregate graphics, including Figure 2 sensitivity
# across N and root-SNR. No external graphics libraries needed.
# =============================================================================

scope_e1_plots <- function(tab,plotDir) {
  scope_make_dir(plotDir)
  signal <- unique(tab$Signal);snr <- unique(tab$SNR)
  tmp <- tab[grepl("_r[0-3]$",tab$Method),]
  Fs <- unique(sub("_r[0-3]$","",tmp$Method))
  amse <- sapply(0:3,function(r)
    vapply(Fs,function(F) tmp$AMSE[match(paste0(F,"_r",r),tmp$Method)],numeric(1)))
  gains <- 100*(amse[,1L]-amse[,2:4,drop=FALSE])/amse[,1L]
  scope_save_figure(file.path(plotDir,sprintf("E1_AMSE_%s_SNR%d",signal,snr)),
    function() {
      op <- par(mar=c(7,5,3,1));on.exit(par(op))
      barplot(t(amse),beside=TRUE,names.arg=Fs,las=2,
              col=gray(seq(.25,.75,length.out=4)),
              ylab="Average MSE (per observation)",
              main=sprintf("Fixed SCOPE / N-SCOPE | %s | root-SNR=%d",signal,snr))
      legend("topright",paste0("r=",0:3),fill=gray(seq(.25,.75,length.out=4)),bty="n")
    },width=11,height=7)
  scope_save_figure(file.path(plotDir,sprintf("E1_GainPercent_%s_SNR%d",signal,snr)),
    function() {
      op <- par(mar=c(7,5,3,1));on.exit(par(op))
      bp <- barplot(t(gains),beside=TRUE,names.arg=Fs,las=2,
               col=gray(seq(.30,.70,length.out=3)),
               ylab="AMSE reduction relative to same-generator r=0 (%)",
               main=sprintf("Neighborhood effect | %s | root-SNR=%d",signal,snr))
      abline(h=0,lty=2);legend("topright",paste0("r=",1:3),
                    fill=gray(seq(.30,.70,length.out=3)),bty="n")
    },width=11,height=7)
}

scope_e2_plots <- function(tab,selection,plotDir) {
  scope_make_dir(plotDir)
  sig <- unique(tab$Signal);snr <- unique(tab$SNR)
  Fs <- unique(selection$Generator)
  vals <- t(vapply(Fs,function(F)
    tab$AMSE[match(paste0(F,c("_FixedSCOPE","_FixedNSCOPE",
                  "_SCOPE_SURE","_NSCOPE_SURE")),tab$Method)],numeric(4)))
  scope_save_figure(file.path(plotDir,sprintf("E2_AMSE_%s_SNR%d",sig,snr)),
    function() {
      op <- par(mar=c(7,5,3,1));on.exit(par(op))
      barplot(t(vals),beside=TRUE,names.arg=Fs,las=2,
         col=gray(seq(.2,.75,length.out=4)),
         main=sprintf("Experiment II | %s | root-SNR=%d",sig,snr),ylab="AMSE")
      legend("topright",c("Fixed r=0","Fixed r=1","SCOPE SURE","N-SCOPE SURE"),
         fill=gray(seq(.2,.75,length.out=4)),bty="n",cex=.8)
    },width=11,height=7)
  # Distribution of selected radius for each generating CDF.
  selection <- subset(selection,Method=="NSCOPE_SURE")
  fr <- sapply(Fs,function(F) {
    as.numeric(table(factor(selection$r[selection$Generator==F],levels=0:3))) /
      sum(selection$Generator==F)
  })
  scope_save_figure(file.path(plotDir,sprintf("E2_RadiusSelection_%s_SNR%d",sig,snr)),
    function() {
      op <- par(mar=c(7,5,3,1));on.exit(par(op))
      barplot(fr,beside=FALSE,names.arg=Fs,las=2,
              col=gray(seq(.20,.80,length.out=4)),ylim=c(0,1),
              main=sprintf("SURE-selected neighborhood radius | %s | SNR=%d",sig,snr),
              ylab="Selection frequency")
      legend("topright",paste0("r=",0:3),
             fill=gray(seq(.20,.80,length.out=4)),bty="n")
    },width=11,height=7)
}

scope_e2_aggregate_plots <- function(rank,selection,plotDir) {
  scope_save_figure(file.path(plotDir,"E2_AggregateRanking"),function() {
     ord <- order(rank$Ratio,decreasing=TRUE)
     op <- par(mar=c(5,15,3,1));on.exit(par(op))
     barplot(rank$Ratio[ord],names.arg=rank$Method[ord],horiz=TRUE,las=1,
             col="gray60",main="Experiment II: mean AMSE / universal hard AMSE",
             xlab="Mean relative AMSE",cex.names=.65)
   },width=12,height=9)
  selection <- subset(selection,Method=="NSCOPE_SURE" & Generator=="Normal")
  if(nrow(selection)==0) return(invisible(NULL))
  sigs <- intersect(scope_conditions,unique(selection$Signal))
  freq <- sapply(sigs,function(s)
    as.numeric(table(factor(selection$r[selection$Signal==s],levels=0:3))) /
      sum(selection$Signal==s))
  scope_save_figure(file.path(plotDir,"E2_RadiusSelection_BySignal"),function() {
    op <- par(mar=c(6,5,3,1));on.exit(par(op))
    barplot(freq,beside=FALSE,names.arg=sigs,las=2,
       col=gray(seq(.20,.80,length.out=4)),ylim=c(0,1),
       ylab="Fraction of replications",main="Normal N-SCOPE SURE: selected radius")
    legend("topright",paste0("r=",0:3),fill=gray(seq(.20,.80,length.out=4)),bty="n")
   },width=10,height=7)
}

scope_e3_plots <- function(tab,mse,plotDir) {
  scope_make_dir(plotDir)
  signal <- unique(tab$Signal);snr <- unique(tab$SNR)
  ord <- order(tab$AMSE,decreasing=TRUE)
  scope_save_figure(file.path(plotDir,sprintf("E3_AMSE_%s_SNR%d",signal,snr)),
    function() {
      op <- par(mar=c(5,15,3,1));on.exit(par(op))
      barplot(tab$AMSE[ord],names.arg=tab$Method[ord],horiz=TRUE,
        col="gray65",las=1,cex.names=.6,
        xlab="AMSE",main=sprintf("All 26 methods | %s | root-SNR=%d",signal,snr))
    },width=12,height=11)
  # Per-replication boxplots shown on the log scale because some methods
  # may be substantially worse than others.
  scope_save_figure(file.path(plotDir,sprintf("E3_MSEBox_%s_SNR%d",signal,snr)),
    function() {
      op <- par(mar=c(13,5,3,1));on.exit(par(op))
      boxplot(as.data.frame(mse),log="y",las=2,outline=FALSE,
              cex.axis=.55,col="gray85",ylab="Replication MSE (log scale)",
              main=sprintf("Paired MSE distributions | %s | root-SNR=%d",signal,snr))
    },width=13,height=9)
}

scope_e3_aggregate_plot <- function(rank,plotDir) {
  scope_save_figure(file.path(plotDir,"E3_AggregateRanking"),function() {
    ord <- order(rank$Ratio,decreasing=TRUE)
    op <- par(mar=c(5,16,3,1));on.exit(par(op))
    barplot(rank$Ratio[ord],names.arg=rank$Method[ord],horiz=TRUE,
       las=1,col="gray60",cex.names=.65,
       xlab="Mean AMSE / universal hard AMSE",
       main="Experiment III: 26-method descriptive ranking")
  },width=13,height=11)
}

scope_plot_N_SNR_sensitivity <- function(root,Ns,snrs,M) {
  if(length(Ns)<2L || length(snrs)<2L) return(invisible(NULL))
  gains <- array(NA_real_,c(length(scope_conditions),length(Ns),length(snrs)))
  for(i in seq_along(Ns)) {
    file <- file.path(root,sprintf("Exp1_N%d_M%d",Ns[i],M),
                      "Experiment1_AMSE_all.csv")
    if(!file.exists(file)) return(invisible(NULL))
    tab <- read.csv(file,stringsAsFactors=FALSE,check.names=FALSE)
    for(s in seq_along(scope_conditions)) for(q in seq_along(snrs)) {
      tt <- subset(tab,Signal==scope_conditions[s] & SNR==snrs[q])
      a <- tt$AMSE[tt$Method=="Normal_r0"]
      b <- tt$AMSE[tt$Method=="Normal_r1"]
      if(length(a)!=1 || length(b)!=1) return(invisible(NULL))
      gains[s,i,q] <- 100*(a-b)/a
    }
  }
  scope_save_figure(file.path(scope_make_dir(file.path(root,"plots")),
                               sprintf("E1_Figure2_Sensitivity_N_SNR_M%d",M)),function() {
    op <- par(mfrow=c(2,2),mar=c(4,5,3,1));on.exit(par(op))
    ltys <- seq_along(snrs); pchs <- c(16,17,15,18)[seq_along(snrs)]
    for(s in seq_along(scope_conditions)) {
      y <- gains[s,,]; y <- matrix(y,nrow=length(Ns))
      yl <- range(c(0,y),finite=TRUE);padding <- max(diff(yl)*.1,2)
      matplot(Ns,y,type="b",lty=ltys,pch=pchs,col=seq_along(snrs),lwd=2,
              ylim=yl+c(-padding,padding),xlab="Sample size N",
              ylab="AMSE reduction (%)",main=scope_conditions[s],xaxt="n")
      axis(1,at=Ns,labels=Ns);abline(h=0,lty=3,col="gray50");grid()
      if(s==1) legend("topright",paste0("Root SNR=",snrs),lty=ltys,pch=pchs,
                      col=seq_along(snrs),bty="n",cex=.8)
    }
    mtext("Fixed Normal N-SCOPE (r=1) versus SCOPE (r=0)",outer=TRUE,
           line=-1,cex=1.0)
  },width=12,height=9)
}
