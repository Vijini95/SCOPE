# =============================================================================
# Experiments I and II. Closely follow the MATLAB mathematical definitions.
# Both use Gaussian N(0,1) noise, N=512/1024/2048, root-SNR=3/5/7, M=1000.
# Results are regenerated in R, not imported from MATLAB checkpoints.
# =============================================================================

scope_e1_names <- function(Fs) c(unlist(lapply(Fs,function(F) paste0(F,"_r",0:3))),
                                   "UniversalHard","UniversalSoft")
scope_e2_names <- function(Fs) c(unlist(lapply(Fs,function(F)
  paste0(F,c("_FixedSCOPE","_FixedNSCOPE","_SCOPE_SURE","_NSCOPE_SURE")))),
  "UniversalHard","UniversalSoft")

scope_condition_setup <- function(signal,N,snr) {
  clean <- scope_root_snr_scale(scope_signal(signal,N),snr)
  filter <- scope_filters(scope_filter_for_signal(signal))
  list(clean=clean,filters=filter,truth=scope_dwt(clean,filter))
}

scope_e1_condition <- function(signal,si,snr,qi,N,M,Fs,cal,outDir,saveEvery=25L) {
  names <- scope_e1_names(Fs); setup <- scope_condition_setup(signal,N,snr)
  path <- file.path(outDir,sprintf("checkpoint_%s_SNR%d.rds",signal,snr))
  set.seed(scope_seed(N,"Exp1",si,snr))
  signature <- list(version="R_E1_20261029_1",N=N,M=M,signal=signal,
                    snr=snr,Fs=Fs,eta=cal$eta)
  init <- list(completed=0L,mse=matrix(NA_real_,M,length(names),dimnames=list(NULL,names)),
               sigma=rep(NA_real_,M))
  state <- scope_run_state(path,init,signature)
  if (state$completed<M) for (it in seq.int(state$completed+1L,M)) {
    y <- setup$clean+rnorm(N)
    d <- scope_dwt(y,setup$filters); sigma <- scope_sigma(d)
    blocks <- scope_neighborhood(d,rGrid=0:3)
    vals <- numeric(length(names)); id <- 0L
    for (F in Fs) for (r in 0:3) {
      id <- id+1L
      est <- scope_apply(d,blocks,sigma,cal,F,5L,r)
      vals[id] <- scope_amse(est,setup$truth)
    }
    vals[length(names)-1L] <- scope_amse(scope_universal(d,sigma,kind="hard"),setup$truth)
    vals[length(names)] <- scope_amse(scope_universal(d,sigma,kind="soft"),setup$truth)
    state$mse[it,] <- vals;state$sigma[it] <- sigma; state$completed <- it
    if (it%%saveEvery==0L || it==M) {
      scope_write_state(state,path)
      cat(sprintf("Exp1 N=%d %-10s SNR=%d: %d/%d\n",N,signal,snr,it,M))
    }
  }
  stopifnot(all(is.finite(state$mse)))
  write.csv(data.frame(Replication=seq_len(M),SigmaHat=state$sigma,state$mse,
                       check.names=FALSE),
            file.path(outDir,sprintf("MSE_replications_%s_SNR%d.csv",signal,snr)),row.names=FALSE)
  mu <- colMeans(state$mse); sd <- apply(state$mse,2,stats::sd)
  res <- data.frame(N=N,Signal=signal,SNR=snr,Method=names,M=M,AMSE=mu,
                    SD_MSE=sd,SE_AMSE=sd/sqrt(M),row.names=NULL)
  for (fi in seq_along(Fs)) for (r in 1:3) {
    i0 <- (fi-1L)*4L+1L; ir <- i0+r
    res$GainPct_vs_sameF_r0[ir] <- 100*(mu[i0]-mu[ir])/mu[i0]
  }
  scope_e1_plots(res,file.path(outDir,"plots"))
  list(summary=res,state=state)
}

scope_run_e1 <- function(N,snrs,M,Fs,root,smoke=FALSE) {
  out <- scope_make_dir(file.path(root,sprintf("Exp1_N%d_M%d",N,M)))
  cal <- scope_get_calibration(file.path(out,"calibration_E1.rds"),N,Fs,5L,0:3)
  calTab <- expand.grid(Generator=Fs,Exponent=5L,Radius=0:3,
                        Level=cal$j0:(as.integer(log2(N))-1L))
  calTab$Eta_j <- cal$targets[match(calTab$Level,cal$j0:(as.integer(log2(N))-1L))]
  calTab$a_NullRisk <- mapply(function(F,k,r,j) cal$a[match(F,Fs),1L,r+1L,j-2L],
                             calTab$Generator,calTab$Exponent,calTab$Radius,calTab$Level)
  write.csv(calTab,file.path(out,"Experiment1_Constants.csv"),row.names=FALSE)
  signals <- if(smoke) c("Doppler","Blocks") else scope_conditions
  res <- list();pair <- list(); i <- 0L
  for (si in seq_along(signals)) for (qi in seq_along(snrs)) {
    sig <- signals[si];snr <- snrs[qi]
    a <- scope_e1_condition(sig,match(sig,scope_conditions),snr,qi,N,M,Fs,cal,out)
    res[[length(res)+1L]] <- a$summary
    for (F in Fs) for (r in 1:3) {
      n0 <- paste0(F,"_r0"); nr <- paste0(F,"_r",r)
      delta <- a$state$mse[,n0]-a$state$mse[,nr]
      pair[[length(pair)+1L]] <- data.frame(N=N,Signal=sig,SNR=snr,Generator=F,
        Radius=r,PairedDeltaAMSE=mean(delta),PairedSE_Delta=sd(delta)/sqrt(M),
        CI95_Lower=mean(delta)-1.96*sd(delta)/sqrt(M),
        CI95_Upper=mean(delta)+1.96*sd(delta)/sqrt(M),
        GainPercent=100*mean(delta)/mean(a$state$mse[,n0]))
    }
  }
  tab <- do.call(rbind,res);tab$RMSE <- sqrt(tab$AMSE)
  write.csv(tab,file.path(out,"Experiment1_AMSE_all.csv"),row.names=FALSE)
  write.csv(do.call(rbind,pair),file.path(out,"Experiment1_PairedGains.csv"),row.names=FALSE)
  ref <- subset(tab,Method=="Normal_r0",c("Signal","SNR","AMSE"))
  rank <- merge(tab,ref,by=c("Signal","SNR"),suffixes=c("","_ref"))
  rank$Relative <- rank$AMSE/rank$AMSE_ref
  rank <- aggregate(Relative~Method,rank,mean)
  rank <- rank[order(rank$Relative),]
  write.csv(rank,file.path(out,"Experiment1_AggregateRanking.csv"),row.names=FALSE)
  cat("Experiment I complete:",out,"\n")
  invisible(tab)
}

scope_e2_condition <- function(signal,si,snr,qi,N,M,Fs,cal,outDir,saveEvery=10L) {
  names <- scope_e2_names(Fs)
  cand <- scope_candidates();setup <- scope_condition_setup(signal,N,snr)
  path <- file.path(outDir,sprintf("checkpoint_%s_SNR%d.rds",signal,snr))
  set.seed(scope_seed(N,"Exp2",si,snr))
  signature <- list(version="R_E2_20261029_1",N=N,M=M,signal=signal,
                    snr=snr,Fs=Fs,eta=cal$eta,candidates=cand)
  init <- list(completed=0L,mse=matrix(NA_real_,M,length(names),dimnames=list(NULL,names)),
               sigma=rep(NA_real_,M),
               selected=matrix(NA_integer_,M,2L*length(Fs)),
               selectedK=matrix(NA_integer_,M,2L*length(Fs)),
               selectedR=matrix(NA_integer_,M,2L*length(Fs)),
               selectedB=matrix(NA_real_,M,2L*length(Fs)))
  state <- scope_run_state(path,init,signature)
  if (state$completed<M) for (it in seq.int(state$completed+1L,M)) {
    y <- setup$clean+rnorm(N);d <- scope_dwt(y,setup$filters)
    sigma <- scope_sigma(d); blocks <- scope_neighborhood(d,rGrid=0:3)
    vals <- numeric(length(names));pos <- 0L
    for (fi in seq_along(Fs)) {
      F <- Fs[fi];scores <- vapply(seq_len(nrow(cand)),function(c)
        scope_sure(blocks,sigma,cal,F,cand$k[c],cand$r[c],cand$b[c]),numeric(1))
      iS <- which(cand$r==0); i0 <- iS[which.min(scores[iS])]; i1 <- which.min(scores)
      scope_assert(scores[i1]<=scores[i0]+1e-7*max(1,abs(scores[i0])),
                   "SURE nesting failed")
      sel <- list(list(k=5L,r=0L,b=1),list(k=5L,r=1L,b=1),
                  as.list(cand[i0,,drop=FALSE]),as.list(cand[i1,,drop=FALSE]))
      for (m in seq_len(4L)) {
        pos <- pos+1L;q <- sel[[m]]
        val <- scope_apply(d,blocks,sigma,cal,F,q$k,q$r,q$b)
        vals[pos] <- scope_amse(val,setup$truth)
      }
      state$selected[it,(2L*fi-1L):(2L*fi)] <- c(i0,i1)
      state$selectedK[it,(2L*fi-1L):(2L*fi)] <- c(cand$k[i0],cand$k[i1])
      state$selectedR[it,(2L*fi-1L):(2L*fi)] <- c(cand$r[i0],cand$r[i1])
      state$selectedB[it,(2L*fi-1L):(2L*fi)] <- c(cand$b[i0],cand$b[i1])
    }
    vals[length(names)-1L] <- scope_amse(scope_universal(d,sigma,kind="hard"),setup$truth)
    vals[length(names)] <- scope_amse(scope_universal(d,sigma,kind="soft"),setup$truth)
    state$mse[it,] <- vals;state$sigma[it] <- sigma;state$completed <- it
    if (it%%saveEvery==0L || it==M) {
      scope_write_state(state,path)
      cat(sprintf("Exp2 N=%d %-10s SNR=%d: %d/%d\n",N,signal,snr,it,M))
    }
  }
  write.csv(data.frame(Replication=seq_len(M),SigmaHat=state$sigma,state$mse,
                       check.names=FALSE),
           file.path(outDir,sprintf("MSE_replications_%s_SNR%d.csv",signal,snr)),row.names=FALSE)
  sel <- list(); for (fi in seq_along(Fs)) for (method in 1:2) {
    idx <- 2L*fi-2L+method
    sel[[length(sel)+1L]] <- data.frame(Replication=seq_len(M),Generator=Fs[fi],
      Method=if(method==1) "SCOPE_SURE" else "NSCOPE_SURE",
      k=state$selectedK[,idx],r=state$selectedR[,idx],b=state$selectedB[,idx])
  }
  write.csv(do.call(rbind,sel),file.path(outDir,sprintf("SelectedParams_%s_SNR%d.csv",signal,snr)),row.names=FALSE)
  mu <- colMeans(state$mse);sd0 <- apply(state$mse,2,sd)
  tab <- data.frame(N=N,Signal=signal,SNR=snr,Method=names,M=M,AMSE=mu,
                    SD_MSE=sd0,SE_AMSE=sd0/sqrt(M),RMSE=sqrt(mu),row.names=NULL)
  scope_e2_plots(tab,do.call(rbind,sel),file.path(outDir,"plots"))
  list(summary=tab,selection=do.call(rbind,sel),state=state)
}

scope_run_e2 <- function(N,snrs,M,Fs,root,smoke=FALSE) {
  out <- scope_make_dir(file.path(root,sprintf("Exp2_N%d_M%d",N,M)))
  cal <- scope_get_calibration(file.path(out,"calibration_E2.rds"),N,Fs,
                               c(3L,5L,7L),0:3)
  signals <- if(smoke) c("Doppler","Blocks") else scope_conditions
  all <- list(); selections <- list();pairs <- list()
  for (signal in signals) for (qi in seq_along(snrs)) {
    snr <- snrs[qi]
    z <- scope_e2_condition(signal,match(signal,scope_conditions),snr,qi,N,M,Fs,cal,out)
    all[[length(all)+1L]] <- z$summary
    ss <- z$selection;ss$N <- N;ss$Signal <- signal;ss$SNR <- snr
    selections[[length(selections)+1L]] <- ss
    for (F in Fs) {
      nm0 <- paste0(F,"_SCOPE_SURE");nm1 <- paste0(F,"_NSCOPE_SURE")
      delta <- z$state$mse[,nm0]-z$state$mse[,nm1]
      pairs[[length(pairs)+1L]] <- data.frame(N=N,Signal=signal,SNR=snr,
        Generator=F,PairedDeltaAMSE=mean(delta),PairedSE=sd(delta)/sqrt(M),
        CI95_Lower=mean(delta)-1.96*sd(delta)/sqrt(M),
        CI95_Upper=mean(delta)+1.96*sd(delta)/sqrt(M))
    }
  }
  tab <- do.call(rbind,all)
  write.csv(tab,file.path(out,"Experiment2_AMSE_all.csv"),row.names=FALSE)
  write.csv(do.call(rbind,pairs),file.path(out,"Experiment2_PairedDifferences.csv"),row.names=FALSE)
  ss <- do.call(rbind,selections)
  write.csv(ss,file.path(out,"Experiment2_SelectedTriples.csv"),row.names=FALSE)
  freq <- aggregate(rep(1,nrow(ss)),by=ss[c("N","Signal","SNR","Generator","Method","k","r","b")],
                    FUN=length); names(freq)[ncol(freq)] <- "Count"
  freq$Frequency <- freq$Count/M
  write.csv(freq,file.path(out,"Experiment2_SelectionFrequencies.csv"),row.names=FALSE)
  ref <- subset(tab,Method=="UniversalHard",c("Signal","SNR","AMSE"))
  rank <- merge(tab,ref,by=c("Signal","SNR"),suffixes=c("","_Hard"))
  rank$Ratio <- rank$AMSE/rank$AMSE_Hard
  rank <- aggregate(Ratio~Method,rank,mean)
  write.csv(rank[order(rank$Ratio),],file.path(out,"Experiment2_AggregateRanking.csv"),row.names=FALSE)
  scope_e2_aggregate_plots(rank,ss,file.path(out,"plots"))
  cat("Experiment II complete:",out,"\n")
  invisible(tab)
}
