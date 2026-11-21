# ============================================================================
# SCOPE_Exp123_ALL_IN_ONE.R
# ANNOTATED single-file, base-R translation of MATLAB SCOPE_Exp123_FINAL.m
# Experiments 1, 2, 3 + automatic PDF and PNG figures. No CRAN dependencies.
# After source(...): SCOPE_Exp123_FINAL("smoke","all")
# All output is placed under getwd()/SCOPE_Exp123_R_results.
# For a modular distribution, see SCOPE_Exp123_FINAL.R and companion modules.
# ============================================================================
.scope_dir <- getwd()



# =============================================================================
# SCOPE / N-SCOPE: shared mathematics (BASE R ONLY)
# Translation of SCOPE_Exp123_FINAL.m, 2026-10-29.
# No wavelet, optimization, data.frame or plotting add-on packages required.
# This port deliberately uses its own R random numbers. R/MATLAB Monte Carlo
# realizations differ; statistical replication, not bitwise numerical identity.
# =============================================================================

scope_assert <- function(ok, msg) if (!isTRUE(all(ok))) stop(msg, call. = FALSE)

scope_filters <- function(name) {
  h <- switch(tolower(name),
    haar = c(1, 1)/sqrt(2),
    db3 = c(.332670552950082616, .806891509311092576,
            .459877502118491570, -.135011020010254589,
            -.085441273882026662, .035226291885709537),
    sym4 = c(-.07576571478950221, -.029635527646002493,
              .497618667632775, .8037387518051321,
              .29785779560530606, -.09921954357663353,
             -.012603967262031304, .032223100604051466),
    stop("Unsupported filter"))
  h <- h/sqrt(sum(h^2))
  # Same QMF convention as MATLAB: fliplr(h .* (-1).^(1:length(h))).
  list(h=h, g=rev(h * (-1)^seq_along(h)))
}

scope_dwt_step <- function(x, filters, shift = 2L) {
  n <- length(x); half <- n %/% 2L
  a <- numeric(half); d <- numeric(half)
  for (r in seq_len(half)) {
    ind <- ((seq_along(filters$h)+2L*r-shift-1L) %% n)+1L
    a[r] <- sum(filters$h * x[ind]); d[r] <- sum(filters$g * x[ind])
  }
  list(a=a,d=d)
}

scope_idwt_step <- function(a,d,filters,shift=2L) {
  half <- length(a); n <- 2L*half
  x <- numeric(n)
  for (r in seq_len(half)) {
    ind <- ((seq_along(filters$h)+2L*r-shift-1L) %% n)+1L
    x[ind] <- x[ind] + filters$h*a[r] + filters$g*d[r]
  }
  x
}

scope_dwt <- function(x, filters, j0=3L) {
  # Coefficient order [2^j0 scaling, j0 detail, ..., J-1 detail].
  n <- length(x); J <- as.integer(round(log2(n)))
  scope_assert(2^J==n && j0>=1L && j0<J, "Invalid dyadic DWT length")
  z <- numeric(n); a <- as.numeric(x)
  for (j in seq(J-1L,j0,by=-1L)) {
    s <- scope_dwt_step(a,filters)
    z[(2^j+1L):2^(j+1L)] <- s$d
    a <- s$a
  }
  z[seq_len(2^j0)] <- a; z
}

scope_idwt <- function(z, filters, j0=3L) {
  n <- length(z); J <- as.integer(round(log2(n)))
  a <- z[seq_len(2^j0)]
  for (j in j0:(J-1L)) {
    a <- scope_idwt_step(a,z[(2^j+1L):2^(j+1L)],filters)
  }
  a
}

scope_filter_for_signal <- function(signal) switch(tolower(signal),
  doppler="sym4",heavisine="sym4",bumps="db3",blocks="haar")

scope_signal <- function(signal,N) {
  # Exact Donoho-Johnstone variants and t=i/N of MATLAB production code.
  t <- seq_len(N)/N
  pos <- c(.10,.13,.15,.23,.25,.40,.44,.65,.76,.78,.81)
  htB <- c(4,5,3,4,5,4.2,2.1,4.3,3.1,5.1,4.2)
  htK <- c(4,-5,3,-4,5,-4.2,2.1,4.3,-3.1,2.1,-4.2)
  w <- c(.005,.005,.006,.010,.010,.030,.010,.010,.005,.008,.005)
  s <- switch(tolower(signal),
    doppler=sqrt(t*(1-t))*sin(2*pi*1.05/(t+.05)),
    heavisine=4*sin(4*pi*t)-sign(t-.30)-sign(.72-t),
    bumps=vapply(t,function(tt) sum(htB/(1+abs((tt-pos)/w))^4),numeric(1)),
    blocks=vapply(t,function(tt) sum(htK*(tt>=pos)),numeric(1)),
    stop("Unknown Donoho-Johnstone signal"))
  as.numeric(s)
}

scope_root_snr_scale <- function(f0,snr) {
  # Uses population standard deviation, MATLAB std(f,1).
  f0 * (snr/sqrt(mean((f0-mean(f0))^2)))
}

scope_fstar <- function(F,u) {
  u <- pmax(0, u)
  pmin(1,pmax(0,switch(tolower(F),
    normal=erf_r(u/sqrt(2)), logistic=tanh(u),
    laplace=-expm1(-u), sech=(2/pi)*atan(sinh(pi*u/2)),
    uniform=pmin(u,1), cauchy=(2/pi)*atan(u),
    stop("Unknown centered CDF"))))
}

erf_r <- function(x) 2*pnorm(x*sqrt(2))-1

scope_fstar_prime <- function(F,u) {
  u <- pmax(u,0)
  switch(tolower(F),normal=sqrt(2/pi)*exp(-u^2/2),
    logistic=1-tanh(u)^2, laplace=exp(-u),
    sech=1/cosh(pi*u/2), uniform=as.numeric(u<1),
    cauchy=2/(pi*(1+u^2)), stop("Unknown centered CDF"))
}

scope_eta <- function(N) {
  # Normal F*, k=5, half attenuation at sqrt(2 log N).
  k <- 5L; a <- qnorm((1+2^(-1/k))/2)/sqrt(2*log(N))
  integrate(function(u) u*scope_fstar("Normal",a*sqrt(u))^(2*k)*dchisq(u,df=1),
            0,Inf,rel.tol=1e-9,subdivisions=500L)$value
}

scope_null_risk <- function(a,F,r,k) {
  m <- 2L*r+1L
  if (a<=0) return(0)
  integrate(function(u) (u/m)*scope_fstar(F,a*sqrt(u/m))^(2*k)*dchisq(u,m),
            0,Inf,rel.tol=1e-7,subdivisions=300L)$value
}

scope_calibrate <- function(N,Fs,kGrid,rGrid,j0=3L) {
  # a[F,k,r,j] solves exact null-risk equation. Cached by the experiment.
  J <- as.integer(round(log2(N))); eta <- scope_eta(N)
  targets <- eta * 2^(J-1L-(j0:(J-1L)))
  scope_assert(max(targets)<1, "Invalid null-risk schedule (coarse target >=1)")
  a <- array(NA_real_,dim=c(length(Fs),length(kGrid),length(rGrid),J-j0),
             dimnames=list(Fs,as.character(kGrid),as.character(rGrid),
                           as.character(j0:(J-1L))))
  for (fi in seq_along(Fs)) for (ki in seq_along(kGrid))
    for (ri in seq_along(rGrid)) for (jj in seq_along(targets)) {
      fn <- function(v) scope_null_risk(v,Fs[fi],rGrid[ri],kGrid[ki])-targets[jj]
      upper <- 1
      while (fn(upper)<0) {upper <- 2*upper;scope_assert(upper<1e6,"Root not bracketed")}
      a[fi,ki,ri,jj] <- uniroot(fn,c(0,upper),tol=1e-9)$root
    }
  list(a=a,eta=eta,targets=targets,N=N,j0=j0,Fs=Fs,kGrid=kGrid,rGrid=rGrid)
}

scope_get_calibration <- function(path,N,Fs,kGrid,rGrid,j0=3L) {
  meta <- list(N=N,Fs=Fs,kGrid=kGrid,rGrid=rGrid,j0=j0)
  if (file.exists(path)) {
    x <- readRDS(path)
    scope_assert(identical(x$meta,meta), "Calibration cache is incompatible")
    return(x$cal)
  }
  cat("Calibrating Gaussian null risk for N=",N,"...\n",sep="")
  cal <- scope_calibrate(N,Fs,kGrid,rGrid,j0)
  saveRDS(list(meta=meta,cal=cal),path)
  cal
}

scope_circshift <- function(x,lag) { # positive lag = MATLAB circshift right
  n <- length(x); x[((seq_len(n)-lag-1L) %% n)+1L]
}

scope_neighborhood <- function(d,j0=3L,rGrid=0:3) {
  N <- length(d); J <- as.integer(round(log2(N)))
  blocks <- vector("list",J-j0)
  for (j in j0:(J-1L)) {
    x <- d[(2^j+1L):2^(j+1L)]; sj <- vector("list",length(rGrid))
    for (ri in seq_along(rGrid)) {
      r <- rGrid[ri]; eng <- x^2
      if (r>0L) for (lag in seq_len(r))
        eng <- eng+scope_circshift(x^2,lag)+scope_circshift(x^2,-lag)
      sj[[ri]] <- sqrt(eng/(2*r+1L))
    }
    blocks[[j-j0+1L]] <- list(j=j,x=x,S=sj)
  }
  blocks
}

scope_apply <- function(d,blocks,sigma,cal,F,k,r,b=1,j0=3L) {
  out <- d
  fi <- match(F,cal$Fs); ki <- match(k,cal$kGrid)
  ri <- match(r,cal$rGrid)
  scope_assert(!anyNA(c(fi,ki,ri)), "Candidate not in calibration grid")
  for (jj in seq_along(blocks)) {
    blk <- blocks[[jj]]; ind <- (2^blk$j+1L):2^(blk$j+1L)
    lam <- b*cal$a[fi,ki,ri,jj]/sigma
    out[ind] <- blk$x*scope_fstar(F,lam*blk$S[[ri]])^k
  }
  out
}

scope_sure <- function(blocks,sigma,cal,F,k,r,b=1) {
  # Plug-in Stein risk: sigmahat and lambda held fixed in partial derivative.
  fi <- match(F,cal$Fs); ki <- match(k,cal$kGrid);ri <- match(r,cal$rGrid)
  m <- 2*r+1; risk <- 0
  for (jj in seq_along(blocks)) {
    blk <- blocks[[jj]]; x <- blk$x; S <- blk$S[[ri]]
    lam <- b*cal$a[fi,ki,ri,jj]/sigma
    u <- lam*S; g <- scope_fstar(F,u);gp <- scope_fstar_prime(F,u)
    w <- g^k
    div <- w + k*lam*x^2/(m*pmax(S,.Machine$double.xmin))*g^(k-1)*gp
    risk <- risk + sum((x*(w-1))^2+2*sigma^2*div-sigma^2)
  }
  scope_assert(is.finite(risk),"Non-finite SURE")
  risk
}

scope_candidates <- function(kGrid=c(3L,5L,7L),rGrid=0:3,
                             bGrid=c(.65,.80,1,1.2,1.5)) {
  # MATLAB traversal order: radius, then k, then multiplier.
  expand.grid(b=bGrid,k=kGrid,r=rGrid,KEEP.OUT.ATTRS=FALSE)
}

scope_select <- function(blocks,sigma,cal,F,cand,rAllowed=cal$rGrid) {
  ids <- which(cand$r %in% rAllowed)
  scores <- vapply(ids,function(i)
    scope_sure(blocks,sigma,cal,F,cand$k[i],cand$r[i],cand$b[i]),numeric(1))
  i <- ids[which.min(scores)]; as.list(cand[i,,drop=FALSE])
}

scope_universal <- function(d,sigma,j0=3L,kind="hard") {
  N <- length(d); out <- d; ind <- (2^j0+1L):N
  T <- sigma*sqrt(2*log(N)); x <- d[ind]
  out[ind] <- if (kind=="hard") x*(abs(x)>T) else sign(x)*pmax(abs(x)-T,0)
  out
}

scope_amse <- function(theta,truth) mean((theta-truth)^2)
scope_sigma <- function(d) median(abs(tail(d,length(d)%/%2L)))/.6744897501960817
scope_conditions <- c("Doppler","HeaviSine","Bumps","Blocks")
scope_generators <- c("Normal","Logistic","Laplace","Sech","Uniform","Cauchy")

scope_make_dir <- function(path) {dir.create(path,recursive=TRUE,showWarnings=FALSE);path}

scope_save_figure <- function(path,draw,width=11,height=7) {
  scope_make_dir(dirname(path))
  grDevices::pdf(paste0(path,".pdf"),width=width,height=height,onefile=FALSE)
  tryCatch(draw(),finally=grDevices::dev.off())
  grDevices::png(paste0(path,".png"),width=round(width*160),height=round(height*160),res=160)
  tryCatch(draw(),finally=grDevices::dev.off())
}

scope_seed <- function(N,experiment,signal_index,snr) {
  # Deterministic across both full and single-SNR runs.
  # R-native pseudorandom generator; not MATLAB rng('twister') replication.
  snr_index <- match(snr,c(3,5,7))
  if(is.na(snr_index)) snr_index <- as.integer(100*snr)+10L
  as.integer((if(experiment=="Exp3") 20271126 else 20261026) +
               N+100L*signal_index+snr_index)
}

scope_run_state <- function(path,init,signature) {
  if (file.exists(path)) {
    x <- readRDS(path)
    scope_assert(identical(x$signature,signature),
                 paste("Checkpoint incompatible:",path))
    assign(".Random.seed",x$rng,envir=.GlobalEnv)
    x
  } else {
    x <- init; x$signature <- signature
    x$rng <- get(".Random.seed",envir=.GlobalEnv)
    x
  }
}

scope_write_state <- function(state,path) {
  state$rng <- get(".Random.seed",envir=.GlobalEnv)
  tmp <- paste0(path,".tmp");saveRDS(state,tmp)
  ok <- file.rename(tmp,path)
  scope_assert(ok,paste("Cannot save checkpoint",path))
}


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


# =============================================================================
# Experiment III -- independent R Monte Carlo with 26 named methods.
# IMPORTANT: published-method-inspired routines in the source MATLAB program
# include explicit 1-D adaptations and provisional approximations. Those
# caveats continue here. The R port is not validation of original author code.
# =============================================================================

scope_e3_names <- c("SCOPE_Fixed_Normal","NSCOPE_Fixed_Normal_r1",
  "SCOPE_SURE_Normal","NSCOPE_SURE_Normal",
  "BAMS_Author_C2_Gamma2","EBayesMedian","RaisedCosineBayes",
  "EpanechnikovBayes","SmoothSCAD_stabilized_author","SureShrink_Hybrid",
  "FDR_BH_q05","EBBlock_1D_proxy","BayesShrink_1D","NeighCoeff",
  "NeighShrink_SURE_1D","SURELET_1D_orthonormal","BiShrink_1D",
  "NeighBlock","SureBlock_hybrid_CZ09","BlockJS",
  "Jeffreys_ABE_FigueiredoNowak","LocalWiener_1D",
  "Context_CYV2000_1D","HMT_1D_EM_proxy","UniversalHard","UniversalSoft")

scope_soft <- function(x,T) sign(x)*pmax(abs(x)-T,0)
scope_bayes_shrink <- function(x,sigma) {
  v <- sqrt(pmax(mean(x^2)-sigma^2,0))
  if (v<=sqrt(.Machine$double.eps)*max(1,sigma)) return(rep(0,length(x)))
  scope_soft(x,sigma^2/v)
}

# Author BAMS: Cscale=2, gamma=2, decomposition depth J-3, tau_j=tau_u sqrt(j).
# The input 'd' is already the author's Wavmat transform in the same ordering.
scope_bams_bayesrule <- function(d,mu,tau,epsi) {
  aa <- sqrt(2*mu); e0 <- exp(-abs(d)*aa); e1 <- exp(-abs(d)/tau)
  de <- .5*aa*e0
  den0 <- 2*tau^2-1/mu
  m <- (tau*e1-e0/aa)/den0
  c0 <- tau^2-1/(2*mu)
  if (d>=0) {
    br <- (tau*c0*d*e1+tau^2/mu*(e0-e1)) /
      (c0*(tau*e1-e0/aa))
  } else {
    br <- (tau*c0*d*e1-tau^2/mu*(e0-e1)) /
      (c0*(tau*e1-e0/aa))
  }
  out <- (1-epsi)*m*br/((1-epsi)*m+epsi*de)
  if (!is.finite(out)) out <- 0
  out
}

scope_bams <- function(d,Cscale=2,gam=2,j0=3L) {
  n <- length(d);J <- as.integer(round(log2(n)))
  fine <- sort(d[(n%/%2L+1L):n]);mm <- length(fine)
  lo <- floor(mm/4);hi <- floor(3*mm/4)
  pseudos <- abs(fine[hi]-fine[lo])/Cscale
  mu <- 1/max(pseudos^2,.Machine$double.xmin)
  tauu <- sqrt(max(stats::sd(d)^2-1/mu,.0001));out <- d
  for (j in seq(J-1L,j0,by=-1L)) {
    epsj <- 1-1/(j-j0+1)^gam
    tauj <- tauu*sqrt(j)
    ind <- (2^j+1L):2^(j+1L)
    out[ind] <- vapply(d[ind],scope_bams_bayesrule,numeric(1),
                       mu=mu,tau=tauj,epsi=epsj)
  }
  out
}

# Chang--Yu--Vetterli context grouping: user's updated 1D algorithm.
# Four neighbors + parent, no intercept in LS, coefficient-level context ranks.
scope_cyv_level <- function(x,parent,sigma) {
  n <- length(x); ax <- abs(x)
  U <- cbind(scope_circshift(ax,1),scope_circshift(ax,-1),
             scope_circshift(ax,2),scope_circshift(ax,-2))
  if (!is.null(parent)) U <- cbind(U,abs(rep(parent,each=2)))
  p <- ncol(U)
  if (n<2*p+2) return(scope_bayes_shrink(x,sigma))
  # QR least squares, pseudoinverse-like fallback via SVD in rank loss.
  w <- tryCatch(qr.solve(U,ax,tol=1e-12),error=function(e) {
    s <- svd(U);tol <- max(dim(U))*max(s$d)*.Machine$double.eps
    s$v %*% ((1/pmax(s$d,tol))*as.vector(crossprod(s$u,ax)))
  })
  Zc <- as.vector(U%*%w);order <- order(Zc,seq_len(n));xs <- x[order]
  L <- min(n,max(50L,ceiling(.02*n)))
  cs <- c(0,cumsum(xs^2));i <- seq_len(n)
  low <- pmin(pmax(i-floor(L/2),1L),n-L+1L)
  high <- low+L-1L
  sigY2 <- (cs[high+1L]-cs[low])/L
  sigX <- sqrt(pmax(sigY2-sigma^2,0));T <- rep(Inf,n)
  T[sigX>0] <- sigma^2/sigX[sigX>0]
  est <- scope_soft(xs,T);out <- numeric(n);out[order] <- est;out
}

scope_cyv <- function(d,sigma,j0=3L) {
  n <- length(d);J <- as.integer(round(log2(n)));out <- d
  for (j in j0:(J-1L)) {
    ind <- (2^j+1L):2^(j+1L)
    parent <- if (j>j0) d[(2^(j-1L)+1L):2^j] else NULL
    out[ind] <- scope_cyv_level(d[ind],parent,sigma)
  }
  out
}

# EBayes posterior median, Gaussian observation and Laplace spike/slab.
scope_laplace_parts <- function(x,a) {
  x <- abs(x)
  Ap <- .5*a*exp(.5*a^2-a*x)*pnorm(x-a)
  Am <- .5*a*exp(.5*a^2+a*x)*pnorm(-x-a)
  Ap[!is.finite(Ap)] <- 0;Am[!is.finite(Am)] <- 0
  list(ph=dnorm(x),g=pmax(Ap+Am,.Machine$double.xmin),Ap=Ap,Am=Am)
}

scope_ebayes_median <- function(x,sigma,a=.5) {
  z <- x/sigma;v <- abs(z); parts <- scope_laplace_parts(v,a)
  tU <- sqrt(2*log(max(length(x),2)))
  tpart <- scope_laplace_parts(tU,a)
  wMin <- tpart$ph/(tpart$ph+tpart$Ap-tpart$Am)
  wMin <- min(max(wMin,1e-8),1-1e-8)
  score <- function(w) sum((parts$g-parts$ph)/pmax((1-w)*parts$ph+w*parts$g,
                                               .Machine$double.xmin))
  if (score(wMin)<=0) w <- wMin else if(score(1-1e-10)>=0) w <- 1-1e-10
  else w <- uniroot(score,c(wMin,1-1e-10),tol=1e-10)$root
  marg <- (1-w)*parts$ph+w*parts$g
  wpost <- w*parts$g/pmax(marg,.Machine$double.xmin)
  isZero <- wpost*parts$Ap/pmax(parts$g,.Machine$double.xmin)<=.5
  out <- numeric(length(x));inds <- which(!isZero & v>0)
  if(length(inds)) {
    xx <- v[inds];wp <- wpost[inds]
    D <- exp(-a*xx)*pnorm(xx-a)+exp(a*xx)*pnorm(-xx-a)
    uu <- pmin(1-1e-14,pmax(1e-14,D*exp(a*xx)/(2*wp)))
    out[inds] <- pmax(xx-a-qnorm(uu),0)
  }
  sign(z)*out*sigma
}

# 64-point Gauss--Legendre quadrature (Golub--Welsch), base R eigen().
scope_gl64 <- local({cache <- NULL;function() {
  if(!is.null(cache)) return(cache)
  n <- 64L;k <- seq_len(n-1L);beta <- k/sqrt(4*k^2-1)
  A <- matrix(0,n,n)
  for(i in k) {A[i,i+1L] <- beta[i];A[i+1L,i] <- beta[i]}
  e <- eigen(A,symmetric=TRUE)
  ord <- order(e$values)
  cache <<- list(x=e$values[ord],w=2*e$vectors[1L,ord]^2)
  cache
}})

scope_quad_posterior <- function(x,sigma,alpha,tau,kind="raised",lambda=1) {
  gl <- scope_gl64();theta <- tau*gl$x; qw <- tau*gl$w
  if(kind=="raised") {
    prior <- (1+cos(pi*theta/tau))/(2*tau)
    E <- exp(-.5*(outer(theta,x,"-")/sigma)^2)
    zero <- exp(-.5*(x/sigma)^2)
  } else {
    prior <- (3/(4*tau))*(1-(theta/tau)^2)
    rate <- sqrt(2*lambda)
    E <- exp(-rate*abs(outer(theta,x,"-")))
    zero <- exp(-rate*abs(x))
  }
  weights <- qw*prior
  slm <- as.vector(crossprod(weights,E))
  sln <- as.vector(crossprod(weights*theta,E))
  (1-alpha)*sln/pmax(alpha*zero+(1-alpha)*slm,.Machine$double.xmin)
}

scope_scad_apply <- function(x,lam,a) {
  if(lam<=0) return(x)
  ab <- abs(x);out <- numeric(length(x))
  mid <- which(ab>lam & ab<a*lam);hi <- which(ab>=a*lam)
  if(length(mid)) {
    s <- (ab[mid]-lam)/((a-1)*lam)
    out[mid] <- x[mid]-sign(x[mid])*lam*cos(.5*pi*s)^2
  }
  out[hi] <- x[hi];out
}

scope_scad_sure <- function(x,sigma,N) {
  # AUTHOR-stabilized variant: SURE over a x lambda, lambda floor = 2 sigma.
  cap <- sigma*sqrt(2*log(N));lamGrid <- c(0,seq(cap/121,cap,length.out=121))
  aGrid <- c(2,2.5,3,3.5);best <- Inf;bestA <- 2;bestT <- 0
  ab <- abs(x);sig2 <- sigma^2
  for(a in aGrid) for(lam in lamGrid) {
    if(lam<=0) {sure <- length(x)*sig2} else {
      low <- ab<=lam; mid <- ab>lam & ab<a*lam;hi <- ab>=a*lam
      F <- numeric(length(x)); F[low] <- x[low]^2-sig2
      if(any(mid)) {
        s <- (ab[mid]-lam)/((a-1)*lam)
        h <- lam*cos(.5*pi*s)^2
        hp <- -pi/(2*(a-1))*sin(pi*s)
        F[mid] <- h^2+sig2-2*sig2*hp
      }
      F[hi] <- sig2; sure <- sum(F)
    }
    if(sure<best) {best <- sure;bestA <- a;bestT <- lam}
  }
  scope_scad_apply(x,min(cap,max(bestT,2*sigma)),bestA)
}

scope_sureshrink <- function(x,sigma) {
  n <- length(x); z <- sort(abs(x)/sigma);cap <- sqrt(2*log(max(n,2)))
  gam <- log2(max(n,2))^1.5/sqrt(n)
  if (mean(z^2)-1<=gam) return(scope_soft(x,sigma*cap))
  ids <- which(z<=cap);cs <- cumsum(z^2)
  # Compute final universal-cap candidate separately to avoid precedence.
  capRisk <- n-2*length(ids)+(if(length(ids)) cs[max(ids)] else 0)+
             (n-length(ids))*cap^2
  risk <- c(n,n-2*ids+cs[ids]+(n-ids)*z[ids]^2,capRisk)
  tvals <- c(0,z[ids],cap)
  scope_soft(x,sigma*tvals[which.min(risk)])
}

scope_rollmean <- function(x,w=7L) {
  r <- (w-1L)%/%2L;v <- numeric(length(x))
  for(lag in -r:r) v <- v+scope_circshift(x,lag)
  v/w
}

scope_neighshrink_sure <- function(x,sigma,N) {
  sig2 <- sigma^2;ref <- sigma*sqrt(2*log(N));best <- Inf;out <- x*0
  for(r in 1:3) {
    e <- x^2
    for(lag in seq_len(r)) e <- e+scope_circshift(x^2,lag)+scope_circshift(x^2,-lag)
    for(b in c(.5,.75,1,1.25)) {
      T <- (b*ref)^2; w <- pmax(1-T/pmax(e,.Machine$double.xmin),0)
      z <- x*w;div <- numeric(length(x));ix <- which(e>T)
      div[ix] <- w[ix]+2*T*x[ix]^2/e[ix]^2
      sure <- sum((z-x)^2+2*sig2*div-sig2)
      if(sure<best) {best <- sure;out <- z}
    }
  }
  out
}

scope_block_js <- function(x,sigma2,N) {
  out <- x*0;L <- max(1L,floor(log(N)))
  for(i in seq(1L,length(x),by=L)) {
    idx <- i:min(length(x),i+L-1L);E <- sum(x[idx]^2)
    out[idx] <- x[idx]*pmax(1-4.50524*length(idx)*sigma2/pmax(E,.Machine$double.xmin),0)
  };out
}

scope_neighblock <- function(x,sigma2,N) {
  n <- length(x);L0 <- max(1L,floor(log(N)/2));L1 <- max(1L,floor(L0/2))
  L <- L0+2*L1;out <- x*0
  for(i in seq(1L,n,by=L0)) {
    idx <- i:min(n,i+L0-1L)
    win <- ((seq(i-L1,i+L0+L1-1L)-1L) %% n)+1L
    E <- sum(x[win]^2)
    out[idx] <- x[idx]*pmax(1-4.50524*L*sigma2/pmax(E,.Machine$double.xmin),0)
  };out
}

scope_ebblock <- function(x,sigma2) {
  n <- length(x);z <- x;L <- 4L
  tau2 <- max(mean(x^2)-sigma2,1e-5*sigma2)
  pi1 <- min(.95,max(.05,mean(x^2>2*sigma2)))
  for(i in seq(1L,n,by=L)) {
    idx <- i:min(i+L-1L,n);a <- x[idx];m <- length(idx)
    b0 <- exp(-sum(a^2)/(2*sigma2))/(sigma2^(m/2))
    b1 <- exp(-sum(a^2)/(2*(sigma2+tau2)))/((sigma2+tau2)^(m/2))
    post <- (pi1*b1)/(pi1*b1+(1-pi1)*b0+.Machine$double.xmin)
    z[idx] <- a*post*tau2/(sigma2+tau2)
  };z
}

scope_surelet <- function(x,sigma) {
  sigma2 <- sigma^2;e <- exp(-x^2/(12*sigma2))
  A <- cbind(x,x*e);div <- cbind(rep(1,length(x)),e*(1-x^2/(6*sigma2)))
  g <- crossprod(A);rhs <- as.vector(crossprod(A,x)-sigma2*colSums(div))
  a <- tryCatch(solve(g,rhs),error=function(e) {
    ss <- svd(g); tol <- 1e-12*max(ss$d)
    as.vector(ss$v %*% ((ifelse(ss$d>tol,1/ss$d,0))*
                        as.vector(crossprod(ss$u,rhs))))
  })
  as.vector(A%*%a)
}

scope_sureblock <- function(x,sigma2) {
  # Cai--Zhou hybrid (R implementation of integer block/threshold search).
  n <- length(x);u <- x/sqrt(sigma2)
  if(n<2) return(x)
  if(mean(u^2)-1<=log2(n)^1.5/sqrt(n))
    return(x*pmax(1-2*log(n)/pmax(u^2,.Machine$double.xmin),0))
  best <- Inf;chosenL <- 1L;chosenLam <- 0
  for(L in seq_len(floor(sqrt(n)))) {
    ix <- split(seq_len(n),ceiling(seq_len(n)/L))
    E <- vapply(ix,function(ids) sum(u[ids]^2),numeric(1))
    sizes <- lengths(ix);lo <- max(L-2,0);hi <- 2*L*log(n)
    # Piecewise quadratic in lambda; evaluate knots and active-set minima.
    cand <- sort(unique(c(lo,hi,E[E>=lo & E<=hi])))
    if(length(cand)>1) for(ii in seq_len(length(cand)-1L)) {
      a <- cand[ii];b <- cand[ii+1L];active <- which(E>(a+b)/2)
      if(length(active)) {
        stationary <- sum((sizes[active]-2)/pmax(E[active],.Machine$double.xmin)) /
                       sum(1/pmax(E[active],.Machine$double.xmin))
        cand <- c(cand,pmin(b,pmax(a,stationary)))
      }
    }
    for(lam in unique(cand)) {
      inactive <- E<=lam;active <- !inactive
      risk <- n+sum(E[inactive]-2*sizes[inactive])+sum(
        (lam^2/pmax(E[active],.Machine$double.xmin))-
          2*lam*(sizes[active]-2)/pmax(E[active],.Machine$double.xmin))
      if (risk<best) {best <- risk;chosenL <- L;chosenLam <- lam}
    }
  }
  z <- numeric(n)
  for(ids in split(seq_len(n),ceiling(seq_len(n)/chosenL))) {
    E <- sum(u[ids]^2)
    z[ids] <- x[ids]*pmax(1-chosenLam/pmax(E,.Machine$double.xmin),0)
  };z
}

# A two-state Gaussian 1D HMT proxy. Upward messages, downward approximate
# posterior updates, empirical variance/transition iterations: NOT exact HMT
# reference software. Trace its approximation status in METHOD_INDEX.csv.
scope_hmt_proxy <- function(d,sigma2,j0=3L) {
  J <- as.integer(round(log2(length(d))));levels <- lapply(j0:(J-1L),function(j)
    d[(2^j+1L):2^(j+1L)])
  tau <- c(.05,max(mean(d[(2^j0+1L):length(d)]^2)-sigma2,.2*sigma2)*4)
  pi0 <- c(.85,.15);A <- matrix(c(.94,.14,.06,.86),2,2)
  gamma <- vector("list",length(levels))
  for(iter in seq_len(15L)) {
    upward <- vector("list",length(levels))
    for(jj in rev(seq_along(levels))) {
      x <- levels[[jj]];v <- sigma2+tau
      B <- rbind(exp(-x^2/(2*v[1])-.5*log(v[1])),
                 exp(-x^2/(2*v[2])-.5*log(v[2])))
      B <- sweep(B,2L,pmax(apply(B,2,max),.Machine$double.xmin),"/")
      if (jj<length(levels)) {
        child <- upward[[jj+1L]];msg <- A%*%child
        B <- B*msg[,seq(1L,ncol(msg),by=2L)]*
               msg[,seq(2L,ncol(msg),by=2L)]
      }
      upward[[jj]] <- sweep(pmax(B,.Machine$double.xmin),2L,
                           pmax(colSums(B),.Machine$double.xmin),"/")
    }
    root <- upward[[1L]]*pi0
    gamma[[1L]] <- sweep(root,2L,pmax(colSums(root),.Machine$double.xmin),"/")
    counts <- matrix(0,2,2)
    if(length(levels)>1) for(jj in 2:length(levels)) {
      B <- upward[[jj]];parent <- ceiling(seq_len(ncol(B))/2)
      old <- gamma[[jj-1L]][,parent,drop=FALSE]
      denom <- A%*%B;xi <- matrix(0,2,ncol(B))
      for(s in 1:2) {
        v <- sweep(B,1L,A[s,],"*")
        v <- sweep(v,2L,old[s,]/pmax(denom[s,],.Machine$double.xmin),"*")
        counts[s,] <- counts[s,]+rowSums(v)
        xi <- xi+v
      }
      gamma[[jj]] <- sweep(xi,2L,pmax(colSums(xi),.Machine$double.xmin),"/")
    }
    pi0 <- pmax(rowMeans(gamma[[1L]]),1e-5);pi0 <- pi0/sum(pi0)
    A <- pmax(counts,1e-6);A <- A/rowSums(A)
    for(s in 1:2) {
      num <- 0;den <- 0
      for(jj in seq_along(levels)) {
        ww <- gamma[[jj]][s,];num <- num+sum(ww*levels[[jj]]^2);den <- den+sum(ww)
      }
      tau[s] <- max(num/pmax(den,1e-6)-sigma2,.01*sigma2)
    }
    tau[2] <- max(tau[2],tau[1]+.05*sigma2)
  }
  out <- d
  for(jj in seq_along(levels)) {
    j <- j0+jj-1L;ind <- (2^j+1L):2^(j+1L)
    w <- gamma[[jj]][1,]*tau[1]/(sigma2+tau[1])+
         gamma[[jj]][2,]*tau[2]/(sigma2+tau[2])
    out[ind] <- levels[[jj]]*w
  }
  out
}

scope_e3_comparators <- function(d,sigma,N) {
  # Returns columns 5--24; order identical to MATLAB table.
  J <- as.integer(round(log2(N)));j0 <- 3L; sig2 <- sigma^2
  out <- matrix(rep(d,20L),nrow=N,ncol=20L)
  out[,1L] <- scope_bams(d,Cscale=2,gam=2,j0=3L)
  # Method 23 uses all levels at once to access parent coefficients.
  out[,19L] <- scope_cyv(d,sigma,j0=3L)
  out[,20L] <- scope_hmt_proxy(d,sig2,j0=3L)
  fine <- d[(N/2+1L):N];fineSD <- max(sd(fine),1e-8)
  epLambda <- 1/fineSD^2+.5*exp(-fineSD/2)
  tauRC <- max(max(abs(d[9L:N])),1e-8)
  # BH-FDR over all details (method 11).
  idx <- (2^j0+1L):N
  p <- 2*pnorm(-abs(d[idx])/sigma);ps <- sort(p)
  K <- length(ps);found <- which(ps<=.05*seq_len(K)/K)
  cut <- if(length(found)) ps[max(found)] else -1
  out[idx,7L] <- d[idx]*(p<=cut)
  for (j in j0:(J-1L)) {
    ids <- (2^j+1L):2^(j+1L);x <- d[ids]
    alpha <- pmin(pmax(1-1/(j-j0+1)^2,0),1-1e-12)
    out[ids,2L] <- scope_ebayes_median(x,sigma,a=.5)
    out[ids,3L] <- scope_quad_posterior(x,sigma,alpha,tauRC,"raised")
    out[ids,4L] <- scope_quad_posterior(x,sigma,alpha,max(abs(x),1e-8),
                                        "epanech",lambda=epLambda)
    out[ids,5L] <- scope_scad_sure(x,sigma,N)
    out[ids,6L] <- scope_sureshrink(x,sigma)
    out[ids,8L] <- scope_ebblock(x,sig2)
    out[ids,9L] <- scope_bayes_shrink(x,sigma)
    eng <- x^2+scope_circshift(x^2,1)+scope_circshift(x^2,-1)
    out[ids,10L] <- x*pmax(1-2*log(N)*sig2/pmax(eng,.Machine$double.xmin),0)
    out[ids,11L] <- scope_neighshrink_sure(x,sigma,N)
    out[ids,12L] <- scope_surelet(x,sigma)
    parent <- if(j==j0) rep(0,length(x)) else
      rep(d[(2^(j-1L)+1L):2^j],each=2L)
    v <- pmax(scope_rollmean(x^2,7L)-sig2,1e-10)
    T <- sqrt(3)*sig2/sqrt(v);R <- sqrt(x^2+parent^2)
    out[ids,13L] <- x*pmax(1-T/pmax(R,.Machine$double.xmin),0)
    out[ids,14L] <- scope_neighblock(x,sig2,N)
    out[ids,15L] <- scope_sureblock(x,sig2)
    out[ids,16L] <- scope_block_js(x,sig2,N)
    out[ids,17L] <- x*pmax(1-3*sig2/pmax(x^2,.Machine$double.xmin),0)
    v <- pmax(scope_rollmean(x^2,7L)-sig2,0)
    out[ids,18L] <- x*v/pmax(v+sig2,.Machine$double.xmin)
  }
  scope_assert(all(is.finite(out)), "Nonfinite Experiment III comparator")
  out
}

scope_e3_condition <- function(signal,si,snr,qi,N,M,cal,outDir,saveEvery=25L) {
  setup <- scope_condition_setup(signal,N,snr)
  cand <- scope_candidates();names <- scope_e3_names
  path <- file.path(outDir,sprintf("checkpoint_%s_SNR%d.rds",signal,snr))
  set.seed(scope_seed(N,"Exp3",si,snr))
  signature <- list(version="R_E3_20261029_1",N=N,M=M,signal=signal,
                    snr=snr,eta=cal$eta,methodNames=names)
  init <- list(completed=0L,mse=matrix(NA_real_,M,26L,dimnames=list(NULL,names)),
               sigma=rep(NA_real_,M),selection=matrix(NA_real_,M,5L))
  state <- scope_run_state(path,init,signature)
  if (state$completed<M) for (it in seq.int(state$completed+1L,M)) {
    y <- setup$clean+rnorm(N);d <- scope_dwt(y,setup$filters)
    sigma <- scope_sigma(d);blocks <- scope_neighborhood(d,rGrid=0:3)
    scores <- vapply(seq_len(nrow(cand)),function(c)
      scope_sure(blocks,sigma,cal,"Normal",cand$k[c],cand$r[c],cand$b[c]),numeric(1))
    iS <- which(cand$r==0); i0 <- iS[which.min(scores[iS])];i1 <- which.min(scores)
    est <- matrix(NA_real_,N,26L)
    est[,1L] <- scope_apply(d,blocks,sigma,cal,"Normal",5L,0L,1)
    est[,2L] <- scope_apply(d,blocks,sigma,cal,"Normal",5L,1L,1)
    est[,3L] <- scope_apply(d,blocks,sigma,cal,"Normal",cand$k[i0],0L,cand$b[i0])
    est[,4L] <- scope_apply(d,blocks,sigma,cal,"Normal",cand$k[i1],cand$r[i1],cand$b[i1])
    est[,5:24] <- scope_e3_comparators(d,sigma,N)
    est[,25L] <- scope_universal(d,sigma,kind="hard")
    est[,26L] <- scope_universal(d,sigma,kind="soft")
    state$mse[it,] <- colMeans((est-setup$truth)^2)
    state$sigma[it] <- sigma
    state$selection[it,] <- c(cand$k[i0],cand$b[i0],cand$k[i1],cand$r[i1],cand$b[i1])
    state$completed <- it
    if (it%%saveEvery==0L || it==M) {
      scope_write_state(state,path)
      cat(sprintf("Exp3 N=%d %-10s SNR=%d: %d/%d\n",N,signal,snr,it,M))
    }
  }
  write.csv(data.frame(Replication=seq_len(M),SigmaHat=state$sigma,state$mse,
                       check.names=FALSE),file.path(outDir,
                         sprintf("MSE_replications_%s_SNR%d.csv",signal,snr)),row.names=FALSE)
  mu <- colMeans(state$mse);dev <- apply(state$mse,2L,sd)
  tab <- data.frame(N=N,Signal=signal,SNR=snr,Method=names,M=M,AMSE=mu,
                    SD_MSE=dev,SE_AMSE=dev/sqrt(M),RMSE=sqrt(mu),row.names=NULL)
  scope_e3_plots(tab,state$mse,file.path(outDir,"plots"))
  list(summary=tab,mse=state$mse)
}

scope_run_e3 <- function(N,snrs,M,root,smoke=FALSE) {
  out <- scope_make_dir(file.path(root,sprintf("Exp3_N%d_M%d",N,M)))
  cal <- scope_get_calibration(file.path(out,"calibration_E3_Normal.rds"),N,
                               "Normal",c(3L,5L,7L),0:3)
  signals <- if(smoke) c("Doppler","Blocks") else scope_conditions
  all <- list();pairs <- list()
  for(signal in signals) for(qi in seq_along(snrs)) {
    snr <- snrs[qi]
    rr <- scope_e3_condition(signal,match(signal,scope_conditions),snr,qi,N,M,cal,out)
    all[[length(all)+1L]] <- rr$summary
    for(k in 1:4) for(oth in 5:24) {
      delta <- rr$mse[,oth]-rr$mse[,k]
      pairs[[length(pairs)+1L]] <- data.frame(N=N,Signal=signal,SNR=snr,
        SCOPE=scope_e3_names[k],Other=scope_e3_names[oth],
        OtherMinusSCOPE=mean(delta),PairedSE=sd(delta)/sqrt(M),
        CI_L=mean(delta)-1.96*sd(delta)/sqrt(M),
        CI_U=mean(delta)+1.96*sd(delta)/sqrt(M))
    }
  }
  tab <- do.call(rbind,all)
  write.csv(tab,file.path(out,"Experiment3_AMSE_all.csv"),row.names=FALSE)
  write.csv(do.call(rbind,pairs),file.path(out,"Experiment3_PairedDifferences.csv"),row.names=FALSE)
  ref <- subset(tab,Method=="UniversalHard",c("Signal","SNR","AMSE"))
  rank <- merge(tab,ref,by=c("Signal","SNR"),suffixes=c("","_Hard"))
  rank$Ratio <- rank$AMSE/rank$AMSE_Hard
  rank <- aggregate(Ratio~Method,rank,mean)
  rank <- rank[order(rank$Ratio),]
  write.csv(rank,file.path(out,"Experiment3_AggregateRanking.csv"),row.names=FALSE)
  scope_e3_aggregate_plot(rank,file.path(out,"plots"))
  scope_method_index(file.path(out,"METHOD_IMPLEMENTATION_NOTES.csv"))
  cat("Experiment III complete:",out,"\n")
  invisible(tab)
}

scope_method_index <- function(file) {
  status <- rep("Published-method-inspired 1D adaptation; independent audit needed",26L)
  status[1:4] <- "SCOPE fixed null-risk calibration / plug-in SURE, as specified"
  status[5] <- "Author BAMS.m: Cscale=2; gamma=2; depth J-3; tau_j=tau_u sqrt(j)"
  status[9] <- "AUTHOR stabilized SmoothSCAD, not plain unrestricted SURE"
  status[c(12,24)] <- "PROVISIONAL 1D proxy; not validated published software"
  status[16] <- "1D orthonormal two-basis SURE-LET adaptation"
  status[23] <- "CYV 1D: 4 neighbors+parent, LS context, 50 or 2% rank window"
  status[25:26] <- "Universal hard / soft wavelet threshold"
  write.csv(data.frame(Index=seq_len(26L),Procedure=scope_e3_names,
                       ValidationStatus=status),file,row.names=FALSE)
}


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


SCOPE_Exp123_FINAL <- function(mode="full",experiment="all",
                               Ns=c(512L,1024L,2048L),snrs=c(3,5,7),M=1000L) {
  # 'smoke' is a small functionality test, never publication evidence.
  mode <- tolower(mode);experiment <- tolower(experiment)
  scope_assert(mode%in%c("full","smoke"),"mode must be full or smoke")
  scope_assert(experiment%in%c("all","exp1","exp2","exp3"),
               "experiment must be all, Exp1, Exp2, or Exp3")
  if(mode=="smoke") {Ns <- 512L;snrs <- 5;M <- 2L}
  Ns <- unique(as.integer(Ns));snrs <- unique(as.numeric(snrs));M <- as.integer(M)
  scope_assert(length(Ns)>0 && all(Ns>=16 & Ns==2^round(log2(Ns))),
               "N values must be powers of two >=16")
  scope_assert(length(snrs)>0 && all(is.finite(snrs)&snrs>0),
               "Root-SNR values must be positive")
  scope_assert(length(M)==1 && M>0, "M must be positive")
  # Freeze R's RNG implementation for reproducible reruns with R >= 3.6.
  RNGkind(kind="Mersenne-Twister",normal.kind="Inversion",sample.kind="Rejection")
  root <- scope_make_dir(file.path(.scope_dir,"SCOPE_Exp123_R_results"))
  cat(sprintf("R port: %s | N=%s | root-SNR=%s | M=%d\n",
        experiment,paste(Ns,collapse=","),paste(snrs,collapse=","),M))
  cat("Noise: iid N(0,1); root SNR=sd(signal,population)/1, NOT power SNR.\n")
  cat("R random-number streams differ from MATLAB; historical .mat checkpoints cannot be reused.\n")
  for(N in Ns) {
    Fs <- if(mode=="smoke") c("Normal","Logistic") else scope_generators
    if(experiment%in%c("all","exp1")) scope_run_e1(N,snrs,M,Fs,root,mode=="smoke")
    if(experiment%in%c("all","exp2")) scope_run_e2(N,snrs,M,Fs,root,mode=="smoke")
    if(experiment%in%c("all","exp3")) scope_run_e3(N,snrs,M,root,mode=="smoke")
  }
  if(mode=="full" && experiment%in%c("all","exp1"))
    scope_plot_N_SNR_sensitivity(root,Ns,snrs,M)
  cat("Completed. Results:",root,"\n")
  invisible(root)
}
