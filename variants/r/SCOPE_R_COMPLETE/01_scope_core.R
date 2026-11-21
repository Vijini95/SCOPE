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
