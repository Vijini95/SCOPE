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
