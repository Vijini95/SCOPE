source("SCOPE_Exp123_ALL_IN_ONE.R")
for (N in c(32L,128L,512L)) for (F in c("haar","db3","sym4")) {
  x <- sin(seq_len(N)*.132) + cos(seq_len(N)*.031)
  filt <- scope_filters(F)
  d <- scope_dwt(x,filt)
  stopifnot(max(abs(scope_idwt(d,filt)-x))<1e-9)
}
for (N in c(512L,1024L,2048L))
  cat(sprintf("eta(%d) = %.8f\n",N,scope_eta(N)))
SCOPE_Exp123_FINAL("smoke","all")
