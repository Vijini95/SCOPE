source("SCOPE_Exp123_FINAL.R")
# Reconstruction/orthogonality check, independent of simulation.
for (N in c(32L,128L,512L)) for (F in c("haar","db3","sym4")) {
  x <- sin(seq_len(N)*.132) + cos(seq_len(N)*.031)
  h <- scope_filters(F)
  d <- scope_dwt(x,h)
  err <- max(abs(scope_idwt(d,h)-x))
  stopifnot(err<1e-9,abs(sum(x^2)-sum(d^2))<1e-8)
  cat(sprintf("DWT test: N=%d, %-4s, error %.3g\n",N,F,err))
}
for (N in c(512L,1024L,2048L))
  cat(sprintf("eta(%d) = %.8f\n",N,scope_eta(N)))
# Small Monte Carlo incl. all 3 experiments and PDFs/PNGs.
SCOPE_Exp123_FINAL("smoke","all")
cat("R smoke test finished. Check the plots subfolders.\n")
