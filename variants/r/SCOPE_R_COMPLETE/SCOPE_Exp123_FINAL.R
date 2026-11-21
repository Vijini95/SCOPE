# =============================================================================
# CANONICAL R ENTRY POINT -- source("SCOPE_Exp123_FINAL.R")
# Use SCOPE_Exp123_FINAL(mode,experiment,Ns,snrs,M).
# Scripts 01--04 are loaded relative to THIS file (not R's working directory).
# Output: SCOPE_Exp123_R_results; checkpoints are .rds, figures .pdf + .png.
# =============================================================================

.scope_this_file <- function() {
  fr <- sys.frames()
  for (i in rev(seq_along(fr))) {
    if (!is.null(fr[[i]]$ofile) && grepl("SCOPE_Exp123_FINAL\\.R$",fr[[i]]$ofile))
       return(normalizePath(fr[[i]]$ofile))
  }
  normalizePath("SCOPE_Exp123_FINAL.R",mustWork=TRUE)
}
.scope_dir <- dirname(.scope_this_file())
for (f in c("01_scope_core.R","04_scope_figures.R",
            "02_scope_experiments12.R","03_scope_experiment3.R")) {
  source(file.path(.scope_dir,f),local=FALSE)
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
