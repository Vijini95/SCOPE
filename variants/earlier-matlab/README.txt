SCOPE / N-SCOPE — CANONICAL MATLAB PACKAGE (2026-10-29)
======================================================

ONE source file: SCOPE_Exp123_FINAL.m.  It contains all transforms,
calibration, 26 Experiment III estimators, plots and checkpoint maintenance.
RUN_SMOKE.m is a short check; RUN_FULL.m is the full three-experiment run.
No external .m dependencies or paid MATLAB toolboxes are required.
MATLAB R2018b+ recommended. MATLAB execution was NOT available in the
creation environment; run RUN_SMOKE.m before publication calculations.

FULL DESIGN
===========
N = [512 1024 2048], root SNR = [3 5 7], 4 Donoho-Johnstone signals,
M=1000 noise draws per signal/N/SNR. White Gaussian noise variance 1.
ROOT SNR=std(f,1)/sigma (SD ratio), NOT E(f^2)/sigma^2.

COMMANDS (MATLAB current folder must contain this .m)
=======================================================
SCOPE_Exp123_FINAL('smoke','all')
SCOPE_Exp123_FINAL('full','Exp1',1024,[3 5 7],1000)
SCOPE_Exp123_FINAL('full','Exp2',1024,[3 5 7],1000)
SCOPE_Exp123_FINAL('full','Exp3',1024,[3 5 7],1000)
SCOPE_Exp123_FINAL('full','all',[512 1024 2048],[3 5 7],1000)

OUTPUT LOCATION AND FIGURES
===========================
SCOPE_Exp123_results/Exp1_N1024_M1000/plots/      -> E1 bars, gains
SCOPE_Exp123_results/Exp2_N1024_M1000/plots/      -> E2 bars, gains, radius selection, aggregate
SCOPE_Exp123_results/Exp3_N1024_M1000/plots/      -> E3 all 24 procedures, box summaries, aggregate
SCOPE_Exp123_results/plots/ -> E1 Figure2 across N and SNR, if all are present
All figures are generated in BOTH PNG and PDF.
All experiment folders also contain MSE CSVs, AMSE summaries and checkpoint
MAT files. Resuming a completed run reloads checkpoints and REBUILDS figures.

AUTHOR BAMS — METHOD #5 (only BAMS used)
========================================
Original user MATLAB BAMS.m conventions:
   Cscale=2, gam=2; depth L=J-3 (NOT L=3), eight scaling coefficients;
   pseudo sigma = (Q3-Q1)/Cscale; mu=1/pseudo sigma^2;
   tau_u=sqrt(max(std(all WT coefficients)^2-1/mu,0.0001));
   tau_j=tau_u*sqrt(j); epsilon_j=1-1/(j-3+1)^2;
   posterior mean from user's supplied bams_bayesrule.
No separate BAMSO reparameterization is included in the 26-method table.

CYV CONTEXT — METHOD #23 (replaces old Context_1D_proxy)
=========================================================
Chang, Yu & Vetterli, IEEE Trans. Image Process. 9(9), 1522-1531 (2000).
This is a 1D orthonormal-wavelet *adaptation*, not the original 2D image
algorithm: four cyclic lateral neighbours plus dyadic parent, LS context
regression, context-sorted local groups L=max(50, ceil(2% of level size))
limited by level size; T=sigma^2/sqrt(max(mean(context-group d^2)-sigma^2,0)).
Small levels fall back to BayesShrink. Output is single-coefficient soft
thresholding, not block shrinkage. Regression uses U\abs(d) when stable.

CHECKPOINTS — DO NOT DELETE EXISTING DATA
==========================================
Exp1 and Exp2 preserve their previous simulation designs and resume files.
Exp3 rejects unknown checkpoints. It recognizes three previous code revisions
and safely converts 26- or 28-column checkpoints into canonical 26 columns.
It backs up each old .mat to *.pre_FINAL.bak BEFORE changing anything.
Replaying all completed Gaussian noise draws verifies each sigma estimate
and the first stored noisy observation; only BAMS #5 and CYV #23 are
recomputed. The other 24 MSE columns are copied unchanged. If verification
fails, conversion aborts without modifying the original checkpoint.
This enables rebuilding Exp3 plots and rankings without recomputing the
other 24 comparators. Do not discard *.pre_FINAL.bak.

LIMITATIONS / VALIDATION
========================
The BAMS and CYV source-level changes are checked structurally against
uploaded sources; final MATLAB numerical execution has NOT been performed
here. Prior audits were incomplete. Other comparison procedures include
explicitly labeled 1D adaptations and provisional proxies: do not describe
them as verified replicas of published programs without further checking.
The Monte Carlo SNR convention must match the paper, tables, and captions.
