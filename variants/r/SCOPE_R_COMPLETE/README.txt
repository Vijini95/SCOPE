SCOPE: ANNOTATED R TRANSLATION + AUTOMATIC FIGURES
=================================================

Original source: uploaded SCOPE_Exp123(1).m, containing the canonical
MATLAB SCOPE_Exp123_FINAL entry point and all three experiments.
Translation: base R; no extra CRAN package, Wavelet Toolbox, or MATLAB needed.

1. CONTENTS

  SCOPE_Exp123_FINAL.R     main entry point: loads all support files
  01_scope_core.R          wavelet transform, Donoho--Johnstone signals,
                           null-risk calibration, SCOPE, N-SCOPE, SURE
  02_scope_experiments12.R Experiments I / II with resuming, CSV and plots
  03_scope_experiment3.R  26-method benchmark, author BAMS, new CYV context
  04_scope_figures.R       plots in BOTH PDF and PNG, including Figure 2
  RUN_SMOKE.R             transform tests, eta check, 2-replicate smoke run
  RUN_FULL.R              full 3-experiment paper design, M=1000

2. HOW TO RUN (from R or RStudio)

Extract all files to a single folder, set that folder as your working
folder in R or RStudio. Run:

  source("RUN_SMOKE.R")

Alternatively:

  source("SCOPE_Exp123_FINAL.R")
  SCOPE_Exp123_FINAL("smoke", "all")
  SCOPE_Exp123_FINAL("full", "Exp1", c(512,1024,2048), c(3,5,7), 1000)
  SCOPE_Exp123_FINAL("full", "Exp2", c(512,1024,2048), c(3,5,7), 1000)
  SCOPE_Exp123_FINAL("full", "Exp3", c(512,1024,2048), c(3,5,7), 1000)

To run only one combination:

  SCOPE_Exp123_FINAL("full", "Exp1", 1024, 5, 1000)

NOTE: N=512/1024/2048 and 'SNR'=3/5/7 are DEFAULT. The model uses
ROOT (standard-deviation-ratio) SNR, not the usual squared-power ratio.
Noise variance remains 1. Data are not standardized per replication.

3. OUTPUTS & FIGURES

All files are generated below:
  SCOPE_Exp123_R_results/
    Exp1_N1024_M1000/
      Experiment1_AMSE_all.csv, Experiment1_PairedGains.csv,
      Experiment1_Constants.csv, Experiment1_AggregateRanking.csv,
      MSE_replications_<signal>_SNR<value>.csv,
      checkpoint_<signal>_SNR<value>.rds,
      plots/E1_AMSE_*.pdf, .png and E1_GainPercent_*.pdf, .png

    Exp2_N1024_M1000/
      Experiment2_AMSE_all.csv, Experiment2_PairedDifferences.csv,
      Experiment2_SelectedTriples.csv,
      Experiment2_SelectionFrequencies.csv, Experiment2_AggregateRanking.csv,
      plots/E2_AMSE_*, E2_RadiusSelection_*, E2_AggregateRanking_*,
            E2_RadiusSelection_BySignal_*, each as PDF and PNG

    Exp3_N1024_M1000/
      Experiment3_AMSE_all.csv, Experiment3_PairedDifferences.csv,
      Experiment3_AggregateRanking.csv, METHOD_IMPLEMENTATION_NOTES.csv,
      plots/E3_AMSE_*, E3_MSEBox_*, E3_AggregateRanking_* (PDF & PNG)

    plots/E1_Figure2_Sensitivity_N_SNR_M1000.pdf/.png
      (created when Exp1 results at 2+ values of N and SNR exist)

Results are NOT silently imported from an older MATLAB simulation.
Checkpoints are distinct .rds R objects and cannot be mixed with .mat.
Running the same R condition again resumes from its saved state.
If you change algorithms, rename or archive its output folder first.

4. EXACTNESS / CAUTIONS

Mathematical SCOPE specifications in Exp1/2 and four SCOPE procedures
in Exp3 follow the supplied MATLAB, including the periodic wavelet filter
phase convention, chi-squared Gaussian null risk, eta(N), SURE derivative,
and level scheduling. BAMS uses the original BAMS.m parameterization:
Cscale=2, gam=2, coarsest level j0=3, depth L=J-3, tau_j=tau_u sqrt(j).
CYV context is the user's updated 1D orthonormal context regression version
(four lateral neighbors, parent, LS weights and context-rank grouping).

Experiment III includes numerous publications-derived MATLAB comparators,
including explicitly PROVISIONAL 1D approximations. Their R translations
are computational ports, not independent validations against the research
papers. The HMT routine is likewise a 1D EM-inspired proxy. Do not label
all competitor rows as exact published-author implementations.

R and MATLAB generate different Gaussian random sequences for a given
seed and have different floating-point integration behavior. Thus the
R AMSE values will not exactly equal the original MATLAB table, even if
the estimator formulas are identical. A dedicated cross-language fixture
comparison is needed to establish coefficient-by-coefficient equivalence.

5. VALIDATION STATUS

Static source/code checks performed during construction. MATLAB and R
interpreters are NOT available in this environment, so the R script has
NOT been executed here. Run RUN_SMOKE.R before a publication-scale run.
In particular, carefully verify advanced Experiment III comparisons before
replacing numbers in the paper. No computed R AMSEs are claimed here.
