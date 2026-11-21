SCOPE / N-SCOPE SURE — Experiment II of Section 6
===================================================

Main file: run_SCOPE_Experiment2.m
Ancillary files: Wavmat.m, MakeONFilterExt.m (supplied by the author)

Requirements: MATLAB with base functions including table/writetable,
  erf, integral, fzero, and matrix multiplication. MATLAB R2018a+
  recommended. NO Wavelet Toolbox, Statistics Toolbox, or WaveLab needed.

HOW TO RUN
  (1) Unzip ALL files into one folder. In MATLAB, cd into that folder.
  (2) Test installation:       run_SCOPE_Experiment2('smoke')
      Smoke uses Normal+Logistic, Doppler+Blocks, SNR 5, M=6; it uses the
      full k/r/b grids but is NOT appropriate for publication.
  (3) Run the full study:      run_SCOPE_Experiment2

FULL DESIGN
  n=1024; Donoho-Johnstone Doppler, HeaviSine, Bumps, Blocks
  Power SNR: 3, 5, 7; independent Gaussian noise variance 1
  M=1000 per signal-SNR cell, same noise for all methods
  Six F: Normal, Logistic, Laplace, Hyperbolic Sech, Uniform, Cauchy
  k={3,5,7}; r={0,1,2,3}; b={0.65,0.8,1,1.2,1.5}
  Null-risk target eta=0.01, eta_j=0.01*2^(9-j)
  lambda_j(F,k,r,b)=b*a(F,k,r,eta_j)/sigmahat
  Sigmahat=MAD of finest wavelet detail coefficients
  Symmlet 4 (8-tap) for Doppler/HeaviSine,
  Daubechies 3 (6-tap) for Bumps, Haar for Blocks
  Orthogonal periodic DWT using the author's Wavmat.m;
  8 scaling coefficients untouched; shrunk detail levels j=3,...,9

EACH GENERATOR HAS 4 OUTPUT RULES
  Fixed SCOPE:  k=5, r=0, b=1
  Fixed N-SCOPE: k=5, r=1, b=1
  SCOPE SURE:   SURE chooses (k,b) from the grid, r fixed at 0
  N-SCOPE SURE: SURE chooses (k,r,b) from the entire grid
  Two additional references: universal HARD and SOFT.

VERY IMPORTANT
  - All selections are global across 7 detail levels for one realization.
  - N-SCOPE SURE includes SCOPE SURE exactly as its r=0 subclass.
  - Derivative includes central RMS differentiation; overlap is allowed.
  - sigmahat is held fixed in Stein differentiation: PLUG-IN SURE.
  - No true-signal tuning and no SURE optimization of the generator F.
  - The same sample noise is used for all F and all 4 rules.
  - For fixed methods, b=1; their calibration matches Experiment I.
  - Seeds equal Experiment I so paired fixed-result checks are possible;
    for final Experiment III, use independent simulations / seeds.

OUTPUT DIRECTORY
  outputs_E2_M1000 (full) or outputs_SMOKE_NOT_FOR_PUBLICATION

KEY OUTPUTS
  Experiment2_AMSE_all.csv         - AMSE/SD/SE, gain vs fixed SCOPE
  Experiment2_PairedDifferences.csv - paired delta AMSE, SE, 95% MC CI
  Experiment2_SelectionFrequencies.csv - marginal selection of k/r/b
  Experiment2_SelectedTriples.csv - joint selected candidate frequencies
  Experiment2_AggregateRanking.csv - aggregate relative AMSE
  Experiment2_ScaleConstants.csv  - 504 noise-calibrated constants
  MSE_replications_*.csv          - replication-level MSE for every rule
  SelectedParams_*.csv            - replication-level SURE selections
  checkpoint_*.mat                - completed simulations; resume enabled
  PARTIAL_*.mat                   - safe resume during a condition
  plots/E2_AMSE_*.png/pdf         - 4 method bar comparisons per signal-SNR
  plots/E2_Gains_*.png/pdf        - percent reduction vs fixed SCOPE
  plots/E2_RadiusSelection_*.png/pdf - SURE selection r=0,1,2,3
  plots/E2_AggregateRelativeAMSE.png/pdf - across 12 conditions
  Experiment2_Results.mat        - configuration and complete results

RESUMING / MODIFICATIONS
  If interrupted, rerun run_SCOPE_Experiment2: it will resume from the
  last completed batch (via PARTIAL checkpoint), skipping completed
  signal/SNR checkpoints. The calibration is cached separately.
  For methodologically different grids/seed/batch sizes, move or delete
  output checkpoints before rerunning. The default M=1000 is required.

If MATLAB reports a syntax or toolbox error, send the exact error message
and line number so that the issue can be addressed reproducibly.
