SCOPE MATLAB -- Rounded-B_N calibration (Section 4)
===================================================

Files:
  SCOPE_Exp123_BN.m       Updated complete Experiments I, II, III
  SCOPE_compare_BN.m      Optional old-versus-new AMSE summary comparison

Section 4 calibration (at each detail level j):
  a_{F,r,k,j} = R_{F,r,k}^{-1}(B_N / 2^j)
  lambda_{j,r}(F,k,b) = b * a_{F,r,k,j} / sigmahat
  For fixed SCOPE/N-SCOPE, set k=5 and b=1.
  For SURE, k={3,5,7}, r={0,1,2,3}, b={.65,.80,1,1.20,1.50}.

Predeclared budgets and coarsest shrunk levels:
       N        B_N     j0
     128       1.75      3
     512       3.5       3
    1024       5         3
    2048       7.5       3
    4096      11.5       4

The normalizing risk inverse is computed once per generator/radius/k/level,
by the original one-dimensional chi-square integral and fzero solver.
No parameter is optimized against the true signal.

Run from the MATLAB Command Window with both .m files in MATLAB's path:

  SCOPE_Exp123_BN('smoke','all')

Or, as an initial full-size test with a single N and a single root SNR:

  SCOPE_Exp123_BN('full','Exp1',1024,5,1000)
  SCOPE_Exp123_BN('full','Exp2',1024,5,1000)
  SCOPE_Exp123_BN('full','Exp3',1024,5,1000)

For all three sample sizes and all root SNRs:

  SCOPE_Exp123_BN('full','all',[512 1024 2048],[3 5 7],1000)

Results are written into a NEW folder next to this file:

  SCOPE_Exp123_results_BNrounded

The original SCOPE_Exp123_results folder is not modified, and original
checkpoints are not reused. New code uses the same deterministic random
seeds as the original, enabling matched-condition comparisons.

If you already have matching original results and completed the new run:

  SCOPE_compare_BN

This writes SCOPE_Exp123_results_BNrounded/BNrounded_versus_original_AMSE.csv
with original AMSE, new AMSE and percentage changes by N/signal/SNR/method.

Important: the source code's SNR is ROOT SNR, std(f,1)/sigma; it is NOT
power SNR E(f^2)/sigma^2. The revised code intentionally retains the
original definition, signal scaling, wavelet filters, seeds and methods.

Validation performed in this environment: source-level structural
checks and independent numerical evaluation of selected Gaussian
null-risk calibration roots. The full MATLAB program was NOT executed
here because MATLAB is not available in the sandbox.
