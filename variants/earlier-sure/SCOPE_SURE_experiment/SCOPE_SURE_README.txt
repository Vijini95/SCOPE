SCOPE plug-in SURE Monte Carlo experiment

Executed results: Python 3 / NumPy / SciPy implementation, not MATLAB.
SCOPE_SURE_test.m is a MATLAB translation supplied for reproduction; not executed here.
Run MATLAB: [Summary, Replications] = SCOPE_SURE_test;
An external Wavmat.m is used if present; otherwise its matrix algorithm is embedded.
No built-in wavelet transform is used. Python uses the same matrix algorithm.

N=1024, R=100, seed=20261024, power SNR=mean(f.^2)/sigma^2=5,
noise iid N(0,1). All distributions and signals share the same noise matrix.
MATLAB and Python random streams differ, so their simulation values will differ.
Filters: Haar, six-tap Daubechies (db3), Sym4, Sym4.
Wavmat depth=7 and shift=2; eight coarse coefficients are untouched.
Sigma estimate: median(abs(finest coefficients))/0.6744897501960817.
The HeaviSine formula follows the earlier runs: 4*sin(4*pi*t)-sign(t-.3)-sign(.72-t).

lambda_j = q_F(k)*c*2^(-alpha*j)/sigmahat.
q_F(k) is the inverse of F*(x)=2F(x)-1 at 2^(-1/k).
Grid fixed before selection: alpha=0:.05:1.5; log2(c)=-3:.05:10.
Mode k5 fixes k=5. Mode select_k searches k in {1,3,5,10}.
Each noisy replication is tuned separately for each distribution.
Selection uses ONLY noisy coefficients and estimated sigma, never truth.
Equal minima choose the first candidate; candidate traversal order can differ by
language, potentially producing different parameter choices for exact ties.

For fixed candidates, the criterion is sum((delta-d)^2 + 2*sigmahat^2*delta' - sigmahat^2).
The constant coarse-level contribution is omitted in minimization.
The scale estimate is held fixed when calculating derivatives, so this is plug-in
SURE, not exactly unbiased SURE for the full data-dependent procedure.
Risk sums for every grid candidate are computed exactly by reusing level risk
curves on an aligned log2 scale grid; no interpolation is used.
No continuous refinement, oracle parameter tuning, or tuning across simulations.
True signals enter only the construction of observations and final MSE evaluation.
AMSE is mean reconstruction MSE across 100 runs; MCSE is std(MSE)/sqrt(100).
BoundaryCount counts runs whose selected alpha or log2(c) hits a grid endpoint.
Selections on boundaries are legitimate for this finite candidate family but
should not be described as unrestricted SURE minima.

Validation: orthogonality checked below 1e-10; analytic derivatives checked against
finite differences for all six distributions and all four k values, away from
Uniform's derivative discontinuity. Final MSE uses reconstructed signals.

Contents:
scope_sure.py: executed Python source
SCOPE_SURE_test.m: MATLAB reproduction code
SCOPE_SURE_summary.csv: 48 rows (4 signals x 6 distributions x 2 modes)
SCOPE_SURE_replications.csv: 4800 rows including selected alpha,c,k and sigmahat
