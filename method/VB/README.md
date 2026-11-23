# Structured variational Bayes for noisy networks

Copy **this entire VB directory** anywhere. Only R, Rcpp, RcppArmadillo and a
C++17-capable R toolchain are required. No simulation code/data, comparator,
saved workspace, compiled binary or external project path is used.

```r
install.packages(c("Rcpp", "RcppArmadillo")) # once
source("method/VB/load.R")                  # change to the copied directory
fit <- vb_fit(X, K_candidates = 1:3, seed = 1L)
fixed <- vb_fit(X, K = 2L, seed = 1L)
fit$K; fit$labels; fit$tau; fit$parameters; fit$elbo
graph <- vb_test(fit, alpha = 0.05)
graph$A
```

X must be finite, symmetric, square, with zero diagonal; one measurement per
unordered pair, with known null N(0,1). Compilation is automatic and uses
private temporary storage. `load_vb()` can be called before a timed fit.
The original `fit_nsbm_vi_spectral()` and `fit_nsbm_vi_select_k()` interfaces
remain available. `vb_fit()` additionally presents a compact, consistent result.

Defaults reproduce the paper: Beta shapes 2/2 for all block pairs; NIG
rou=0, kappa=0.01, delta=3, xi=2; allocation gamma=1; initializer epsilon=0.1,
25 k-means starts, two 15% perturbations, 15 screening sweeps, 300 continuation
sweeps, tolerance 1e-6, minimum 3 sweeps. Selection uses maximum fixed-K ELBO,
with no log prior on K and no log(K!) addition. `...` accepts the original
initializer/iteration controls and `hyperparameters` can be overridden.

Edges use the original analytic-mFDR plug-in q-values, with strict q<alpha.
`vb_test()` returns a symmetric adjacency matrix plus q-values, local l-values
and their explicit i-then-j unordered pair indices. `parameters$sigma2` is
variance, not standard deviation. Labels are row-wise argmax of tau; candidate
K need not equal the number of occupied hard labels.

This is one implementation, not a collection of old versions. A single explicit
`count_pmf_rebuild` switch surrounds the existing positive-convolution numerical
fallback: TRUE (default) enables the validated stability protection; FALSE
reproduces batches that used polynomial division without that fallback.
All CAVI equations, ELBO terms, initialization routines and graph-testing rules
are unchanged. Simulation seeds select the historically used switch per cell.
No ARI, plotting, Slurm or competing method is included here.
