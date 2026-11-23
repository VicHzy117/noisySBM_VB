# Simulation replication

R/Rcpp replication of the current main article and simulation Supplement.
Only the paper's simulations are included: no real data, network illustrations,
historical reports, intermediate recovery jobs, Python or precompiled objects.

## Design and methods

Fourteen settings, n=50,100,200,400, and 100 replicates per setting/size:

| Settings | Latent graph |
| --- | --- |
| S1/S2 | Balanced two-community, Q=(0.8,0.1;0.1,0.8), signal mean 2/1 |
| S3/S9 | Three near-duplicate communities, weaker/stronger block-specific means |
| S4/S10/S11/S12 | Balanced bipartite, (cross probability,mean)=(0.65,1)/(0.65,2)/(0.25,2)/(0.25,1) |
| S5/S13 | Uniformly chosen star center, exactly n-1 edges, mean 1/2 |
| S6/S14 | Preferential attachment, mean 1/2 |
| S7/S8 | Balanced two-community, Q=(0.55,0.35;0.35,0.55), mean 2/1 |

Null measurements are N(0,1); signal variance is 1. Exact Q and mean matrices,
community proportions, names and display order live in
`data/generate/settings.R`; `settings.csv` is its human-readable export.
Three-community sizes differ by at most one, with the remainder assigned to
the first groups. PA starts with 0.20n nodes and ER probability 0.50; redraw
until at least floor(0.75*0.20n+0.5) root nodes have positive degree, attach
each new node to that many distinct existing nodes with degree-proportional
weights, and randomly permute node indices after growth.

The paper numbering is used everywhere. In particular S1/S2 use the validated
strong-contrast data (means 2/1), not the superseded 0.65/0.25 connectivity.

- **VB:** structured Bayesian VI, fixed-K candidates selected by maximum ELBO.
- **VEM:** noisySBM 0.1.4, Gauss01, native ICL selection.
- **Greedy:** noisysbmGGM 0.1.2.3 NSBM component, NIG and native merge.
- **SC:** the original local Sun-Cai-style known-null implementation, unchanged.

Candidate K is 1:3 for S1/S2/S7/S8 and 1:4 otherwise, identically bounded for
all three models. VEM uses nbOfTau=2, requested nbOfPointsPerTau=2 (Gauss01
internally sets this to 1), maxNbOfPasses=2, minNbOfPasses=1, one core.
Greedy uses nbOfZ=12, Nbrepet=5, sigma0=1, one core, NIG, native merge.
VB defaults are documented in `method/VB/README.md`. SC retains gamma=0.1 in
its API; the archived implementation's frequency grid does not use gamma.

## Contents

```
data/generate/       authoritative settings, confirmed seeds, generator
data/datasets/       56 compressed RDS files; 100 paired replicates per file
method/VB/          standalone R interface and Rcpp/RcppArmadillo implementation
simulation/         comparator interfaces, metrics, unified runner and Slurm
results/data/       paper replicate metrics, fit diagnostics/times, reference values
results/code/       R-only figures, tables and numerical validation
results/figures/    Figure 2, supplement curves, Figure 3/supplement timing pages
results/tables/     unrounded summaries, two-decimal mean-ARI LaTeX, validation
tests/              executable independent-copy/data/method/plot checks
```

Each dataset contains X and true A as symmetric matrices; Z is retained for
community-truth settings and NULL for star/PA. Data are not duplicated across
methods. Read `data/datasets/S1/n050.rds` as a list of 100 numbered replicates.
RDS files use xz compression: 1.36 GiB total, with the largest file about
76.4 MiB (every file is below 100 MiB). Use Git LFS for the dataset files if
preferred: `git lfs track "data/datasets/**/*.rds"` before adding them to Git.
The actual datasets are included in this repository, not replaced by download
placeholders.

## Dependencies and environment

```r
install.packages(c("Rcpp", "RcppArmadillo", "noisySBM", "ggplot2", "digest"))
```

Use noisySBM 0.1.4 for historical-fit replication (for example with
`remotes::install_version("noisySBM", version="0.1.4")`). VB alone needs only
Rcpp, RcppArmadillo and a C++17 toolchain; no GGM preprocessing dependencies
are required.

Third-party components: the unchanged NSBM backend from noisysbmGGM 0.1.2.3
by Valentin Kilian and Fanny Villers is included in `simulation/Greedy/`,
with its copyright notices and [GPL-2 license](simulation/Greedy/LICENSE)
retained. VEM uses the external noisySBM package. SC is the existing local
Sun-Cai-style implementation, not an identical upstream release. These
comparators are separate from the standalone VB method.

Confirmed cluster environment: R 4.4.2, module `R/4.4.2-gfbf-2024a`, Slurm,
one CPU per task, OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=
BLIS_NUM_THREADS=1. Retained logs report Linux x86_64, Ubuntu 18.04.6,
FlexiBLAS/OpenBLAS, LAPACK 3.12.0, noisySBM 0.1.4, RcppArmadillo 15.2.6-1,
Rcpp 1.1.1-1.1, ggplot2 4.0.3. The exact original compiler build was not
retained; it is not inferred from the module name. Local validation versions
are recorded by `tests/validate.R`. Optimizer results can depend on platform,
R/package versions and floating-point arithmetic; bitwise-identical refits on
arbitrary machines are not promised. Published results are retained separately.

## Seeds and timing

`data/generate/seeds.csv` and its identical RDS version distinguish data,
separate observation-noise streams where applicable, and VB/VEM/Greedy
initialization seeds. SC is deterministic. These seeds were checked against
the original generators/configurations and archived data; no missing seed is
imputed. Within a setting, replicate seeds are deterministic; some paired
settings intentionally share graph/noise streams. Slurm order does not affect
any seed. VB candidate K uses the original seed+1009*K offset internally.

The generator preserves the original normal-draw ordering, including the
separate n=200 null-then-signal stream in the affected SBM settings. Small
macOS/Linux libm rounding differences are allowed only up to 1e-12 in X;
A and applicable Z must match exactly. Archived X is always used for fitting.

Times are **wall-clock elapsed seconds**, not CPU time. Initialization and
candidate search/native merge are included. VB also computes plug-in q-values
inside its fitting call. Backend loading/compilation, data input/output,
metric evaluation, aggregation and plotting are excluded. VEM/Greedy graph
inference is timed separately. Historical fit times span multiple batches,
including older n=200 runs, and are not hardware-matched reruns. Community
diagnostic recovery never replaces the original paper fit times.

## Four examples

Run the examples from this directory (the command-line scripts also resolve
their own root when invoked by absolute path).

**1. Fit VB independently.** Copy `method/VB/` to use it outside this project.

```r
source("method/VB/load.R")
X <- readRDS("data/datasets/S1/n050.rds")[[1]]$X # or your own symmetric X
fit <- vb_fit(X, K_candidates=1:3, seed=1L)
fit$K; fit$labels; fit$tau; fit$parameters; fit$elbo
Ahat <- vb_test(fit, alpha=0.05)$A
```

**2. Read and independently regenerate one replicate; run a method.**

```r
source("data/generate/generate.R")
d <- read_replicate("S1", 50, 1)
g <- generate_replicate("S1", 50, 1)
stopifnot(identical(d$A,g$A), identical(d$Z,g$Z), max(abs(d$X-g$X))<1e-12)
```

```bash
Rscript simulation/run_simulation.R --setting=S1 --n=50 --replicate=1 --method=VB
```

Choose VB, VEM, Greedy or SC. One fit generates results at all five levels
0.01,0.05,0.10,0.15,0.25. A single keyed task file is saved; rerunning the same
task reuses its verified output, never appends duplicate rows. Refit outputs
are separate from the retained paper inputs. To collect reruns:
`Rscript simulation/collect_results.R results/data/reruns results/data/new_summary`.

**3. Submit the complete Slurm array.**

```bash
mkdir -p results/data/reruns/logs
sbatch --array=1-22400%28 simulation/slurm/run_array.slurm
```

The same entry point accepts `--task-id=1`. It maps deterministically to the
explicit setting/n/replicate/method grid, without using the task ID as a seed.

**4. Generate all paper simulation figures and tables without fitting.**

```bash
Rscript results/code/reproduce.R
```

This reads only saved result inputs plus the shared setting configuration.
Figure 2 is `results/figures/TDR_FDR_main_six_settings.pdf`; S7-S14 curves
are `S07_TDR_FDR.png` through `S14_TDR_FDR.png`. Figure 3 is page 1 of
`fitting_time_s1_s14.pdf`; supplement timing results are pages 2-14.
The boxplots use raw-second type-7 quartiles/whiskers before the log transform,
with no horizontal offsets. `results/tables/mean_ari.tex` uses two decimals;
`mean_ari.csv` retains unrounded means. Numerical checks compare every
FDR/TDR point, mean ARI and timing median with the current paper inputs.

## Statistical definitions and validation

For each replicate: TP/FP count correctly/incorrectly discovered edges;
FDP=FP/max(discoveries,1), TDR=TP/(TP+FN). Zero discoveries have FDP=0.
Empirical FDR is mean replicate FDP, **not** pooled FP/discoveries. All five
alphas reuse a fit. VB/VEM/Greedy report K and ARI where true communities
exist; SC has no K/ARI. Candidate K is not silently replaced by the number
of nonempty hard labels. Valid nonconverged fits are not filtered out.

```bash
Rscript tests/validate.R
```

The checks copy VB into temporary storage and fit in a clean R session;
regenerate a replicate in every setting/size cell; run all four methods and
metrics on one paired replicate and test idempotence; then copy only the saved
results, plotting code and setting configuration to another temporary root
and regenerate all figures/tables without datasets or method code. Original
archived data were also regenerated/checked during assembly; the numerical
audit is `data/generate/verification.rds`. No historical reports or recovery
machinery are needed to use this directory.

The delivered checks passed in a fresh local R session: all four methods'
TP/FP matched the retained paper inputs at all five levels for the tested
replicate. `tests/validation_results.csv` records the checks. Verify copied
files with `Rscript tests/check_integrity.R`. Intentional regeneration changes
derived file bytes; renew the manifest with `--write` only after validating
those changes.
