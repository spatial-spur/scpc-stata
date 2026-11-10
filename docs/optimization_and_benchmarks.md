# SCPC numerical methods and validation

The working covariance is `Sigma(c)[i,j] = exp(-c * distance[i,j])`, calibrated to the requested average off-diagonal correlation. A safeguarded Newton solve uses a data-scaled bracket and checks attainable correlations, including co-located observations.

Rejection probabilities use an adaptive integral after `u = sin(theta)`. Critical values return the upper root bracket and are checked against Student-t special cases. The conditional covariance retains `W'W` even when the weights are not orthonormal.

Exact mode stores the distance matrix and one temporary covariance at a time; it does not retain the full covariance grid. The automatic exact/approximate cutoff is 4,500 locations. Approximate mode uses sampling without replacement, the same geographic metric in every kernel, and a single normalization/permutation across basis-size adaptation.

Run `do tests/test_numerical.do` from the repository root in licensed Stata. This chat's Linux environment has no Stata executable, so Stata execution and benchmark timings remain unverified. Run `benchmarks/benchmark_scpc.do` to collect fresh measurements; tighter numerical tolerances can change runtime and upstream output.

A finite covariance grid and sampled large-N covariance are numerical approximations. Root precision alone does not prove size control over the entire continuous covariance class or eliminate Monte Carlo error.
