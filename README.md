# `scpc`: High-Performance Spatial Correlation Robust Inference (v2.1.0)

[![GitHub Release](https://img.shields.io/badge/version-2.1.0-blue.svg)](https://github.com/spatial-spur/scpc-stata)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/Stata-13%2B-brightgreen.svg)]()

This repository contains the high-performance, pure-Mata modernization of the **SCPC (Spatial Correlation Principal Components)** Stata package developed by **Ulrich K. Müller and Mark W. Watson**, maintained by **Sascha O. Becker, P. David Boll, and Hans-Joachim Voth** as part of the [`spatial-spur`](https://github.com/spatial-spur) econometrics ecosystem.

Version 2.1.0 delivers a **2.0× to 2.5× end-to-end performance improvement** while remaining **100% pure Mata** (no external compilers, no C plugins, no dependencies).

---

## Performance Overview

| Sample Size ($N$) | Mode | Upstream Baseline (`ukmueller/SCPC`) | `scpc` v2.1.0 (Pure Mata) | Speedup | Numerical Precision |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **$N = 200$** | Exact | 0.179 s | **0.094 s** | **1.90×** | Exact ($0.000$) |
| **$N = 722$** | Exact | 0.580 s | **0.288 s** | **2.01×** | Exact ($0.000$) |
| **$N = 1,500$** | Exact | 2.424 s | **1.146 s** | **2.12×** | Exact ($0.000$) |
| **$N = 3,000$** | Exact | 11.718 s | **5.559 s** | **2.11×** | Exact ($0.000$) |
| **$N = 5,000$** | Approx (Large-$N$) | 40.567 s | **21.161 s** | **1.92×** | Exact ($0.000$) |
| **$N = 10,000$** | Approx (Large-$N$) | 42.871 s | **20.652 s** | **2.08×** | Exact ($0.000$) |

See [**`docs/optimization_and_benchmarks.md`**](docs/optimization_and_benchmarks.md) for full algorithmic breakdowns and mathematical details.

---

## Installation

Install the package directly into Stata:

```stata
// Uninstall previous versions if present
cap ado uninstall scpc

// Install v2.1.0 from spatial-spur
net install scpc, from("https://raw.githubusercontent.com/spatial-spur/scpc-stata/dev-perf-fixes/src")
```

---

## Quick Start

`scpc` runs as a postestimation command following `regress`, `ivregress`, `areg`, `logit`, or `probit` with `robust` or `cluster()` options. Spatial coordinates must be stored in variables beginning with `s_` (e.g. `s_1`, `s_2`):

```stata
sysuse auto, clear

// Define spatial coordinates
gen s_1 = runiform() * 10
gen s_2 = runiform() * 10

// Estimate regression model
regress price weight length foreign, robust

// Compute SCPC inference and critical values
scpc, cvs
```

### Spherical (Latitude / Longitude) Coordinates:
```stata
// s_1 is latitude in degrees, s_2 is longitude in degrees
scpc, latlong
```

### Explicit Computation Modes:
```stata
// Force exact spectral decomposition
scpc, exact

// Force Large-N Nyström approximation
scpc, approx

// Legacy PRNG mode for exact backward replication of 2021 results
scpc, approx legacyprng
```

---

## What's New in Version 2.1.0

1. **64-Point Gauss-Legendre Quadrature & Brent's Root Finder ($10^{-8}$):** Replaced coarse 20-point bisection with adaptive Dekker-Brent root-finding, eliminating step-like jitter in critical values.
2. **Unconditional Critical Value Caching:** Precomputes and caches the theoretical lower bound $c_{\text{uncond}}$, eliminating redundant calculations across regressors.
3. **Vectorized PRNG in Mata:** Vectorized Monte Carlo index draws in Large-$N$ Nyström approximation mode.
4. **Zero Macro Serialization Overhead:** Streamlined data flow between Mata and Stata memory.
5. **Configurable Thresholds:** Explicit `exact`, `approx`, and `threshold()` options give users full control over when the Large-$N$ approximation triggers.

---

## References

- **Müller, Ulrich K., and Mark W. Watson.** 2021. "Spatial Correlation Robust Inference." *Econometrica* 89 (6): 2901–2927.
- **Müller, Ulrich K., and Mark W. Watson.** 2023. "Spatial Correlation Robust Inference in Linear Regression and Panel Models." *Journal of Business & Economic Statistics* 41 (4): 1050–1064.
- **Becker, Sascha O., P. David Boll, and Hans-Joachim Voth.** 2026. "Spatial Regressions and Spatial Unit Roots in Stata." *The Stata Journal* (forthcoming).
