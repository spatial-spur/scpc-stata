# `scpc`: Spatial Correlation Robust Inference

This Stata package implements the C-SCPC method described in [Müller and Watson (2021)](http://www.princeton.edu/~umueller/SCPC.pdf) and SCPC method [Müller and Watson (2022)](http://www.princeton.edu/~umueller/SpatialRegression.pdf) for the construction of confidence intervals that account for many forms of spatial correlation.
It is implemented as a postestimation command that can be used after the STATA commands "regress", "ivregress", "areg" , "logit" or "probit" as long as these are used with the standard error option "robust" or "cluster". If "cluster" is chosen, then the method assumes that all observations in a cluster are at the same spatial location, and corrects for potential spatial correlation between clusters. 

## Install a local checkout

Run this in Stata with the absolute path to this checkout.

```stata
net install scpc, replace from("/path/to/spatial-spur/scpc-stata/src")
```

The checked-in development plugin targets macOS arm64 and was tested with
StataNow 19.5. Build the plugin for your machine from a checkout of the shared
Rust workspace. From the parent directory containing both packages, run
`make -C spatial-rust/crates/spatial-stata` after fetching the locked Cargo
dependencies. The Makefile has macOS and Linux targets. Restart Stata after
replacing a plugin. Use `scpc, backend(fallback)` to select the Mata
implementation directly. The default `backend(auto)` uses Mata if the matching
plugin cannot load. Check `e(backend)` and `e(backend_reason)` after the call.

Conditional inference supports unweighted `regress` and `ivregress 2sls`
fits. For `ivregress liml`, use `scpc, uncond backend(fallback)`.

# References

Müller, Ulrich K and Mark W. Watson (2021). "Spatial Correlation Robust Inference". Working Paper. https://www.princeton.edu/~umueller/SHAR.pdf.

Müller, Ulrich K and Mark W. Watson (2022). "Spatial Correlation Robust Inference in Linear Regression and Panel Models". Working Paper. https://www.princeton.edu/~umueller/SptialRegression.pdf.
