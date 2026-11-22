# `scpc`: Spatial Correlation Robust Inference

This Stata package implements the C-SCPC method described in [Müller and Watson (2021)](http://www.princeton.edu/~umueller/SCPC.pdf) and SCPC method [Müller and Watson (2022)](http://www.princeton.edu/~umueller/SpatialRegression.pdf) for the construction of confidence intervals that account for many forms of spatial correlation.
It is implemented as a postestimation command that can be used after the STATA commands "regress", "ivregress", "areg" , "logit" or "probit" as long as these are used with the standard error option "robust" or "cluster". If "cluster" is chosen, then the method assumes that all observations in a cluster are at the same spatial location, and corrects for potential spatial correlation between clusters. 

## Install a local checkout

This development checkout requires Stata 16 or later. Run this in Stata with
the absolute path to the checkout.

```stata
net install scpc, replace from("/path/to/spatial-spur/scpc-stata/src")
```

The development checkout bundles a plugin for Apple Silicon, tested with
StataNow 19.5. `net install` installs the plugin on Apple Silicon and the
Mata implementation on all platforms. For a source build, run
`make -C spatial-rust/crates/spatial-stata` from the shared workspace after
fetching the locked Cargo dependencies. The Makefile supports macOS and Linux.
Use `adopath ++` with the `scpc-stata/src` directory to load your locally built
plugin from the checkout. Restart Stata after replacing a plugin. Use
`scpc, backend(fallback)` to select the Mata implementation directly.
The default `backend(auto)` uses Mata if the matching
plugin cannot load. Check `e(backend)` and `e(backend_reason)` after the call.

Conditional inference supports unweighted `regress` and `ivregress 2sls`
fits. For 2SLS, `k()` counts estimable coefficients. Conditional 2SLS with
active factor or time-series terms uses Mata. For `ivregress liml`,
use `scpc, uncond backend(fallback)`.

# References

Müller, Ulrich K and Mark W. Watson (2021). "Spatial Correlation Robust Inference". Working Paper. https://www.princeton.edu/~umueller/SHAR.pdf.

Müller, Ulrich K and Mark W. Watson (2022). "Spatial Correlation Robust Inference in Linear Regression and Panel Models". Working Paper. https://www.princeton.edu/~umueller/SptialRegression.pdf.
