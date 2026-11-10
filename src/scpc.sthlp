{smcl}
{* *! version 2.1.0  10nov2026}{...}
{vieweralsosee "[R] regress" "help regress"}{...}
{vieweralsosee "[R] ivregress" "help ivregress"}{...}
{viewerjumpto "Syntax" "scpc##syntax"}{...}
{viewerjumpto "Description" "scpc##description"}{...}
{viewerjumpto "Options" "scpc##options"}{...}
{viewerjumpto "Performance Optimizations" "scpc##optimizations"}{...}
{viewerjumpto "Stored results" "scpc##results"}{...}
{viewerjumpto "Authors and References" "scpc##references"}{...}
{title:Title}

{p2colset 5 16 18 2}{...}
{p2col :{cmd:scpc} {hline 2}}High-Performance Pure Mata Spatial Correlation Principal Components (SCPC) Inference{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:scpc} [{cmd:,} {it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt :{opt cvs}}compute and display two-sided conditional critical values at 32%, 10%, 5%, and 1% significance levels{p_end}
{synopt :{opt k(#)}}number of coefficients to report SCPC inference for (default is 10 or total number of regressors){p_end}
{synopt :{opt avc(#)}}target average pairwise spatial correlation; default is {cmd:avc(0.03)}{p_end}
{synopt :{opt latlong}}treat location coordinates {cmd:s_1} (latitude) and {cmd:s_2} (longitude) as coordinates on a sphere{p_end}
{synopt :{opt uncond}}report only unconditionally valid critical values; allows commands beyond {cmd:regress} and {cmd:ivregress}{p_end}

{syntab:Performance & Numerical Options (v2.1 Pure Mata)}
{synopt :{opt exact}}force exact spatial spectral decomposition; eliminates all Monte Carlo sampling variance{p_end}
{synopt :{opt approx}}force Large-N Nyström approximation{p_end}
{synopt :{opt threshold(#)}}sample size threshold below which exact computation is automatically chosen; default is {cmd:threshold(4500)}{p_end}
{synopt :{opt cgridfac(#)}}grid expansion multiplier for spatial scale grid {it:c}; default is {cmd:cgridfac(1.2)}{p_end}
{synopt :{opt legacyprng}}use legacy 32-bit linear congruential generator and 40-point Gauss-Legendre quadrature for a legacy deterministic stream{p_end}
{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd}
{cmd:scpc} conducts spatial correlation robust inference following the Spatial Correlation Principal Components (SCPC) method
of Müller and Watson (2022, 2023). Location coordinates must be stored in variables prefixed with {cmd:s_} ({cmd:s_1}, {cmd:s_2}, ...).
Observations belonging to the same cluster (if {cmd:vce(cluster ...)} was specified) must have identical coordinates.

{pstd}
Version 2.1 is implemented in {bf:100% pure Mata} with zero external binary or compiled plugin dependencies, so no separate plugin binary is required; platform-specific execution still needs testing.

{marker optimizations}{...}
{title:Performance & Mathematical Optimizations}

{pstd}
Numerical calibration uses adaptive angular quadrature, an upper critical-value bracket,
and a safeguarded spatial-scale solver. Exact mode retains one covariance matrix at a time;
it does not cache the full grid. Approximate mode samples without replacement and retains
{it:W'W} for nonorthonormal conditional weights. Numerical tolerances do not account for
finite-grid or Monte Carlo approximation error. Historical benchmark claims and exact
upstream parity are not correctness guarantees.

{marker options}{...}
{title:Options}

{phang}
{opt cvs} computes and reports the two-sided critical values of the SCPC t-tests at the 32%, 10%, 5%, and 1% levels, stored in {cmd:e(scpccvs)}.

{phang}
{opt k(#)} specifies the number of regressor coefficients for which SCPC inference is computed. The default is 10 or all coefficients, whichever is smaller.

{phang}
{opt avc(#)} specifies the target average pairwise spatial correlation. The default is {cmd:avc(0.03)}. Must be strictly between 0.0001 and 0.99.

{phang}
{opt latlong} indicates that {cmd:s_1} contains latitude and {cmd:s_2} contains longitude in degrees. Geodesic distances are computed on the unit sphere using a numerically clamped 3D chord formula.

{phang}
{opt exact} forces the exact eigenvalue decomposition algorithm regardless of sample size {it:N}.

{phang}
{opt approx} forces the Large-{it:N} Nyström approximation regardless of sample size {it:N}.

{phang}
{opt threshold(#)} sets the sample size {it:N} cutoff between exact and approximate computation. Default is 4,500.

{phang}
{opt legacyprng} activates the legacy 32-bit linear congruential generator (LCG) and 40-point Gauss-Legendre quadrature for a legacy deterministic stream; corrected calculations need not reproduce historical output.

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:scpc} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(scpcstats)}}inference table containing columns: Coef, Std_Err, t, P>|t|, 95% Conf Lower, 95% Conf Upper{p_end}
{synopt:{cmd:e(scpccvs)}}two-sided critical values at 32%, 10%, 5%, and 1% levels (if {cmd:cvs} option specified){p_end}
{p2colreset}{...}

{marker references}{...}
{title:References}

{phang}
Müller, Ulrich K., and Mark W. Watson. 2022. "Spatial Correlation Robust Inference."
{it:Econometrica} 90 (6): 2901–2935.
{p_end}
