{smcl}
{* *! version 1.1.0 Jan 2022}{...}
{title:Title}

{pstd}
{hi:scpc} {hline 2} Spatial Correlation Robust Inference via SCPC.


{title:Syntax}

{p 8 16 2}
{cmd:scpc} [{cmd:,} {it:options}]
{p_end}

{synoptset 26}{...}
{synopthdr}
{synoptline}
{synopt :{opt avc(#)}}overrides default value of 0.03 for the maximal average pairwise correlation (must be between 0.001 and 0.99){p_end}
{synopt :{opt uncond}}overrides default of computing critical values conditional on regressors to compute "unconditional" SCPC inference as described in {help scpc##mainpaper:Müller and Watson (2021)}{p_end}
{synopt :{opt latlong}}if present, expects s_1 to contain latitude and s_2 to be longitude and uses implied great-circle distance; if not present, computes Euclidian distance between coordinates defines by s_* variables{p_end}
{synopt :{opt k(#)}}overrides default value of 10 for the number of coefficients SCPC t-statistics are computed for{p_end}
{synopt :{opt cvs}}prints a table of two-sided critical values of SCPC t-statistic of level 32%, 10%, 5% and 1%{p_end}
{synopt :{opt backend(auto|rust|fallback)}}choose the computational backend. {opt auto} uses the matching Rust plugin when available and Mata otherwise. {opt rust} requires the matching plugin. Set the default with {cmd:SPATIAL_BACKEND} as a Stata global or environment variable.{p_end}
{synopt :{opt nthreads(#)}}set the number of Rust threads. {opt nthreads(1)} runs Rust with one thread.{p_end}
{synopt :{opt exact}}use exact SCPC.{p_end}
{synopt :{opt approx}}use the Nyström approximation.{p_end}
{synopt :{opt seed(#)}}set the seed for the Nyström approximation. Default is {opt seed(1)}.{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
This Stata package implements the Spatial Correlation Principal Components (SCPC) method described in {help scpc##mainpaper:Müller and Watson (2022, 2023)} for the construction of confidence intervals that account for many forms of spatial correlation.
The {it:scpc} command expects the locations of the observations to be stored in the variables s_*. If the option {it:latlong} is present, then s_1 is interpreted as the latitude, and s_2 as the longitude of the observation, and distances between observations are computed with the great-circle formula. If the option {it:latlong} is not present, then all {it:p} variables beginning with s_ are interpreted as points in {it:p}-dimensional space, and distances are computed with the Euclidian norm. 
scpc is implemented as a postestimation command that can be used after the unweighted Stata commands {manhelp regress R:regress}, {manhelp ivregress R:ivregress}, {manhelp areg R:areg}, {manhelp logit R:logit} or {manhelp probit R:probit},
as long as these are used with the standard error option {it:robust} or {it:cluster} (see {manhelp vce_option R:vce option}). If the estimation uses the {it:cluster} option,
then {it:scpc} corrects for spatial correlations between clusters, assuming that all observations within a cluster share the same location.
By default, {cmd:scpc} computes critical values conditional on the realized regressors, as described in Müller and Watson (2023). Conditional inference supports {manhelp regress R:regress} and {cmd:ivregress 2sls}. Use {it:uncond} for unconditional critical values after {manhelp areg R:areg}, {manhelp logit R:logit}, or {manhelp probit R:probit}. After {cmd:ivregress liml}, use {it:uncond} with the Mata backend.
By default, {it:scpc} conducts spatially robust inference under the assumption that the largest average pairwise correlation between the observations / clusters is no larger than 0.03. This default can be overridden by the
{it:avc} option. Note that computation times increase for smaller values of {it:avc}. To make the algorithm faster, SCPC inference is computed only for the first 10 coefficients. This can be changed by the {it:k} option. 
By default, SCPC uses the Nyström approximation for 4,500 or more locations. A sampled covariance that fails the positive semidefinite check produces an error. {opt exact} uses the full covariance grid.

 

{marker options}{...}
{title:Options}

{phang}
{opt avc(#)} overrides default value of 0.03 for the maximal average pairwise correlation; see {help scpc##mainpaper:Müller and Watson (2021)} for details.

{phang}
{opt uncond} overrides the default of computing conditional critical values in regression; see {help scpc##mainpaper:Müller and Watson (2021)} for details.

{phang}
{opt latlong} if present, expects s_1 to contain latitude and s_2 to be longitude and uses implied great-circle distance; if not present, computes Euclidian distance between coordinates defines by variables beginning with "s_".

{phang} 
{opt cvs} prints a table of two-sided critical values of SCPC t-statistic of level 32%, 10%, 5% and 1%.

{phang} 
{opt k(#)} overrides default value of 10 for the number of coefficients for which SCPC inference is computed.

{phang}
{opt backend(auto|rust|fallback)} selects the computational backend. {opt auto} uses the matching Rust plugin when eligible and falls back to Mata otherwise. {opt rust} requires the matching plugin. {opt fallback} runs in Mata. The Stata global or environment variable {cmd:SPATIAL_BACKEND} sets the default when present. Restart Stata after replacing an incompatible plugin that was already loaded.

{phang}
{opt nthreads(#)} sets the number of threads used by the Rust backend.

{phang}
{opt exact} forces exact SCPC regardless of sample size.

{phang}
{opt approx} uses the Nyström approximation. An indefinite sampled covariance raises an error rather than falling back.

{phang}
{opt seed(#)} sets the nonnegative 32-bit integer seed used by the Nyström approximation (default {cmd:1}).

{marker examples}{...}
{title:Examples}

{phang}{cmd:. sysuse auto}{p_end}
{phang}{cmd:. gen s_1=rnormal(0,1)}{p_end}
{phang}{cmd:. gen s_2=rnormal(0,1)}{p_end}
{phang}{cmd:. regress mpg weight length, robust}{p_end}
{phang}{cmd:. scpc}{p_end}
{phang}{cmd:. scpc ,avc(0.03)}{space 6}(equivalent to above command){p_end}
{phang}{cmd:. scpc ,k(4)}(equivalent to above command){p_end}
{phang}{cmd:. scpc ,avc(0.01) cvs k(1)}{p_end}
{phang}{cmd:. gen clust=round(rnormal(0,10),1)}{p_end}
{phang}{cmd:. regress mpg weight length, cluster(clust)}{p_end}
{phang}{cmd:. scpc}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:scpc} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(scpcstats)}}Matrix of SCPC inference results (same as printed table){p_end}
{synopt:{cmd:e(scpccvs)}}Matrix of SCPC critical values (same as printed table, only computed under {it:cvs} option){p_end}

{p2col 5 15 19 2: Macros}{p_end}
{synopt:{cmd:e(backend)}}computational backend used{p_end}
{synopt:{cmd:e(backend_reason)}}reason for using Mata, when applicable{p_end}
{synopt:{cmd:e(method)}}exact or approximate SCPC{p_end}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:e(q)}}number of selected low-frequency spatial cosine/eigenvector averages{p_end}
{synopt:{cmd:e(c0)}}calibrated kernel scale parameter corresponding to {opt avc}{p_end}
{synopt:{cmd:e(cv)}}95% SCPC critical value used for confidence intervals{p_end}
{synopt:{cmd:e(seed)}}seed recorded for the SCPC run{p_end}
{synopt:{cmd:e(tied)}}indicator of a tied eigenspace reported by Rust{p_end}
{p2colreset}{...}

{pstd}
Additionally, {cmd:scpc} preserves all macros and scalars of the estimated model in memory.



{marker authors}{...}
{title:Author}

{pstd}
Ulrich K. Müller{break}
Princeton University{break}


{marker disclaimer}{...}
{title:Disclaimer}

{pstd}
This software is provided "as is", without warranty of any kind.
If you have suggestions or want to report problems, please create a new issue in the {browse "https://github.com/ukmueller/SCPC/issues":project repository} or contact the project maintainer.


{marker references}{...}
{title:References}

{marker mainpaper}{...}
{phang}Müller, Ulrich K. and Mark W. Watson "Spatial Correlation Robust Inference", Econometrica 90 (2022), 2901–2935. {browse "https://www.princeton.edu/~umueller/SHAR.pdf"}.

{marker followpaper}{...}
{phang}Müller, Ulrich K. and Mark W. Watson "Spatial Correlation Robust Inference in Linear Regression and Panel Models", Journal of Business and Economic Statistics 41 (2023), 1050–1064. {browse "https://www.princeton.edu/~umueller/SpatialRegression.pdf"}.


