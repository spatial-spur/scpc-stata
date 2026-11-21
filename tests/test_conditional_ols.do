* Run from the shared workspace root or pass its path as the first argument.
clear all
set more off
version 14
args repo_root
if `"`repo_root'"' == "" local repo_root "."
capture confirm file "`repo_root'/scpc-stata/src/scpc.ado"
if _rc == 0 adopath ++ "`repo_root'/scpc-stata/src"
which scpc
quietly findfile scpc.ado
local scpc_source `"`r(fn)'"'

input s_1 x
0 -1
1 -4
2  2
4  2
6  1
end
generate double y = x + (s_1 == 0)

* Reference cutoffs are from Gaussian quadrature under the §5.1 null.
* The tolerance allows for SCPC's 0.001 bisection interval.
quietly regress y x, vce(robust)
quietly scpc, avc(.5) k(1) cvs backend(fallback)
matrix cv = e(scpccvs)
assert abs(cv[1,3] - 17.625169) < .002

quietly regress y x, vce(robust) noconstant
quietly scpc, avc(.5) k(1) cvs backend(fallback)
matrix cv = e(scpccvs)
assert abs(cv[1,3] - 9.973326) < .002

quietly regress y x, vce(robust)
quietly scpc, avc(.5) k(1) cvs backend(rust)
matrix cv = e(scpccvs)
assert abs(cv[1,3] - 17.625169) < .002

quietly regress y x, vce(robust) noconstant
quietly scpc, avc(.5) k(1) cvs backend(rust)
matrix cv = e(scpccvs)
assert abs(cv[1,3] - 9.973326) < .002

* Oblique IV exercises both projection paths and no-exogenous score extraction.
generate double z = x + .01*s_1
foreach cons in constant noconstant {
    local opt = cond("`cons'" == "noconstant", "noconstant", "")
    local iv_ref = cond("`cons'" == "noconstant", 9.798657, 17.432446)
    quietly ivregress 2sls y (x=z), vce(robust) `opt'
    quietly scpc, avc(.5) k(1) cvs backend(fallback)
    matrix fallback_cv = e(scpccvs)
    assert abs(fallback_cv[1,3] - `iv_ref') < .002
    quietly ivregress 2sls y (x=z), vce(robust) `opt'
    quietly scpc, avc(.5) k(1) cvs backend(rust)
    assert "`e(backend)'" == "rust"
    matrix native_cv = e(scpccvs)
    assert abs(native_cv[1,3] - `iv_ref') < .002
    mata: assert(max(abs(st_matrix("fallback_cv") - st_matrix("native_cv"))) < 1e-8)
}

* Changing units of a regressor must not change the 2SLS cutoff.
generate double scaled_x = 1e-8 * x
quietly ivregress 2sls y (scaled_x = z), vce(robust)
foreach backend in rust fallback {
    quietly scpc, avc(.5) k(1) cvs backend(`backend')
    matrix cv = e(scpccvs)
    assert abs(cv[1,3] - 17.432446) < .002
}

* Check public parity at a condition number Stata's 2SLS fit still accepts.
preserve
clear
set obs 8
generate double s_1 = _n - 1
generate double w = (_n == 1) + 0.1 * (_n == 4)
generate double z2 = (_n == 2) + 0.1 * (_n == 5)
generate double x2 = w + 1e-4 * z2 + 0.5 * (_n == 3)
generate double y = w + 2 * x2 + 0.3 * sin(_n * 1.3)
quietly ivregress 2sls y w (x2 = z2), vce(robust) noconstant
assert e(rank) == 2
quietly scpc, avc(.5) k(2) cvs backend(rust)
assert "`e(backend)'" == "rust"
matrix ill_rs = e(scpcstats)
matrix ill_rc = e(scpccvs)
quietly scpc, avc(.5) k(2) cvs backend(fallback)
assert "`e(backend)'" == "fallback"
matrix ill_fs = e(scpcstats)
matrix ill_fc = e(scpccvs)
mata: assert(max(abs(st_matrix("ill_rc") - st_matrix("ill_fc"))) < 1e-10)
mata: assert(max(abs(st_matrix("ill_rs") - st_matrix("ill_fs")) :/ (1 :+ abs(st_matrix("ill_fs")))) < 5e-8)
restore

* Zero unclustered score sign follows >= 0 -> +1 on both Rust and Mata backends.
preserve
clear
input double s_1 double x double y
0 -2 -3
1 -1 -1
2  0  2
3  1  1
4  2  1
end
quietly regress y x, vce(robust)
quietly scpc, avc(.5) k(1) cvs backend(rust)
assert "`e(backend)'" == "rust"
matrix zs_r = e(scpcstats)
matrix zs_rcvs = e(scpccvs)
quietly scpc, avc(.5) k(1) cvs backend(fallback)
assert "`e(backend)'" == "fallback"
matrix zs_f = e(scpcstats)
matrix zs_fcvs = e(scpccvs)
mata: assert(max(abs(st_matrix("zs_r") - st_matrix("zs_f"))) < 1e-10)
mata: assert(max(abs(st_matrix("zs_rcvs") - st_matrix("zs_fcvs"))) < 1e-10)
restore

expand 2
bysort s_1: replace x = x + cond(_n == 1, -.25, .25)
replace y = x + (s_1 == 0)
replace z = x + .01*s_1
quietly ivregress 2sls y (x=z), vce(cluster s_1)
quietly scpc, avc(.5) k(1) cvs backend(fallback)
matrix fallback_cv = e(scpccvs)
assert abs(fallback_cv[1,3] - 17.178540) < .002
quietly ivregress 2sls y (x=z), vce(cluster s_1)
quietly scpc, avc(.5) k(1) cvs backend(rust)
assert "`e(backend)'" == "rust"
matrix native_cv = e(scpccvs)
assert abs(native_cv[1,3] - 17.178540) < .002
mata: assert(max(abs(st_matrix("fallback_cv") - st_matrix("native_cv"))) < 1e-8)

quietly regress y x, vce(cluster s_1)
quietly scpc, avc(.5) k(1) cvs backend(fallback)
matrix cv = e(scpccvs)
assert abs(cv[1,3] - 17.380977) < .002

quietly regress y x, vce(cluster s_1) noconstant
quietly scpc, avc(.5) k(1) cvs backend(fallback)
matrix cv = e(scpccvs)
assert abs(cv[1,3] - 9.902023) < .002

quietly regress y x, vce(cluster s_1)
quietly scpc, avc(.5) k(1) cvs backend(rust)
assert "`e(backend)'" == "rust"
matrix cv = e(scpccvs)
assert abs(cv[1,3] - 17.380977) < .002

quietly regress y x, vce(cluster s_1) noconstant
quietly scpc, avc(.5) k(1) cvs backend(rust)
assert "`e(backend)'" == "rust"
matrix cv = e(scpccvs)
assert abs(cv[1,3] - 9.902023) < .002

* Invalid sampled covariance is an error even in auto (uncond and cond). Exact stays available.
clear
set obs 8
generate double s_1 = _n - 1
generate double s_2 = sin(s_1)
generate double x = cos(2 * s_1)
generate double y = cos(s_1) + 0.5 * x
quietly regress y, robust
foreach backend in rust auto fallback {
    capture noisily scpc, approx avc(.95) uncond backend(`backend')
    assert _rc == 198
    quietly scpc, exact avc(.95) uncond backend(`backend')
    assert "`e(method)'" == "exact"
}
quietly regress y x, robust
foreach backend in rust auto fallback {
    capture noisily scpc, approx avc(.95) backend(`backend')
    assert _rc == 198
}

* Coordinates and cluster representatives come from the fitted sample only.
clear
set seed 23
set obs 30
generate cl = ceil(_n / 3)
generate byte selected = mod(_n - 1, 3) != 0
generate double s_1 = cond(selected, cl, -9999)
generate double s_2 = cos(cl)
generate double x = rnormal()
generate double y = x + rnormal()
quietly regress y x if selected, vce(cluster cl)
quietly scpc, latlong cvs backend(rust)
matrix native_cv = e(scpccvs)
quietly scpc, latlong cvs backend(fallback)
matrix fallback_cv = e(scpccvs)
mata: assert(max(abs(st_matrix("fallback_cv") - st_matrix("native_cv"))) < 1e-8)
replace s_2 = 170 in 2
foreach backend in rust auto fallback {
    capture noisily scpc, latlong backend(`backend')
    assert _rc == 999
}
replace s_1 = 91 in 2
foreach backend in rust auto fallback {
    capture noisily scpc, latlong backend(`backend')
    assert _rc == 198
}

* Regression weights have no SCPC profile on either backend.
generate double w = 1 + mod(_n, 2)
quietly regress y x [aw=w] if selected, vce(robust)
foreach backend in rust auto fallback {
    capture noisily scpc, uncond backend(`backend')
    assert _rc == 198
}

* The conditional projection is defined for 2SLS. Preserve unconditional LIML.
clear
set seed 958
set obs 30
generate double s_1 = _n
generate double z1 = rnormal()
generate double z2 = rnormal()
generate double u = rnormal()
generate double x = z1 + z2 + u
generate double y = x + u + rnormal()
quietly ivregress liml y (x = z1 z2), vce(robust)
foreach backend in rust auto fallback {
    capture noisily scpc, backend(`backend')
    assert _rc == 198
}
quietly scpc, uncond backend(fallback)
assert "`e(backend)'" == "fallback"

* Equivalent IV designs give the same results with redundant, factor, and lagged terms.
clear
set seed 203
set obs 30
generate double s_1 = _n
generate long time = _n
generate byte g = 1 + mod(_n, 3)
generate double d2 = g == 2
generate double d3 = g == 3
generate double z1 = d2 + .3 * rnormal()
generate double z2 = d3 + .3 * rnormal()
generate double x = z1 + rnormal()
generate double xcopy = x
generate double y = x + d2 + rnormal()

quietly ivregress 2sls y (x = z1 z2), vce(robust)
quietly scpc, avc(.5) k(2) cvs backend(fallback)
matrix ref_stats = e(scpcstats)
matrix ref_cvs = e(scpccvs)
quietly ivregress 2sls y (x xcopy = z1 z2), vce(robust)
foreach backend in rust fallback {
    quietly scpc, avc(.5) k(2) cvs backend(`backend')
    matrix stats = e(scpcstats)
    matrix cv = e(scpccvs)
    assert rowsof(stats) == 2
    mata: assert(max(abs(st_matrix("stats") - st_matrix("ref_stats"))) < 1e-8)
    mata: assert(max(abs(st_matrix("cv") - st_matrix("ref_cvs"))) < 1e-8)
}

quietly ivregress 2sls y (d2 d3 = z1 z2), vce(robust)
quietly scpc, avc(.5) k(3) cvs backend(fallback)
matrix ref_stats = e(scpcstats)
matrix ref_cvs = e(scpccvs)
quietly ivregress 2sls y (i.g = z1 z2), vce(robust)
quietly scpc, avc(.5) k(3) cvs backend(auto)
assert "`e(backend)'" == "fallback"
matrix stats = e(scpcstats)
matrix cv = e(scpccvs)
assert rowsof(stats) == 3
mata: assert(max(abs(st_matrix("stats") - st_matrix("ref_stats"))) < 1e-8)
mata: assert(max(abs(st_matrix("cv") - st_matrix("ref_cvs"))) < 1e-8)

quietly tsset time
generate double lag_x = L.x
generate double lag_z1 = L.z1
generate double lag_z2 = L.z2
quietly ivregress 2sls y (lag_x = lag_z1 lag_z2), vce(robust)
quietly scpc, avc(.5) k(2) cvs backend(fallback)
matrix ref_stats = e(scpcstats)
matrix ref_cvs = e(scpccvs)
quietly ivregress 2sls y (L.x = L.z1 L.z2), vce(robust)
quietly scpc, avc(.5) k(2) cvs backend(auto)
assert "`e(backend)'" == "fallback"
matrix stats = e(scpcstats)
matrix cv = e(scpccvs)
mata: assert(max(abs(st_matrix("stats") - st_matrix("ref_stats"))) < 1e-8)
mata: assert(max(abs(st_matrix("cv") - st_matrix("ref_cvs"))) < 1e-8)

* Cluster grouping must preserve lag values even when clusters cross time order.
preserve
generate byte cl = 1 + mod(time, 5)
replace s_1 = cl
quietly ivregress 2sls y (lag_x = lag_z1 lag_z2), vce(cluster cl)
quietly scpc, avc(.5) k(2) cvs backend(fallback)
matrix ref_stats = e(scpcstats)
matrix ref_cvs = e(scpccvs)
quietly ivregress 2sls y (L.x = L.z1 L.z2), vce(cluster cl)
quietly scpc, avc(.5) k(2) cvs backend(auto)
assert "`e(backend)'" == "fallback"
matrix stats = e(scpcstats)
matrix cv = e(scpccvs)
mata: assert(max(abs(st_matrix("stats") - st_matrix("ref_stats"))) < 1e-8)
mata: assert(max(abs(st_matrix("cv") - st_matrix("ref_cvs"))) < 1e-8)

* Clustered OLS with time-series lag operator L.x must also preserve time order.
quietly tsset time
quietly regress y lag_x, vce(cluster cl)
quietly scpc, avc(.5) k(2) cvs backend(fallback)
matrix ref_ols_stats = e(scpcstats)
matrix ref_ols_cvs = e(scpccvs)
quietly tsset time
quietly regress y L.x, vce(cluster cl)
quietly scpc, avc(.5) k(2) cvs backend(auto)
assert "`e(backend)'" == "fallback"
matrix ols_stats = e(scpcstats)
matrix ols_cv = e(scpccvs)
mata: assert(max(abs(st_matrix("ols_stats") - st_matrix("ref_ols_stats"))) < 1e-8)
mata: assert(max(abs(st_matrix("ols_cv") - st_matrix("ref_ols_cvs"))) < 1e-8)
restore

* Unconditional ivregress liml with omitted factor base level must compact omitted columns.
quietly ivregress liml y (d2 d3 = z1 z2), vce(robust)
quietly scpc, avc(.5) k(3) uncond cvs backend(fallback)
matrix ref_liml_stats = e(scpcstats)
quietly ivregress liml y (i.g = z1 z2), vce(robust)
quietly scpc, avc(.5) k(3) uncond cvs backend(auto)
matrix liml_stats = e(scpcstats)
assert rowsof(liml_stats) == 3
mata: assert(max(abs(st_matrix("liml_stats") - st_matrix("ref_liml_stats"))) < 1e-8)

* A zero intercept slot in the stored fit must leave the noconstant design intact.
program append_omitted_constant, eclass
    tempname b covariance bread
    matrix `b' = e(b), 0
    matrix `covariance' = (e(V), 0 \ 0, 0)
    matrix `bread' = (e(V_modelbased), 0 \ 0, 0)
    matrix colnames `b' = x _cons
    matrix colnames `covariance' = x _cons
    matrix rownames `covariance' = x _cons
    matrix colnames `bread' = x _cons
    matrix rownames `bread' = x _cons
    ereturn repost b=`b' V=`covariance', rename resize
    ereturn matrix V_modelbased = `bread'
end
quietly ivregress 2sls y (x = z1 z2), vce(robust) noconstant
quietly scpc, avc(.5) cvs backend(fallback)
matrix ref_stats = e(scpcstats)
matrix ref_cvs = e(scpccvs)
append_omitted_constant
foreach backend in rust fallback {
    quietly scpc, avc(.5) cvs backend(`backend')
    matrix stats = e(scpcstats)
    matrix cv = e(scpccvs)
    assert rowsof(stats) == 1
    mata: assert(max(abs(st_matrix("stats") - st_matrix("ref_stats"))) < 1e-8)
    mata: assert(max(abs(st_matrix("cv") - st_matrix("ref_cvs"))) < 1e-8)
}

* Load the Mata definitions directly to exercise conditioning beyond the
* range where Stata's 2SLS fit keeps both columns.
clear all
quietly do `"`scpc_source'"'
mata:
    Q = (.5, .5, .5 \ .5, -.5, .5 \ .5, .5, -.5 \ .5, -.5, -.5)
    Z = Q[., 1..2]
    X = Q[., 1], Q[., 1] + 1e-7 * Q[., 2] + Q[., 3]
    V = scpc_iv_projector(X, Z)
    assert(max(abs(cross(X, V) - I(2))) < 1e-8)
    assert(max(abs(cross(Q[., 3], V))) < 1e-8)
    X = Q[., 1], Q[., 3] + 1e-13 * Q[., 2]
end
capture mata: scpc_iv_projector(X, Z)
assert _rc == 999

* A projected design uniformly at roundoff scale is unidentified at the cutoff.
mata: Q = I(4); Z = Q[., 1..2]; X = Q[., 3] + 1e-14 * Q[., 1]
capture mata: scpc_iv_projector(X, Z)
assert _rc == 999
mata: X = (Q[., 3] + 1e-14 * Q[., 1]), (Q[., 4] + 1e-14 * Q[., 2])
capture mata: scpc_iv_projector(X, Z)
assert _rc == 999

* Redundant and zero instrument columns can make the instrument matrix wide.
mata:
    Q = I(4)
    Z = 1e-150 * Q[., 1], 1e150 * Q[., 2], Q[., 1], Q[., 2], J(4, 1, 0)
    X = Q[., 1], Q[., 2] + .5 * Q[., 3]
    V = scpc_iv_projector(X, Z)
    assert(max(abs(cross(X, V) - I(2))) < 1e-10)
end

* Selected-design rank uses the same relative SVD threshold in every host.
mata: x = (-1\1\-1\1); z = (-1\-1\1\1)
foreach unit in 1e-200 1 1e200 {
    mata: assert(cols(scpc_design_basis((x, `unit'*(x+1e-9*z)), 1)) == 2)
    foreach gap in 1e-11 1e-13 {
        capture mata: scpc_design_basis((x, `unit'*(x+`gap'*z)), 1)
        assert _rc == 999
    }
}

* Strict Gaussian rejection at, and immediately next to, a one-factor tie.
mata:
    r = 1 - 5e-15
    assert(getrp((1,r \ r,1), 1) == .5)
    assert(getrp((1,0,r \ 0,0,0 \ r,0,1), 1) == .5)
    loading = (5,3,4)'
    assert(getrp(loading * loading', 1) == 0)
    assert(getrp(loading * loading', 1-1e-15) == 1)
    for (i=-900; i<=900; i=i+900) {
        scale = scpc_power_two(i)
        a = 2401 * scale
        b = 49 * scale
        assert(getrp((a,b \ b,scale), 49) == 0)
        b = b * (1 - epsilon(1))
        assert(getrp((a,b \ b,scale), 49) == .5)
    }
    near = loading * loading'
    near[1,1] = 25 + 2^(-48)
    assert(abs(getrp(near,1) - .5) < 2e-7)
    assert(getrp((1,-1 \ -1,1), 1) == 0)
end
capture mata: getrp((0,1 \ 1,0), 1)
assert _rc == 3498

display "Conditional OLS cutoff tests passed."
