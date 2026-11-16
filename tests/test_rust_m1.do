clear all
set more off
version 14

local repo_root "/Users/pauldavidboll/Github/spatial-spur"
adopath ++ "`repo_root'/scpc-stata/src"

// Helper to compare Rust vs Fallback outputs within 1e-10
capture program drop check_rust_vs_fallback
program define check_rust_vs_fallback
    syntax , [cvs latlong uncond k(real 10) avc(real 0.03)]
    local opts "avc(`avc') k(`k')"
    if "`cvs'" != "" local opts "`opts' cvs"
    if "`latlong'" != "" local opts "`opts' latlong"
    if "`uncond'" != "" local opts "`opts' uncond"

    scpc, `opts' backend(rust)
    assert "`e(backend)'" == "rust"
    local r_q = e(q)
    local r_c0 = e(c0)
    local r_cv = e(cv)
    tempname r_stats r_cvs f_stats f_cvs
    matrix `r_stats' = e(scpcstats)
    if "`cvs'" != "" {
        matrix `r_cvs' = e(scpccvs)
    }

    scpc, `opts' backend(fallback)
    assert "`e(backend)'" == "fallback"
    local f_q = e(q)
    local f_c0 = e(c0)
    local f_cv = e(cv)
    matrix `f_stats' = e(scpcstats)
    if "`cvs'" != "" {
        matrix `f_cvs' = e(scpccvs)
    }

    assert `r_q' == `f_q'
    assert abs(`r_c0' - `f_c0') < 1e-10
    assert abs(`r_cv' - `f_cv') < 1e-10

    tempname d_stats d_cvs
    mata: st_numscalar("`d_stats'", max(abs(st_matrix("`r_stats'") - st_matrix("`f_stats'"))))
    di as text "  stats max diff = " %9.3e scalar(`d_stats') " (q=" `r_q' ", c0=" %12.9g `r_c0' ", cv=" %12.9g `r_cv' ")"
    assert scalar(`d_stats') < 1e-10

    if "`cvs'" != "" {
        mata: st_numscalar("`d_cvs'", max(abs(st_matrix("`r_cvs'") - st_matrix("`f_cvs'"))))
        di as text "  cvs   max diff = " %9.3e scalar(`d_cvs')
        assert scalar(`d_cvs') < 1e-10
    }
end

// =========================================================================
// Test 1: Direct ping + repeated scpc calls (plugin remains loaded)
// =========================================================================
di as text "=== Test 1: ping + repeated scpc calls ==="
program _test_ping_plug, plugin using("scpc_rust.plugin")
plugin call _test_ping_plug, ping
program drop _test_ping_plug

import delimited "`repo_root'/verification/data/bench_200.csv", clear
drop s_1_clustered s_2_clustered
regress y x1 x2, robust

scpc, cvs
assert "`e(backend)'" == "rust"
scpc
assert "`e(backend)'" == "rust"
scpc, uncond cvs
assert "`e(backend)'" == "rust"

// =========================================================================
// Test 2: bench_200.csv parity (Euclidean & latlong, cond & uncond, cvs)
// =========================================================================
di as text "=== Test 2: bench_200.csv Euclidean & latlong parity ==="
import delimited "`repo_root'/verification/data/bench_200.csv", clear
drop s_1_clustered s_2_clustered
regress y x1 x2, robust

di as text "--- bench_200 Euclidean conditional + cvs ---"
check_rust_vs_fallback, cvs
di as text "--- bench_200 Euclidean conditional (no cvs) ---"
check_rust_vs_fallback
di as text "--- bench_200 Euclidean unconditional + cvs ---"
check_rust_vs_fallback, uncond cvs
di as text "--- bench_200 Euclidean unconditional (no cvs) ---"
check_rust_vs_fallback, uncond

// Swap s_1=lat, s_2=lon (CSV has s_1=lon, s_2=lat)
rename s_1 tmp_lon
rename s_2 s_1
rename tmp_lon s_2
di as text "--- bench_200 latlong conditional + cvs ---"
check_rust_vs_fallback, latlong cvs
di as text "--- bench_200 latlong unconditional (no cvs) ---"
check_rust_vs_fallback, latlong uncond

// Double-precision load cross-check vs pinned M1 exact_scpc values
clear
set type double
import delimited "`repo_root'/verification/data/bench_200.csv", clear
set type float
drop s_1_clustered s_2_clustered
regress y x1 x2, robust
scpc, backend(rust) cvs
assert e(q) == 19
assert abs(e(c0) - 3.500409719601247) < 1e-10
assert abs(e(cv) - 2.1008365544083087) < 1e-10
matrix pinned = (1.5473280763431876, 0.09650902048915827, 16.032989128897157, 4.3785289611173374e-12, 1.3375098664953151, 1.75714628619106 \ ///
                 -0.6589528089605742, 0.10738682321940547, -6.13625386435211, 8.27465855135604e-06, -0.8878059409227403, -0.430099676998408 \ ///
                 1.4834167407175964, 0.11454779987636342, 12.95019845268713, 1.5882043004635851e-10, 1.236953655047813, 1.7298798263873798)
tempname d_pin
mata: st_numscalar("`d_pin'", max(abs(st_matrix("e(scpcstats)") - st_matrix("pinned"))))
di as text "  bench_200 double vs pinned M1 max diff = " %9.3e scalar(`d_pin')
assert scalar(`d_pin') < 1e-10

// =========================================================================
// Test 3: bench_1000.csv parity (Euclidean & latlong, cond & uncond, cvs)
// =========================================================================
di as text "=== Test 3: bench_1000.csv Euclidean & latlong parity ==="
import delimited "`repo_root'/verification/data/bench_1000.csv", clear
drop s_1_clustered s_2_clustered
regress y x1 x2, robust

di as text "--- bench_1000 Euclidean conditional + cvs ---"
check_rust_vs_fallback, cvs
di as text "--- bench_1000 Euclidean unconditional (no cvs) ---"
check_rust_vs_fallback, uncond

rename s_1 tmp_lon
rename s_2 s_1
rename tmp_lon s_2
di as text "--- bench_1000 latlong conditional + cvs ---"
check_rust_vs_fallback, latlong cvs
di as text "--- bench_1000 latlong unconditional + cvs ---"
check_rust_vs_fallback, latlong uncond cvs

// =========================================================================
// Test 4: nthreads(1) vs nthreads(4)
// =========================================================================
di as text "=== Test 4: nthreads(1) vs nthreads(4) ==="
import delimited "`repo_root'/verification/data/bench_200.csv", clear
drop s_1_clustered s_2_clustered
regress y x1 x2, robust

scpc, backend(rust) nthreads(1) cvs
assert "`e(backend)'" == "rust"
matrix t1_stats = e(scpcstats)
matrix t1_cvs = e(scpccvs)
local t1_q = e(q)
local t1_c0 = e(c0)
local t1_cv = e(cv)

scpc, backend(rust) nthreads(4) cvs
assert "`e(backend)'" == "rust"
assert e(q) == `t1_q'
assert abs(e(c0) - `t1_c0') < 1e-12
assert abs(e(cv) - `t1_cv') < 1e-12
tempname d_nt_stats d_nt_cvs
mata: st_numscalar("`d_nt_stats'", max(abs(st_matrix("e(scpcstats)") - st_matrix("t1_stats"))))
mata: st_numscalar("`d_nt_cvs'", max(abs(st_matrix("e(scpccvs)") - st_matrix("t1_cvs"))))
di as text "  nthreads(1) vs (4) diffs: stats=" %9.3e scalar(`d_nt_stats') " cvs=" %9.3e scalar(`d_nt_cvs')
assert scalar(`d_nt_stats') < 1e-12
assert scalar(`d_nt_cvs') < 1e-12

// =========================================================================
// Test 5: Edge cases (noconstant, intercept-only, k(1), if subset, d>=10)
// =========================================================================
di as text "=== Test 5: Edge cases ==="
// 5a: noconstant
regress y x1 x2, robust noconstant
check_rust_vs_fallback, cvs

// 5b: intercept-only
regress y, robust
check_rust_vs_fallback, cvs

// 5c: k(1) subset of coefficients
regress y x1 x2 x3, robust
check_rust_vs_fallback, k(1) cvs

// 5d: if subsetting + missing coord outside e(sample)
preserve
replace s_1 = . in 180/200
regress y x1 x2 if id <= 150, robust
check_rust_vs_fallback, cvs
restore

// 5e: d = 10 coordinates in non-alphabetical / non-numerical dataset order
preserve
gen s_10 = x3 * 0.1
forvalues j = 3/9 {
    gen s_`j' = x1 * 0.05 * `j'
}
regress y x1 x2, robust
check_rust_vs_fallback, cvs
restore

// =========================================================================
// Test 6: $SPATIAL_BACKEND global override & backend() precedence
// =========================================================================
di as text "=== Test 6: SPATIAL_BACKEND global & option precedence ==="
regress y x1 x2, robust
global SPATIAL_BACKEND "fallback"
scpc
assert "`e(backend)'" == "fallback"
scpc, backend(rust)
assert "`e(backend)'" == "rust"
global SPATIAL_BACKEND "rust"
scpc
assert "`e(backend)'" == "rust"
scpc, backend(fallback)
assert "`e(backend)'" == "fallback"
global SPATIAL_BACKEND ""

// =========================================================================
// Test 7: Fallback and backend(rust) rejection on unsupported M1 cases
// =========================================================================
di as text "=== Test 7: Unsupported Milestone 1 cases ==="
// 7a: Clustered OLS (use cluster-constant coordinates s_1_clustered, s_2_clustered)
import delimited "`repo_root'/verification/data/bench_200.csv", clear
drop s_1 s_2
rename s_1_clustered s_1
rename s_2_clustered s_2
regress y x1 x2, cluster(cluster_id)
scpc, uncond
assert "`e(backend)'" == "fallback"
capture scpc, uncond backend(rust)
assert _rc == 198

// 7b: ivregress 2sls
import delimited "`repo_root'/verification/data/bench_200.csv", clear
drop s_1_clustered s_2_clustered
ivregress 2sls y x2 (x1 = z1), robust
scpc
assert "`e(backend)'" == "fallback"
capture scpc, backend(rust)
assert _rc == 198

// 7c: Omitted/collinear e(b) drop (backend(rust) rejects; master Mata getcv() loops on zero-variance omitted var)
gen x_col = 2 * x1
regress y x1 x2 x_col, robust
capture scpc, backend(rust)
assert _rc == 198

// 7d: Active sample >= 4500
clear
set obs 4500
gen y = rnormal()
gen x1 = rnormal()
gen s_1 = rnormal()
gen s_2 = rnormal()
regress y x1, robust
capture scpc, backend(rust)
assert _rc == 198

// 7e: Invalid option arguments
capture scpc, backend(bogus)
assert _rc == 198
capture scpc, nthreads(0)
assert _rc == 198

// 7f: Missing plugin behavior in isolated ado directory
clear all
tempfile tmpado
local tmpdir "/tmp/scpc_no_plug_test"
shell mkdir -p "`tmpdir'"
copy "`repo_root'/scpc-stata/src/scpc.ado" "`tmpdir'/scpc.ado", replace
adopath - "`repo_root'/scpc-stata/src"
adopath ++ "`tmpdir'"
import delimited "`repo_root'/verification/data/bench_200.csv", clear
drop s_1_clustered s_2_clustered
regress y x1 x2, robust
scpc
assert "`e(backend)'" == "fallback"
capture scpc, backend(rust)
assert _rc == 198
shell rm -rf "`tmpdir'"

di as text "ALL MILESTONE 1 STATA TESTS PASSED!"
