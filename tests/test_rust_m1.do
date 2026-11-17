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
    assert abs(`r_c0' - `f_c0') <= 1e-10 * max(1, abs(`f_c0'))
    assert abs(`r_cv' - `f_cv') < 1e-10

    tempname d_stats d_cvs
    mata: R = st_matrix("`r_stats'"); F = st_matrix("`f_stats'"); ///
          D = abs(F) :* (abs(F) :> 1) :+ (abs(F) :<= 1); ///
          st_numscalar("`d_stats'", max(abs(R - F) :/ D))
    di as text "  stats max rel/abs diff = " %9.3e scalar(`d_stats') " (q=" `r_q' ", c0=" %12.9g `r_c0' ", cv=" %12.9g `r_cv' ")"
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
plugin call _test_ping_plug, ping_v2
plugin call _test_ping_plug, ping_v3
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

// Double-precision load cross-check vs pinned M1 score statistics.
// The pinned conditional intervals use the historical weight order.
clear
set type double
import delimited "`repo_root'/verification/data/bench_200.csv", clear
set type float
drop s_1_clustered s_2_clustered
regress y x1 x2, robust
scpc, backend(rust) cvs
assert e(q) == 19
tempname avc_chk d_pin
mata: S = st_data(., ("s_1", "s_2")); v = J(0, 1, 0); for(i=1;i<rows(S);i++) v = v \ sqrt(rowsum((S[i+1::rows(S),.] :- S[i,.]):^2)); st_numscalar("`avc_chk'", mean(exp(-st_numscalar("e(c0)") * v)))
di as text "  bench_200 double achieved avc = " %18.15f scalar(`avc_chk')
assert abs(scalar(`avc_chk') - 0.03) < 1e-12
assert abs(e(c0) - 3.500409719601247) < 1e-3
assert abs(e(cv) - 2.1008365544083087) < 1e-3
matrix pinned = (1.5473280763431876, 0.09650902048915827, 16.032989128897157 \ ///
                 -0.6589528089605742, 0.10738682321940547, -6.13625386435211 \ ///
                 1.4834167407175964, 0.11454779987636342, 12.95019845268713)
mata: st_numscalar("`d_pin'", max(abs(st_matrix("e(scpcstats)")[.,1..3] - st_matrix("pinned"))))
di as text "  bench_200 double vs pinned M1 max diff = " %9.3e scalar(`d_pin')
assert scalar(`d_pin') < 1e-3
check_rust_vs_fallback, cvs

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
capture adopath - "`tmpdir'"
shell rm -rf "`tmpdir'"
adopath ++ "`repo_root'/scpc-stata/src"

// Re-define check_rust_vs_fallback after Test 7f's clear all
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
    assert abs(`r_c0' - `f_c0') <= 1e-10 * max(1, abs(`f_c0'))
    assert abs(`r_cv' - `f_cv') < 1e-10

    tempname d_stats d_cvs
    mata: R = st_matrix("`r_stats'"); F = st_matrix("`f_stats'"); ///
          D = abs(F) :* (abs(F) :> 1) :+ (abs(F) :<= 1); ///
          st_numscalar("`d_stats'", max(abs(R - F) :/ D))
    di as text "  stats max rel/abs diff = " %9.3e scalar(`d_stats') " (q=" `r_q' ", c0=" %12.9g `r_c0' ", cv=" %12.9g `r_cv' ")"
    assert scalar(`d_stats') < 1e-10

    if "`cvs'" != "" {
        mata: st_numscalar("`d_cvs'", max(abs(st_matrix("`r_cvs'") - st_matrix("`f_cvs'"))))
        di as text "  cvs   max diff = " %9.3e scalar(`d_cvs')
        assert scalar(`d_cvs') < 1e-10
    }
end

// =========================================================================
// Test 8: Scratch ownership, colliding user state, order preservation, & e(scpccvs)
// =========================================================================
di as text "=== Test 8: Scratch ownership & fresh e(scpccvs) ==="

capture program drop make_colliding_user_state
program define make_colliding_user_state
    // Scramble observation order so sortpreserve / Mata sorting must restore it
    tempvar perm
    gen double `perm' = mod(id * 73 + 19, _N)
    sort `perm' id
    drop `perm'
    assert id[1] != 1

    // Add colliding user variables
    gen double scpc_user_keep = 1000 + id * 1.5
    gen double scpc_score42 = -200.25 + id * 0.75
    gen double scpc_score1 = sin(id)
    gen double scpc_r = id * 2 + 3
    gen double scpc_rx = id * 3 - 5
    gen double scpc_x = cos(id)
    gen double scpc_xs = -id
    gen double scpc_fs_1 = id + 0.125

    // Add a user held estimate named oreg
    qui regress y x3
    capture _estimates drop oreg
    _estimates hold oreg

    // Add colliding user matrices (exact 2x2 non-trivial matrices)
    matrix b = (11, 12 \ 13, 14)
    matrix coef = (21, 22 \ 23, 24)
    matrix scpc_Bread = (31, 32 \ 33, 34)
    matrix scpc_Bread0 = (41, 42 \ 43, 44)
    matrix scpctab = (51, 52 \ 53, 54)
    matrix scpc_cvs = (61, 62 \ 63, 64)
    matrix scpcstats = (71, 72 \ 73, 74)
    matrix scpc_cvs_c = (81, 82 \ 83, 84)
end

capture program drop assert_user_state_preserved
program define assert_user_state_preserved
    syntax, exp_vars(string) exp_cmd(string) exp_depvar(string) exp_n(integer) exp_samp(varname)
    unab got_vars : _all
    assert `"`got_vars'"' == `"`exp_vars'"'

    // Dataset values and exact observation order match snapshot in Mata
    mata: assert(st_data(., tokens(st_local("exp_vars"))) == __t_data)

    // User held estimate oreg survived untouched
    qui _estimates dir
    assert strpos(" `r(names)' ", " oreg ") > 0

    // Colliding user matrices survived untouched
    mata: assert(st_matrix("b") == (11, 12 \ 13, 14))
    mata: assert(st_matrix("coef") == (21, 22 \ 23, 24))
    mata: assert(st_matrix("scpc_Bread") == (31, 32 \ 33, 34))
    mata: assert(st_matrix("scpc_Bread0") == (41, 42 \ 43, 44))
    mata: assert(st_matrix("scpctab") == (51, 52 \ 53, 54))
    mata: assert(st_matrix("scpc_cvs") == (61, 62 \ 63, 64))
    mata: assert(st_matrix("scpcstats") == (71, 72 \ 73, 74))
    mata: assert(st_matrix("scpc_cvs_c") == (81, 82 \ 83, 84))

    // Underlying fitted model in e() is preserved
    assert `"`e(cmd)'"' == `"`exp_cmd'"'
    assert `"`e(depvar)'"' == `"`exp_depvar'"'
    assert e(N) == `exp_n'
    mata: assert(st_matrix("e(b)") == __t_eb)
    mata: assert(st_matrix("e(V)") == __t_eV)
    mata: assert(st_matrix("e(V_modelbased)") == __t_eVm)
    tempvar cur_esamp
    gen byte `cur_esamp' = e(sample)
    assert `cur_esamp' == `exp_samp'
end

// 8a: Exact conditional OLS on subset: Rust + Fallback + backend/cvs toggling on same fit
import delimited "`repo_root'/verification/data/bench_200.csv", clear
drop s_1_clustered s_2_clustered
make_colliding_user_state

regress y x1 x2 if id <= 150, robust
assert e(df_r) == 147
tempvar samp_ols
gen byte `samp_ols' = e(sample)
unab ols_vars : _all
mata: __t_data = st_data(., tokens(st_local("ols_vars")))
mata: __t_eb = st_matrix("e(b)")
mata: __t_eV = st_matrix("e(V)")
mata: __t_eVm = st_matrix("e(V_modelbased)")

// Call 1: Rust with cvs
scpc, cvs backend(rust)
assert "`e(backend)'" == "rust"
confirm matrix e(scpccvs)
assert rowsof(e(scpccvs)) == 3 & colsof(e(scpccvs)) == 4
assert e(df_r) == 147
assert_user_state_preserved, exp_vars("`ols_vars'") exp_cmd("regress") exp_depvar("y") exp_n(150) exp_samp(`samp_ols')

// Call 2: Same fit, Fallback without cvs -> stale e(scpccvs) must be cleared
scpc, backend(fallback)
assert "`e(backend)'" == "fallback"
capture confirm matrix e(scpccvs)
assert _rc != 0
assert e(df_r) == 147
assert_user_state_preserved, exp_vars("`ols_vars'") exp_cmd("regress") exp_depvar("y") exp_n(150) exp_samp(`samp_ols')

// Call 3: Same fit, Fallback with cvs -> fresh e(scpccvs)
scpc, cvs backend(fallback)
assert "`e(backend)'" == "fallback"
confirm matrix e(scpccvs)
assert rowsof(e(scpccvs)) == 3 & colsof(e(scpccvs)) == 4
assert_user_state_preserved, exp_vars("`ols_vars'") exp_cmd("regress") exp_depvar("y") exp_n(150) exp_samp(`samp_ols')

// Call 4: Same fit, Rust without cvs -> clears fallback's e(scpccvs)
scpc, backend(rust)
assert "`e(backend)'" == "rust"
capture confirm matrix e(scpccvs)
assert _rc != 0
assert_user_state_preserved, exp_vars("`ols_vars'") exp_cmd("regress") exp_depvar("y") exp_n(150) exp_samp(`samp_ols')

// Call 5: Same fit, Rust with k(1) cvs -> 1x4 e(scpccvs), then Call 6 without cvs -> cleared
scpc, k(1) cvs backend(rust)
assert "`e(backend)'" == "rust"
confirm matrix e(scpccvs)
assert rowsof(e(scpccvs)) == 1 & colsof(e(scpccvs)) == 4
assert_user_state_preserved, exp_vars("`ols_vars'") exp_cmd("regress") exp_depvar("y") exp_n(150) exp_samp(`samp_ols')

scpc, backend(rust)
assert "`e(backend)'" == "rust"
capture confirm matrix e(scpccvs)
assert _rc != 0
assert_user_state_preserved, exp_vars("`ols_vars'") exp_cmd("regress") exp_depvar("y") exp_n(150) exp_samp(`samp_ols')

// Verify user held estimate oreg itself kept its original coefficients
_estimates unhold oreg
assert "`e(cmd)'" == "regress" & "`e(depvar)'" == "y"
mata: assert(cols(st_matrix("e(b)")) == 2)
_estimates drop _all

// 8b: Clustered fallback + logit conditional error after cluster sort & score creation
import delimited "`repo_root'/verification/data/bench_200.csv", clear
drop s_1 s_2
rename s_1_clustered s_1
rename s_2_clustered s_2
make_colliding_user_state

regress y x1 x2, cluster(cluster_id)
tempvar samp_cl
gen byte `samp_cl' = e(sample)
unab cl_vars : _all
mata: __t_data = st_data(., tokens(st_local("cl_vars")))
mata: __t_eb = st_matrix("e(b)")
mata: __t_eV = st_matrix("e(V)")
mata: __t_eVm = st_matrix("e(V_modelbased)")

scpc, cvs backend(fallback)
assert "`e(backend)'" == "fallback"
confirm matrix e(scpccvs)
assert rowsof(e(scpccvs)) == 3 & colsof(e(scpccvs)) == 4
assert_user_state_preserved, exp_vars("`cl_vars'") exp_cmd("regress") exp_depvar("y") exp_n(200) exp_samp(`samp_cl')

scpc, backend(fallback)
assert "`e(backend)'" == "fallback"
capture confirm matrix e(scpccvs)
assert _rc != 0
assert_user_state_preserved, exp_vars("`cl_vars'") exp_cmd("regress") exp_depvar("y") exp_n(200) exp_samp(`samp_cl')

// Error after cluster sort, scpc_setscores, and setOmsWfin on fallback path
gen byte y_bin = (y > 0)
logit y_bin x1 x2, vce(cluster cluster_id)
tempvar samp_log
gen byte `samp_log' = e(sample)
unab log_vars : _all
mata: __t_data = st_data(., tokens(st_local("log_vars")))
mata: __t_eb = st_matrix("e(b)")
mata: __t_eV = st_matrix("e(V)")
mata: __t_eVm = st_matrix("e(V_modelbased)")

capture scpc, backend(fallback)
assert _rc == 999
assert_user_state_preserved, exp_vars("`log_vars'") exp_cmd("logit") exp_depvar("y_bin") exp_n(200) exp_samp(`samp_log')

// 8c: IV 2SLS fallback + Rust runtime error after scratch score/matrix creation
import delimited "`repo_root'/verification/data/bench_200.csv", clear
drop s_1_clustered s_2_clustered
make_colliding_user_state

ivregress 2sls y x2 (x1 = z1), robust
tempvar samp_iv
gen byte `samp_iv' = e(sample)
unab iv_vars : _all
mata: __t_data = st_data(., tokens(st_local("iv_vars")))
mata: __t_eb = st_matrix("e(b)")
mata: __t_eV = st_matrix("e(V)")
mata: __t_eVm = st_matrix("e(V_modelbased)")

scpc, cvs backend(fallback)
assert "`e(backend)'" == "fallback"
confirm matrix e(scpccvs)
assert rowsof(e(scpccvs)) == 3 & colsof(e(scpccvs)) == 4
assert_user_state_preserved, exp_vars("`iv_vars'") exp_cmd("ivregress") exp_depvar("y") exp_n(200) exp_samp(`samp_iv')

scpc, backend(fallback)
assert "`e(backend)'" == "fallback"
capture confirm matrix e(scpccvs)
assert _rc != 0
assert_user_state_preserved, exp_vars("`iv_vars'") exp_cmd("ivregress") exp_depvar("y") exp_n(200) exp_samp(`samp_iv')

// Response rescaling (y_tiny = 1e-12 * y) succeeds & preserves q, t, p while scaling SE & CI
regress y x1 x2, robust
scpc, backend(rust) cvs
tempname orig_stats orig_cvs
matrix `orig_stats' = e(scpcstats)
matrix `orig_cvs' = e(scpccvs)
local orig_q = e(q)
local orig_c0 = e(c0)
local orig_cv = e(cv)

gen double y_tiny = 1e-12 * y
regress y_tiny x1 x2, robust
tempvar samp_tiny
gen byte `samp_tiny' = e(sample)
unab tiny_vars : _all
mata: __t_data = st_data(., tokens(st_local("tiny_vars")))
mata: __t_eb = st_matrix("e(b)")
mata: __t_eV = st_matrix("e(V)")
mata: __t_eVm = st_matrix("e(V_modelbased)")

scpc, backend(rust) cvs
assert "`e(backend)'" == "rust"
assert e(q) == `orig_q'
assert abs(e(c0) - `orig_c0') < 1e-12
assert abs(e(cv) - `orig_cv') < 1e-12
tempname tiny_stats tiny_cvs d_tiny_coef d_tiny_se d_tiny_t d_tiny_p d_tiny_ci d_tiny_cvs
matrix `tiny_stats' = e(scpcstats)
matrix `tiny_cvs' = e(scpccvs)
mata: O = st_matrix("`orig_stats'"); T = st_matrix("`tiny_stats'"); ///
      st_numscalar("`d_tiny_coef'", max(abs(T[.,1] - 1e-12 * O[.,1]))); ///
      st_numscalar("`d_tiny_se'", max(abs(T[.,2] - 1e-12 * O[.,2]))); ///
      st_numscalar("`d_tiny_t'", max(abs(T[.,3] - O[.,3]))); ///
      st_numscalar("`d_tiny_p'", max(abs(T[.,4] - O[.,4]))); ///
      st_numscalar("`d_tiny_ci'", max(abs(T[.,5..6] - 1e-12 * O[.,5..6]))); ///
      st_numscalar("`d_tiny_cvs'", max(abs(st_matrix("`tiny_cvs'") - st_matrix("`orig_cvs'"))))
assert scalar(`d_tiny_coef') < 1e-22
assert scalar(`d_tiny_se') < 1e-22
assert scalar(`d_tiny_t') < 1e-10
assert scalar(`d_tiny_p') < 1e-10
assert scalar(`d_tiny_ci') < 1e-22
assert scalar(`d_tiny_cvs') < 1e-10
assert_user_state_preserved, exp_vars("`tiny_vars'") exp_cmd("regress") exp_depvar("y_tiny") exp_n(200) exp_samp(`samp_tiny')

// Fallback also handles y_tiny and matches Rust
scpc, backend(fallback) cvs
assert "`e(backend)'" == "fallback"
assert e(q) == `orig_q'
tempname fb_tiny_stats d_rt_se d_rt_t d_rt_p
matrix `fb_tiny_stats' = e(scpcstats)
mata: R = st_matrix("`tiny_stats'"); F = st_matrix("`fb_tiny_stats'"); ///
      st_numscalar("`d_rt_se'", max(abs(R[.,2] - F[.,2]))); ///
      st_numscalar("`d_rt_t'", max(abs(R[.,3] - F[.,3]))); ///
      st_numscalar("`d_rt_p'", max(abs(R[.,4] - F[.,4])))
assert scalar(`d_rt_se') < 1e-22
assert scalar(`d_rt_t') < 1e-10
assert scalar(`d_rt_p') < 1e-10
assert_user_state_preserved, exp_vars("`tiny_vars'") exp_cmd("regress") exp_depvar("y_tiny") exp_n(200) exp_samp(`samp_tiny')

// Error inside plugin AFTER scpc_setscores and temporary result matrices exist (identical coordinates)
replace s_1 = 1
replace s_2 = 1
regress y x1 x2, robust
tempvar samp_id
gen byte `samp_id' = e(sample)
unab id_vars : _all
mata: __t_data = st_data(., tokens(st_local("id_vars")))
mata: __t_eb = st_matrix("e(b)")
mata: __t_eV = st_matrix("e(V)")
mata: __t_eVm = st_matrix("e(V_modelbased)")

capture scpc, backend(rust)
assert _rc == 198
assert_user_state_preserved, exp_vars("`id_vars'") exp_cmd("regress") exp_depvar("y") exp_n(200) exp_samp(`samp_id')

capture scpc, backend(fallback)
assert _rc == 999
assert_user_state_preserved, exp_vars("`id_vars'") exp_cmd("regress") exp_depvar("y") exp_n(200) exp_samp(`samp_id')

_estimates drop _all
mata: mata drop __t_data __t_eb __t_eV __t_eVm

// =========================================================================
// Test 9: Coordinate exactness, tied eigenspace fallback, & rescaling parity
// =========================================================================
di as text "=== Test 9: Coordinate exactness, tied fallback, & rescaling ==="

// 9a: Lone s_10 (and gapped s_1 s_10) must not abbreviate to s_1
import delimited "`repo_root'/verification/data/bench_200.csv", clear
drop s_1 s_2 s_1_clustered s_2_clustered
gen double s_10 = x1
regress y x1 x2, robust

capture scpc
assert _rc == 999
capture scpc, backend(fallback)
assert _rc == 999
capture scpc, backend(rust)
assert _rc == 198

gen double s_1 = x2
capture scpc
assert _rc == 999
capture scpc, backend(fallback)
assert _rc == 999
capture scpc, backend(rust)
assert _rc == 198

// 9b: Tied eigenspace automatic fallback on first 20 rows of bench_200.csv
import delimited "`repo_root'/verification/data/bench_200.csv", clear
drop s_1_clustered s_2_clustered
preserve
keep in 1/20
regress y x1 x2, robust

scpc, cvs
assert "`e(backend)'" == "fallback"
assert e(q) == 7
assert e(tied) == 0
tempname auto_stats auto_cvs
matrix `auto_stats' = e(scpcstats)
matrix `auto_cvs' = e(scpccvs)
local auto_c0 = e(c0)
local auto_cv = e(cv)

scpc, cvs backend(fallback)
assert "`e(backend)'" == "fallback"
assert e(q) == 7
assert e(tied) == 0
assert abs(e(c0) - `auto_c0') < 1e-12
assert abs(e(cv) - `auto_cv') < 1e-12
tempname d_tie_stats d_tie_cvs
mata: st_numscalar("`d_tie_stats'", max(abs(st_matrix("e(scpcstats)") - st_matrix("`auto_stats'"))))
mata: st_numscalar("`d_tie_cvs'", max(abs(st_matrix("e(scpccvs)") - st_matrix("`auto_cvs'"))))
assert scalar(`d_tie_stats') < 1e-12
assert scalar(`d_tie_cvs') < 1e-12

scpc, cvs backend(rust)
assert "`e(backend)'" == "rust"
assert e(q) == 7
assert e(tied) == 1
restore

// Also check tied eigenspace fallback with `if _n <= 20` subsetting
preserve
regress y x1 x2 if _n <= 20, robust
scpc
assert "`e(backend)'" == "fallback"
assert e(q) == 7
assert e(tied) == 0
restore

// Subsequent untied call clears e(tied) back to 0
regress y x1 x2, robust
scpc
assert "`e(backend)'" == "rust"
assert e(tied) == 0

// 9c: Coordinate rescaling (s_1 * 1e-6, s_2 * 1e-6) parity & invariance
preserve
recast double s_1 s_2
scpc, backend(rust) cvs
tempname base_stats base_cvs
matrix `base_stats' = e(scpcstats)
matrix `base_cvs' = e(scpccvs)
local base_q = e(q)
local base_c0 = e(c0)
local base_cv = e(cv)

replace s_1 = s_1 * 1e-6
replace s_2 = s_2 * 1e-6
regress y x1 x2, robust
check_rust_vs_fallback, cvs
scpc, backend(rust) cvs
assert e(q) == `base_q'
assert abs(e(c0) * 1e-6 - `base_c0') / `base_c0' < 1e-10
assert abs(e(cv) - `base_cv') < 1e-10
tempname d_cs_stats d_cs_cvs
mata: st_numscalar("`d_cs_stats'", max(abs(st_matrix("e(scpcstats)") - st_matrix("`base_stats'"))))
mata: st_numscalar("`d_cs_cvs'", max(abs(st_matrix("e(scpccvs)") - st_matrix("`base_cvs'"))))
assert scalar(`d_cs_stats') < 1e-10
assert scalar(`d_cs_cvs') < 1e-10
restore

// 9d: Small-unit regressor (x1 * 1e-10, smallest decade retained by Stata regress) parity & scale invariance
regress y x1 x2, robust
scpc, backend(rust) cvs
matrix `base_stats' = e(scpcstats)
matrix `base_cvs' = e(scpccvs)
local base_q = e(q)

preserve
gen double x1_tiny = x1 * 1e-10
regress y x1_tiny x2, robust
check_rust_vs_fallback, cvs
scpc, backend(rust) cvs
assert e(q) == `base_q'
tempname xs_stats xs_cvs d_xs_c1 d_xs_se1 d_xs_tp d_xs_rest d_xs_cvs
matrix `xs_stats' = e(scpcstats)
matrix `xs_cvs' = e(scpccvs)
mata: B = st_matrix("`base_stats'"); X = st_matrix("`xs_stats'"); ///
      st_numscalar("`d_xs_c1'", abs(X[1,1] * 1e-10 - B[1,1])); ///
      st_numscalar("`d_xs_se1'", abs(X[1,2] * 1e-10 - B[1,2])); ///
      st_numscalar("`d_xs_tp'", max(abs(X[1,3..4] - B[1,3..4]))); ///
      st_numscalar("`d_xs_rest'", max(abs(X[2..3,.] - B[2..3,.]))); ///
      st_numscalar("`d_xs_cvs'", max(abs(st_matrix("`xs_cvs'") - st_matrix("`base_cvs'"))))
assert scalar(`d_xs_c1') < 1e-10
assert scalar(`d_xs_se1') < 1e-10
assert scalar(`d_xs_tp') < 1e-10
assert scalar(`d_xs_rest') < 1e-10
assert scalar(`d_xs_cvs') < 1e-10
restore

// 9e: Co-location floor (attainable & unattainable avc with duplicate coordinates)
preserve
keep if inlist(id, 1, 2, 4, 5, 8, 9, 10, 11)
replace s_1 = s_1[1] in 2
replace s_2 = s_2[1] in 2
regress y x1 x2, robust

// 1 co-located pair out of 28 -> floor = 1/28 = 0.0357; avc(0.03) is unattainable
capture scpc, avc(0.03) backend(rust)
assert _rc == 198
capture scpc, avc(0.03) backend(fallback)
assert _rc == 999
capture scpc, avc(0.03)
assert _rc == 198

// avc(0.20) exceeds floor; minavc_eff handles cmax above floor cleanly and matches
check_rust_vs_fallback, avc(0.20) cvs
scpc, avc(0.20)
assert "`e(backend)'" == "rust"
assert e(tied) == 0
restore

di as text "ALL MILESTONE 1 STATA TESTS PASSED!"
