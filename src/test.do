* test.do: Comprehensive test suite for scpc v2.1.0 in scpc-stata
clear all
set more off

// Dynamically add the directory containing this do-file to adopath
local curdir : pwd
adopath ++ "`curdir'"

di as text _n "=========================================================="
di as text "TEST 1: Basic OLS regression with scpc default (Pure Mata mode)"
di as text "=========================================================="
sysuse auto, clear
set seed 12345
gen s_1 = runiform() * 10
gen s_2 = runiform() * 10

regress price weight length foreign
scpc, cvs
matrix list e(scpcstats)
matrix list e(scpccvs)

di as text _n "=========================================================="
di as text "TEST 2: Explicit exact mode (exact option)"
di as text "=========================================================="
scpc, cvs exact
matrix list e(scpcstats)

di as text _n "=========================================================="
di as text "TEST 3: Spherical Lat/Long coordinates"
di as text "=========================================================="
replace s_1 = 30 + runiform() * 15 // Latitude
replace s_2 = -100 + runiform() * 20 // Longitude
scpc, latlong
matrix list e(scpcstats)

di as text _n "=========================================================="
di as text "TEST 4: Large-N Nyström approximation (approx option)"
di as text "=========================================================="
scpc, approx
matrix list e(scpcstats)

di as text _n "=========================================================="
di as text "TEST 5: Seed reproducibility under Large-N approx"
di as text "=========================================================="
set seed 98765
scpc, approx
matrix S1 = e(scpcstats)

set seed 98765
scpc, approx
matrix S2 = e(scpcstats)

mata: st_numscalar("diff", max(abs(st_matrix("S1") - st_matrix("S2"))))
di as text "Max absolute difference between identical seeds: " %9.2e scalar(diff)
assert diff < 1e-12

di as text _n "=========================================================="
di as text "TEST 6: Legacy PRNG mode for backward compatibility"
di as text "=========================================================="
scpc, approx legacyprng
matrix list e(scpcstats)

di as text _n "=========================================================="
di as text "TEST 7: Custom cgridfac and threshold options"
di as text "=========================================================="
scpc, cgridfac(1.1) threshold(50)
matrix list e(scpcstats)

di as text _n "=========================================================="
di as text "TEST 8: Clustered standard errors"
di as text "=========================================================="
gen state = rep78
replace state = 1 if missing(state)
// Ensure coordinates are constant within cluster
bysort state: replace s_1 = s_1[1]
bysort state: replace s_2 = s_2[1]

regress price weight length, vce(cluster state)
scpc
matrix list e(scpcstats)

di as text _n "=========================================================="
di as text "TEST 9: IV regression (ivregress 2sls)"
di as text "=========================================================="
drop state
replace s_1 = runiform() * 10
replace s_2 = runiform() * 10
ivregress 2sls price length (weight = trunk headroom gear_ratio)
scpc
matrix list e(scpcstats)

di as text _n "=========================================================="
di as text "ALL TESTS COMPLETED SUCCESSFULLY!"
di as text "=========================================================="
