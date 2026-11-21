clear all
set more off
version 14
args repo_root backend
if `"`repo_root'"' == "" local repo_root "."
if "`backend'" == "" local backend "fallback rust"
adopath ++ "`repo_root'/scpc-stata/src"

import delimited using "`repo_root'/verification/data/bench_2000.csv", asdouble clear
drop s_1_clustered s_2_clustered
tempname full_stats full_cvs
foreach engine of local backend {
    local options approx avc(.03) seed(1) cvs nthreads(1) backend(`engine')
    quietly regress y x1 x2, vce(robust)
    quietly scpc, `options' k(2)
    assert "`e(backend)'" == "`engine'" & "`e(method)'" == "approx"
    matrix `full_stats' = e(scpcstats)
    matrix `full_cvs' = e(scpccvs)

    // x2 moves from the second coefficient to the only reported coefficient.
    quietly regress y x2 x1, vce(robust)
    quietly scpc, `options' k(1)
    assert "`e(backend)'" == "`engine'" & "`e(method)'" == "approx"
    mata: assert(mreldif(st_matrix("`full_stats'")[2,.], st_matrix("e(scpcstats)")) < 1e-10)
    mata: assert(mreldif(st_matrix("`full_cvs'")[2,.], st_matrix("e(scpccvs)")) < 1e-10)
    display "Nystrom coefficient order and subset passed (`engine')."
}
