* Run from src so Stata loads scpc.ado from this checkout.
clear all
set more off
which scpc

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

expand 2
bysort s_1: replace x = x + cond(_n == 1, -.25, .25)
replace y = x + (s_1 == 0)

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

display "Conditional OLS cutoff tests passed."
