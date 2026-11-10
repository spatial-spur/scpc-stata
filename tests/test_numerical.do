* Run from the repository root: do tests/test_numerical.do
clear all
set more off
do src/scpc.ado
mata:
latlongflag = 0
legacyprng_flag = 0
setGQxw(0)
s = (1::5), (1::5):^2
D = getdistmat(s)
assert(max(abs(D - getdistmat(s :+ 1e9))) < 1e-10)
assert(max(abs(getcrossdistmat(s, s) - D)) < 1e-10)
mixed = (0, 0 \ 1e9, 0 \ 1e9 + 1, 0)
assert(abs(getdistmat(mixed)[2,3] - 1) < 1e-10)
latlongflag = 1
near = (35, -70 \ 35 + 1e-7, -70)
assert(getdistmat(near)[1,2] > 0)
latlongflag = 0
for (q = 1; q <= 10; q++) {
    for (c = 0; c <= 3; c = c + 0.1) {
        assert(abs(getrp(I(q + 1), c) - 2 * ttail(q, sqrt(q) * c)) < 1e-8)
    }
}
for (i = 1; i <= 20; i++) {
    r = (1::100), (1::100)
    jumble_s(r, 20)
    assert(rows(uniqrows(r)) == 100)
}
for (scale = 1e-12; scale <= 1e12; scale = scale * 1e12) {
    c = getc0_from_mat(D * scale, 0.03)
    assert(abs((sum(exp(-c * D * scale)) - rows(D)) / (rows(D) * (rows(D) - 1)) - 0.03) < 1e-10)
}
capM = 2000
cgridfac = 1.2
W = J(5, 1, 2), (1::5)
Oms = LNgetOms(s, 1e6, 2e6, W)
assert(max(abs(Oms[1].mat - cross(W, W))) < 1e-10)
end
sysuse auto, clear
set seed 42
gen s_1 = runiform()
gen s_2 = runiform()
regress price weight length, robust
scpc, exact cvs
assert e(scpcstats)[1,4] >= 0 & e(scpcstats)[1,4] <= 1
