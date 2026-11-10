* benchmark_scpc.do: Benchmark suite for scpc v2.1.0 across sample sizes
clear all
set more off

local curdir : pwd
adopath ++ "`curdir'/../src"

local n_list "200 722 1500 3000 5000 10000"
local k 3

di as text _n "=========================================================================================="
di as text "BENCHMARK SUITE: scpc v2.1.0 (Pure Mata)"
di as text "=========================================================================================="

foreach n of local n_list {
	drop _all
	set obs `n'
	set seed 42
	
	gen s_1 = runiform() * 10
	gen s_2 = runiform() * 10
	
	forval j = 1/`k' {
		gen x`j' = rnormal()
	}
	gen u = rnormal()
	gen y = 1.0 + 2.0 * x1 - 1.5 * x2 + u
	
	qui regress y x1-x`k', robust
	
	// Benchmark auto mode
	timer clear 1
	timer on 1
	qui scpc, cvs
	timer off 1
	qui timer list 1
	local t_auto = r(t1)
	
	local mode "Exact"
	if (`n' >= 10000) local mode "Approx (Large-N)"

	di as result "N = " %6.0f `n' " | Mode: " "`mode'" " | Runtime: " %7.3f `t_auto' " s"
}

di as text _n "Benchmark completed successfully."
