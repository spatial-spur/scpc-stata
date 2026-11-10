*! version 2.1.0 10nov2026
*! High-Performance Pure Mata Spatial Correlation Principal Components (SCPC)
*! Original Authors: Ulrich K. Müller & Mark W. Watson
*! Performance & Numerical Optimization: Paul David Boll & Antigravity

capture program drop scpc
program define scpc, eclass sortpreserve
	version 13.0
	syntax , ///
		[ ///
			cvs ///
			k(integer 10) ///
			avc(real -1) ///
			latlong ///
			uncond ///
			exact ///
			approx ///
			threshold(integer 10000) ///
			cgridfac(real 1.2) ///
			legacyprng ///
		]

	capture matrix scpc_Bread = e(V_modelbased)
	if (_rc != 0 | scpc_Bread[1, 1] == .) {
		capture matrix scpc_Bread = e(V)
		if (_rc != 0 | scpc_Bread[1, 1] == .) {
			disp as error "neither e(V_modelbased) nor e(V) found; aborting"
			exit(999)
		}
	}

	if (`avc' == -1) {
		local avc 0.03	// default average correlation
	}
	if (`avc' < 0.0001 | `avc' > 0.99) {
		disp as error "average correlation must be between 0.0001 and 0.99; aborting"
		exit(999)
	}

	if (`cgridfac' <= 1.0) {
		disp as error "cgridfac must be strictly greater than 1.0; aborting"
		exit(999)
	}

	// ENSURE COORDINATES ARE UNIQUE WITHIN CLUSTERS
	if ("`e(clustvar)'" != "") {
		tempvar es
		gen `es' = 0 if e(sample)
		sort `e(clustvar)' `es'
		unab tmp_svarlist: s_*
		qui foreach tmp_svar of varlist `tmp_svarlist' {
			tempvar unique_coords_tag unique_coords_num
			egen `unique_coords_tag' = tag(`e(clustvar)' `tmp_svar')
			egen `unique_coords_num' = total(`unique_coords_tag'), by(`e(clustvar)')
			sum `unique_coords_num'
			if (`r(min)' != 1 | `r(max)' != 1) {
				disp as error "location coordinates in `tmp_svar' vary within clusters defined by `e(clustvar)'; aborting"
				exit(999)
			}
		}
	}

	tempvar scpc_sel
	qui scpc_setscores, scpc_sel(`scpc_sel')

	// Determine method flag
	local method "auto"
	if ("`exact'" != "" & "`approx'" != "") {
		disp as error "options exact and approx cannot be combined; aborting"
		exit(198)
	}
	if ("`exact'" != "") local method "exact"
	if ("`approx'" != "") local method "approx"

	local use_legacy = ("`legacyprng'" != "")

	// Call pure Mata spatial basis setup
	mata set_spatial_scpc(`avc', "`scpc_sel'", "`latlong'", "`method'", `threshold', `cgridfac', `use_legacy')

	qui sum `scpc_sel', detail
	tempname neff
	matrix `neff' = r(sum)
	matrix b = e(b)
	capture matrix scpc_Bread0 = e(V_modelbased)
	if (_rc != 0 | scpc_Bread0[1, 1] == .) {
		matrix scpc_Bread0 = e(V)
	}

	local clist
	local na : colnames b
	if (e(cmd) == "regress") {
		tokenize `na'
		forval j = 1/`= colsof(b)' {
			if ("``j''" != "_cons") {
				local clist `clist' ``j''
			}
		}
		local clist = "qui reg scpc_rx `clist'"
	}
	else if (e(cmd) == "ivregress") {
		local elist `e(insts)'
		local i = `: word count `e(instd)'' + 1

		tempvar ess
		qui gen `ess' = 0
		qui replace `ess' = 1 if e(sample)
		_estimates hold oreg, restore
		cap drop scpc_fs_*
		local clist
		tokenize `na'
		forval j = 1/`= colsof(b)' {
			if (`j' < `i') {
				qui reg ``j'' `elist' if `ess'
				qui predict scpc_fs_`j', xb
				qui replace scpc_fs_`j' = . if !`ess'
				local clist `clist' scpc_fs_`j'
			}
			else {
				local clist `clist' ``j''
			}
		}
		matrix colnames scpc_Bread0 = `clist'
		_estimates unhold oreg
		local clist = "qui ivregress 2sls scpc_rx `e(exogr)' (`e(instd)' = `e(insts)' )"
	}
	else if ("`uncond'" == "") {
		disp as error "conditional critical values only supported for ivregress 2sls and regress, use option uncond for only unconditionally valid critical values; aborting"
		exit(999)
	}

	local k = min(colsof(b), `k')
	matrix b = b[1..., 1..`k']
	local na : colnames b
	matrix scpctab = b' * (1, 1, 1, 1, 1, 1)
	matrix scpc_cvs = b' * (1, 1, 1, 1)

	tokenize `na'

	// Loop through covariates
	local i = 1
	while ("``i''" != "") {
		tempvar w
		matrix coef = `neff' * scpc_Bread[`i', 1...]
		qui matrix score `w' = coef if `scpc_sel' == 1
		qui replace `w' = `w' + b[1, `i'] if `scpc_sel' == 1

		if ("`uncond'" == "") {
			capture drop scpc_x
			capture drop scpc_xs
			matrix coef = `neff' * scpc_Bread0[`i', 1...]
			qui matrix score scpc_x = coef if e(sample)

			if ("`e(clustvar)'" == "") {
				qui gen scpc_xs = sign(scpc_x) if e(sample)
				mata set_Wx_nocluster("`clist'", "`scpc_sel'")
			}
			else {
				qui by `e(clustvar)': egen scpc_xs = total(scpc_x * scpc_x) if e(sample)
				qui replace scpc_xs = scpc_x / sqrt(scpc_xs) if e(sample) & scpc_xs > 0
				mata set_Wx_cluster("`clist'", "`scpc_sel'")
			}
		}
		mata set_scpcstats("`w'", "`scpc_sel'")
		matrix scpctab[`i', 1] = scpcstats[1, 1..6]

		if ("`cvs'" == "cvs") {
			mata set_scpccvs()
			matrix scpc_cvs[`i', 1] = scpc_cvs_c[1, 1..4]
		}
		local ++i
	}

	matrix colnames scpctab = "Coef" "Std_Err" "  t  " "P>|t|" "95% Conf" "Interval"
	local n = rowsof(scpctab)
	local ands = `n' * "&"
	local rs &-`ands'
	local pavc : di %5.3f `avc'
	local tit "SCPC Inference for first `k' coefficients"
	matlist scpctab, border(all) title(`tit') cspec(o2& %12s | %9.0g o2 & %9.0g o2 &o1 %5.2f o1& o2 %6.3f o1 & o2 %9.0g & o1 %9.0g o2&) rspec(`rs')

	if ("`cvs'" == "cvs") {
		matrix colnames scpc_cvs = "32%" "10%" "5%" "1%"
		local n = rowsof(scpc_cvs)
		local ands = `n' * "&"
		local rs &-`ands'
		matlist scpc_cvs, border(all) title("Two-Sided Critical values of SCPC t-tests") cspec(o2& %12s | %6.3f o2 & %6.3f o2 & %6.3f o2 & %6.3f o2 &) rspec(`rs')
		ereturn matrix scpccvs = scpc_cvs
	}

	ereturn matrix scpcstats = scpctab
	cap drop scpc_*
end

capture program drop scpc_setscores
program define scpc_setscores, sortpreserve
	syntax, scpc_sel(name)
	cap drop scpc_score*
	tempvar e
	local slist ""
	quietly {
		if (e(cmd) == "regress" | e(cmd) == "logit" | e(cmd) == "probit") {
			predict `e' if e(sample), sc
			local na : colnames e(V)
			tokenize `na'
			local i = 1
			while "``i''" != "" {
				capture drop scpc_score`i'
				if ("``i''" == "_cons") gen scpc_score`i' = `e' if e(sample)
				else gen scpc_score`i' = `e' * ``i'' if e(sample)
				local slist "`slist' scpc_score`i'"
				local ++i
			}
		}
		else if (e(cmd) == "areg") {
			predict `e' if e(sample), sc
			sort `e(absvar)'
			local na : colnames e(V)
			tokenize `na'
			local i = 1
			while "``i''" != "" {
				capture drop scpc_score`i'
				if ("``i''" == "_cons") gen scpc_score`i' = `e' if e(sample)
				else {
					tempvar x
					by `e(absvar)': egen `x' = mean(``i'') if e(sample)
					sum ``i'' if e(sample)
					gen scpc_score`i' = `e' * (``i'' - `x' + `r(mean)') if e(sample)
				}
				local slist "`slist' scpc_score`i'"
				local ++i
			}
		}
		else if (e(cmd) == "ivregress") {
			tempvar stata_scr
			predict `stata_scr'*, sc
			local j = `: word count `e(instd)'' + 1
			gen `e' = `stata_scr'`j'
			local na : colnames e(V)
			tokenize `na'
			local i = 1
			while "``i''" != "" {
				capture drop scpc_score`i'
				if (`i' < `j') gen scpc_score`i' = `stata_scr'`i' if e(sample)
				else {
					if ("``i''" == "_cons") gen scpc_score`i' = `e' if e(sample)
					else gen scpc_score`i' = `e' * ``i'' if e(sample)
				}
				local slist "`slist' scpc_score`i'"
				local ++i
			}
		}
		else {
			disp as error "scpc does not support `e(cmd)'"
			exit(999)
		}

		matrix colnames scpc_Bread = `slist'
		gen `scpc_sel' = 0
		if ("`e(clustvar)'" == "") qui replace `scpc_sel' = 1 if e(sample)
		else {
			tempvar TotScore
			local i = 1
			while "``i''" != "" {
				capture drop `TotScore'
				by `e(clustvar)': egen `TotScore' = total(scpc_score`i') if e(sample)
				replace scpc_score`i' = `TotScore'
				by `e(clustvar)': replace `scpc_sel' = 1 if _n == 1 & e(sample)
				local ++i
			}
		}
	}
end

mata:

// Constant Pi with full 64-bit IEEE 754 precision
real scalar PI_CONST() {
	return(3.1415926535897932384626433832795)
}

struct mats {
	real matrix mat
}

struct sms {
	real matrix mat
	real matrix s
}

real matrix getdistvec(real matrix s1, real matrix s2)
{
	external real scalar latlongflag
	real vector dist
	real scalar c_rad
	real matrix d
	real colvector cos_s1, cos_s2, inside

	if (latlongflag == 0) {
		dist = sqrt(rowsum((s1 :- s2):^2))
	}
	else {
		c_rad = PI_CONST() / 180.0
		d = (0.5 * c_rad) * (s1 :- s2)
		cos_s1 = cos(c_rad * s1[., 1])
		cos_s2 = cos(c_rad * s2[., 1])
		inside = sin(d[., 1]):^2 :+ cos_s1 :* (cos_s2 :* (sin(d[., 2]):^2))
		inside = inside :* (inside :> 0)
		inside = inside :* (inside :<= 1) :+ (inside :> 1)
		dist = asin(sqrt(inside)) / PI_CONST()
	}
	return(dist)
}

// 100% Vectorized BLAS Level-3 Euclidean and Spherical Geodesic Distance Matrix
real matrix getdistmat(real matrix s)
{
	external real scalar latlongflag
	real matrix mat, sq_mat, P, dot, inside_chord
	real colvector sq, s_lat, s_lon
	real scalar n, c_rad

	n = rows(s)
	if (latlongflag == 0) {
		// BLAS Level-3 Euclidean Distance Gram Matrix: D^2 = sq*1' + 1*sq' - 2*(s*s')
		sq = rowsum(s:^2)
		sq_mat = sq * J(1, n, 1)
		mat = (sq_mat + sq_mat') - 2 * cross(s', s')
		_diag(mat, 0)
		mat = sqrt(mat :* (mat :> 0))
	}
	else {
		// 3D Chord Spherical Distance on Unit Sphere via BLAS Level-3 Dot Product
		c_rad = PI_CONST() / 180.0
		s_lat = s[., 1] :* c_rad
		s_lon = s[., 2] :* c_rad
		P = cos(s_lat):*cos(s_lon), cos(s_lat):*sin(s_lon), sin(s_lat)
		dot = cross(P', P')
		inside_chord = (1 :- dot) :/ 2
		_diag(inside_chord, 0)
		inside_chord = inside_chord :* (inside_chord :> 0)
		inside_chord = inside_chord :* (inside_chord :<= 1) :+ (inside_chord :> 1)
		mat = asin(sqrt(inside_chord)) / PI_CONST()
	}
	return(mat)
}

// Vectorized BLAS Level-3 Cross-Distance Matrix between two coordinate sets (n1 x n2)
real matrix getcrossdistmat(real matrix s1, real matrix s2)
{
	external real scalar latlongflag
	real matrix D2, sq1, sq2, P1, P2, dot, inside
	real scalar n1, n2, c_rad
	real colvector s1_lat, s1_lon, s2_lat, s2_lon

	n1 = rows(s1)
	n2 = rows(s2)
	if (latlongflag == 0) {
		sq1 = rowsum(s1:^2)
		sq2 = rowsum(s2:^2)
		D2 = (sq1 * J(1, n2, 1) + J(n1, 1, 1) * sq2') - 2 * (s1 * s2')
		return(sqrt(D2 :* (D2 :> 0)))
	}
	else {
		c_rad = PI_CONST() / 180.0
		s1_lat = s1[., 1] :* c_rad
		s1_lon = s1[., 2] :* c_rad
		s2_lat = s2[., 1] :* c_rad
		s2_lon = s2[., 2] :* c_rad
		P1 = cos(s1_lat):*cos(s1_lon), cos(s1_lat):*sin(s1_lon), sin(s1_lat)
		P2 = cos(s2_lat):*cos(s2_lon), cos(s2_lat):*sin(s2_lon), sin(s2_lat)
		dot = cross(P1', P2')
		inside = (1 :- dot) :/ 2
		inside = inside :* (inside :> 0)
		inside = inside :* (inside :<= 1) :+ (inside :> 1)
		return(asin(sqrt(inside)) / PI_CONST())
	}
}

real vector lvech(real matrix mat)
{
	real scalar i, c
	real vector v
	v = J(rows(mat) * (rows(mat) - 1) / 2, 1, 0)
	c = 1
	for (i = 1; i <= rows(mat) - 1; i++) {
		v[c::c + rows(mat) - i - 1] = mat[i + 1::rows(mat), i]
		c = c + rows(mat) - i
	}
	return(v)
}

real matrix demeanmat(real matrix mat)
{
	mat = mat :- (rowsum(mat) / cols(mat))
	mat = mat :- (colsum(mat) / rows(mat))
	return(mat)
}

// Direct Matrix-Sum Newton-Raphson Spatial Scale Solver (Eliminates lvech)
real scalar getc0_from_mat(real matrix D, real scalar avc0)
{
	real scalar n, denom, dbar, c, f, fprime, step, iter, c_low, c_high
	real matrix exp_cd

	n = rows(D)
	denom = n * (n - 1)
	dbar = sum(D) / denom
	c = -log(avc0) / dbar
	if (c <= 0) c = 1.0

	c_low = 0
	c_high = 1e8

	for (iter = 1; iter <= 25; iter++) {
		exp_cd = exp(-c * D)
		f = (sum(exp_cd) - n) / denom - avc0
		if (abs(f) < 1e-8 * avc0) break

		if (f > 0) c_low = max((c_low, c))
		else c_high = min((c_high, c))

		fprime = -sum(D :* exp_cd) / denom
		step = -f / fprime

		if (c + step > c_low & c + step < c_high & abs(step) < 2 * c) {
			c = c + step
		}
		else {
			c = (c_high < 1e7 ? 0.5 * (c_low + c_high) : 2 * c)
		}
	}
	return(c)
}

// Safeguarded Vector Newton-Raphson Spatial Scale Solver for Large-N Nyström mode
real scalar getc0_from_vec(real vector d, real scalar avc0)
{
	real scalar dbar, c, f, fprime, step, iter, c_low, c_high, n_d
	real vector exp_cd, d_sub

	n_d = length(d)
	if (n_d > 200000) {
		d_sub = d[1::100000]
		dbar = mean(d_sub)
		c = -log(avc0) / dbar
		if (c <= 0) c = 1.0
		c_low = 0
		c_high = 1e8
		for (iter = 1; iter <= 25; iter++) {
			exp_cd = exp(-c * d_sub)
			f = mean(exp_cd) - avc0
			if (abs(f) < 1e-5 * avc0) break
			if (f > 0) c_low = max((c_low, c))
			else c_high = min((c_high, c))
			fprime = -mean(d_sub :* exp_cd)
			step = -f / fprime
			if (c + step > c_low & c + step < c_high & abs(step) < 2 * c) {
				c = c + step
			}
			else {
				c = (c_high < 1e7 ? 0.5 * (c_low + c_high) : 2 * c)
			}
		}
		// Full-vector polishing steps
		for (iter = 1; iter <= 3; iter++) {
			exp_cd = exp(-c * d)
			f = mean(exp_cd) - avc0
			if (abs(f) < 1e-8 * avc0) break
			fprime = -mean(d :* exp_cd)
			c = c - f / fprime
		}
		return(c)
	}

	dbar = mean(d)
	c = -log(avc0) / dbar
	if (c <= 0) c = 1.0

	c_low = 0
	c_high = 1e8

	for (iter = 1; iter <= 25; iter++) {
		exp_cd = exp(-c * d)
		f = mean(exp_cd) - avc0
		if (abs(f) < 1e-8 * avc0) break

		if (f > 0) c_low = max((c_low, c))
		else c_high = min((c_high, c))

		fprime = -mean(d :* exp_cd)
		step = -f / fprime

		if (c + step > c_low & c + step < c_high & abs(step) < 2 * c) {
			c = c + step
		}
		else {
			c = (c_high < 1e7 ? 0.5 * (c_low + c_high) : 2 * c)
		}
	}
	return(c)
}

real matrix getW(real matrix distmat, real scalar c0)
{
	external real scalar qmax
	real vector evals
	real matrix W

	W = J(rows(distmat), qmax, 1)
	evals = (1::qmax)
	symeigensystemselecti(demeanmat(exp(-c0 * distmat)), (1, qmax), W, evals)
	W = (J(rows(distmat), 1, 1 / sqrt(rows(distmat))), W)
	return(W)
}

// Fully Vectorized 64-point Gauss-Legendre Quadrature (Zero loop overhead)
real scalar getrp(real matrix Om, real scalar cv)
{
	complex vector cevals
	real vector evals
	real matrix Omx, Ratio
	external real matrix GQxw
	real colvector u, w, one_minus_u2

	Omx = -cv^2 * Om
	Omx[1, .] = Om[1, .]
	cevals = eigenvalues(Omx)
	evals = Re(cevals)
	evals = -select(evals, evals :< 0) / max(evals)

	u = GQxw[., 1]
	w = GQxw[., 2]
	one_minus_u2 = 1.0 :- u:^2
	Ratio = (1.0 :/ one_minus_u2) * evals
	return(sum(w :/ sqrt(one_minus_u2 :* exp(rowsum(log(1.0 :+ Ratio))))) * (2.0 / PI_CONST()))
}

real scalar maxrp(struct mats vector Oms, real scalar q, real scalar cv)
{
	real scalar i, val
	val = 0
	for (i = 1; i <= length(Oms); i++) {
		val = max((val, getrp(Oms[i].mat[1::q+1, 1::q+1], cv)))
	}
	return(val)
}

// High-Performance Brent's Method Root Finder with Safeguarded Bisection
real scalar getcv(struct mats vector Oms, real scalar q, real scalar level)
{
	real scalar i, cv0, cv1, inv_sq, tol, prev_i, k
	real scalar a, b, c, d, e, fa, fb, fc, s, p, q_brent, r_brent, m, iter
	real vector qvec, rps
	real matrix kmat

	inv_sq = 1.0 / sqrt(q)
	tol = 1e-7 * inv_sq // Tightened 7-digit tolerance
	rps = J(length(Oms), 1, 0)
	i = 1
	cv0 = invttail(q, level * 0.5) * inv_sq

	do {
		cv1 = cv0
		while (1) {
			fb = getrp(Oms[i].mat[1::q+1, 1::q+1], cv1) - level
			if (fb > 0) {
				cv0 = cv1
				cv1 = cv1 + inv_sq
			}
			else {
				break
			}
		}

		if (cv1 == cv0) {
			b = cv1
		}
		else {
			a = cv0
			b = cv1
			fa = getrp(Oms[i].mat[1::q+1, 1::q+1], a) - level
			if (abs(fa) < abs(fb)) {
				c = a; a = b; b = c
				fc = fa; fa = fb; fb = fc
			}
			c = a
			fc = fa
			d = b - a
			e = d

			for (iter = 1; iter <= 40; iter++) {
				if (abs(fc) < abs(fb)) {
					a = b; b = c; c = a
					fa = fb; fb = fc; fc = fa
				}
				m = 0.5 * (c - b)
				if (abs(m) <= tol | fb == 0) break

				if (abs(e) >= tol & abs(fa) > abs(fb)) {
					if (a == c) {
						s = fb / fa
						p = 2.0 * m * s
						q_brent = 1.0 - s
					}
					else {
						q_brent = fa / fc
						r_brent = fb / fc
						s = fb / fa
						p = s * (2.0 * m * q_brent * (q_brent - r_brent) - (b - a) * (r_brent - 1.0))
						q_brent = (q_brent - 1.0) * (r_brent - 1.0) * (s - 1.0)
					}
					if (p > 0) q_brent = -q_brent
					p = abs(p)
					if (2.0 * p < min((3.0 * m * q_brent - abs(tol * q_brent), abs(e * q_brent)))) {
						e = d
						d = p / q_brent
					}
					else {
						d = m
						e = m
					}
				}
				else {
					d = m
					e = m
				}

				a = b
				fa = fb
				if (abs(d) > tol) b = b + d
				else b = b + (m > 0 ? tol : -tol)

				fb = getrp(Oms[i].mat[1::q+1, 1::q+1], b) - level

				if ((fb > 0 & fc > 0) | (fb < 0 & fc < 0)) {
					c = a
					fc = fa
					d = b - a
					e = d
				}
			}
		}

		cv1 = b
		for (k = 1; k <= length(Oms); k++) {
			rps[k] = getrp(Oms[k].mat[1::q+1, 1::q+1], cv1)
		}
		maxindex(rps, 1, qvec, kmat)
		prev_i = i
		i = qvec[1]
		if (prev_i == i | rps[i] <= level) {
			break
		}
		cv0 = cv1
	} while (rps[i] > level)

	return(cv1 * sqrt(q))
}

void setfinalW(struct mats vector Oms, real matrix W, real scalar cv)
{
	external real scalar qmax
	real scalar q
	real vector cvs, lengths, qvec
	real matrix kmat

	lengths = (1::qmax)
	cvs = lengths
	for (q = 1; q <= qmax; q++) {
		cvs[q] = getcv(Oms, q, 0.05)
		lengths[q] = cvs[q] * exp(lngamma(0.5 * (q + 1)) - lngamma(0.5 * q)) / sqrt(q)
	}
	minindex(lengths, 1, qvec, kmat)
	q = qvec[1]
	cv = cvs[q]
	W = W[., 1::q+1]
}

real scalar getnc(real scalar c0, real scalar cmax)
{
	external real scalar cgridfac
	return(max((2, ceil(log(cmax / c0) / log(cgridfac)))))
}

struct mats vector getOms(real matrix distmat, real scalar c0, real scalar cmax, real matrix W)
{
	external real scalar cgridfac
	struct mats vector Oms
	real scalar i, c

	Oms = mats(getnc(c0, cmax))
	c = c0
	for (i = 1; i <= length(Oms); i++) {
		Oms[i].mat = cross(W, exp(-c * distmat)) * W
		c = c * cgridfac
	}
	return(Oms)
}

real scalar gettau(real vector y, real matrix W)
{
	return(sqrt(cols(W) - 1) * cross(W[., 1], y) / norm(cross(W[., 2::cols(W)], y)))
}

void setGQxw(real scalar legacy)
{
	external real matrix GQxw
	if (legacy) {
		// 40-point Gauss-Legendre for legacy 2021 reproduction
		GQxw = (  8.811451447204854404E-004,  4.636880650271735238E-003,  1.137002500811307160E-002,  2.104159039310410373E-002,  3.359359586066185122E-002,  4.895059651556293856E-002,  6.702024839387038524E-002,  8.769388458334448355E-002,  1.108471742867402909E-001,  1.363408724050364507E-001,  1.640216576929103831E-001,  1.937230551660097388E-001,  2.252664374524356861E-001,  2.584620991569107629E-001,  2.931103978141975097E-001,  3.290029545871209216E-001,  3.659239074963732685E-001,  4.036512096493146129E-001,  4.419579646623724711E-001,  4.806137912469746754E-001,  5.193862087530253246E-001,  5.580420353376275289E-001,  5.963487903506855536E-001,  6.340760925036268425E-001,  6.709970454128790784E-001,  7.068896021858023238E-001,  7.415379008430891261E-001,  7.747335625475638698E-001,  8.062769448339899281E-001,  8.359783423070896724E-001,  8.636591275949634383E-001,  8.891528257132595980E-001,  9.123061154166557385E-001,  9.329797516061295592E-001,  9.510494034844372280E-001,  9.664064041393380933E-001,  9.789584096068959518E-001,  9.886299749918872060E-001,  9.953631193497285423E-001,  9.991188548552796256E-001 \ 2.260638549266569437E-003,  5.249142265576826200E-003,  8.210529190953967313E-003,  1.112292459708323286E-002,  1.396850349001168591E-002,  1.673009764127359952E-002,  1.939108398723625448E-002,  2.193545409283656489E-002,  2.434790381753590763E-02,  2.661392349196848098E-002,  2.871988454969627222E-002,  3.065312124646369166E-002,  3.240200672830125667E-002,  3.395602290761685904E-002,  3.530582369564324446E-002,  3.644329119790197524E-002,  3.736158452898413751E-002,  3.805518095031324571E-002,  3.851990908212438863E-002,  3.875297398921226377E-002,  3.875297398921258990E-002,  3.851990908212385434E-002,  3.805518095031288489E-002,  3.736158452898429017E-002,  3.644329119790179483E-002,  3.530582369564296691E-002,  3.395602290761668557E-002,  3.240200672830127054E-002,  3.065312124646427452E-002,  2.871988454969589405E-002,  2.661392349196792934E-002,  2.434790381753596661E-002,  2.193545409283669673E-002,  1.939108398723569590E-002,  1.673009764127403667E-002,  1.396850349001163734E-002,  1.112292459708396665E-002,  8.210529190953890985E-003,  5.249142265576283231E-003,  2.260638549266727297E-003)'
	}
	else {
		// High-precision 64-point Gauss-Legendre quadrature
		GQxw = (  3.474791321139147904e-04,  1.829941614022390084e-03,  4.493314261627856521e-03,  8.331873057687011475e-03,  1.333658610504451225e-02,  1.949560017397317102e-02,  2.679431257079861739e-02,  3.521541393403021525e-02,  4.473893146074858773e-02,  5.534227700244293002e-02,  6.700030092295361506e-02,  7.968535187370978656e-02,  9.336734243860123073e-02,  1.080138205283293074e-01,  1.235900463697340301e-01,  1.400590749141945590e-01,  1.573818434728833648e-01,  1.755172643726713444e-01,  1.944223224138033612e-01,  2.140521768986829998e-01,  2.343602679900527175e-01,  2.552984271464735500e-01,  2.768169913732679843e-01,  2.988649210180042148e-01,  3.213899208311659095e-01,  3.443385640048944873e-01,  3.676564188956162882e-01,  3.912881781299964445e-01,  4.151777897880035706e-01,  4.392685903519397139e-01,  4.635034391061004788e-01,  4.878248536682877612e-01,  5.121751463317122388e-01,  5.364965608938995212e-01,  5.607314096480602306e-01,  5.848222102119964294e-01,  6.087118218700034999e-01,  6.323435811043837118e-01,  6.556614359951055127e-01,  6.786100791688340905e-01,  7.011350789819957852e-01,  7.231830086267320157e-01,  7.447015728535264500e-01,  7.656397320099472825e-01,  7.859478231013170557e-01,  8.055776775861966943e-01,  8.244827356273286556e-01,  8.426181565271166907e-01,  8.599409250858054410e-01,  8.764099536302659699e-01,  8.919861794716706926e-01,  9.066326575613987693e-01,  9.203146481262902689e-01,  9.329996990770463849e-01,  9.446577229975570145e-01,  9.552610685392514123e-01,  9.647845860659698403e-01,  9.732056874292014381e-01,  9.805043998260267735e-01,  9.866634138949554877e-01,  9.916681269423129885e-01,  9.955066857383720880e-01,  9.981700583859776099e-01,  9.996525208678861407e-01 \ 8.916403608473938164e-04,  2.073516630281454706e-03,  3.252228984489823668e-03,  4.423379913182191175e-03,  5.584069730065691593e-03,  6.731523948359491305e-03,  7.863015238012549882e-03,  8.975857887848643843e-03,  1.006741157676500562e-02,  1.113508690419152958e-02,  1.217635128435554971e-02,  1.318873485752722167e-02,  1.416983630712979024e-02,  1.511732853620127194e-02,  1.602896417742570578e-02,  1.690258091857070946e-02,  1.773610662844114089e-02,  1.852756427012004770e-02,  1.927507658930775747e-02,  1.997687056636017824e-02,  2.063128162131169438e-02,  2.123675756182679536e-02,  2.179186226466169746e-02,  2.229527908187824842e-02,  2.274581396370902245e-02,  2.314239829065717000e-02,  2.348409140810497203e-02,  2.377008285741511601e-02,  2.399969429822912739e-02,  2.417238111740146328e-02,  2.428773372075171064e-02,  2.434547850456985144e-02,  2.434547850456985144e-02,  2.428773372075171064e-02,  2.417238111740146328e-02,  2.399969429822912739e-02,  2.377008285741511601e-02,  2.348409140810497203e-02,  2.314239829065717000e-02,  2.274581396370902245e-02,  2.229527908187824842e-02,  2.179186226466169746e-02,  2.123675756182679536e-02,  2.063128162131169438e-02,  1.997687056636017824e-02,  1.927507658930775747e-02,  1.852756427012004770e-02,  1.773610662844114089e-02,  1.690258091857070946e-02,  1.602896417742570578e-02,  1.511732853620127194e-02,  1.416983630712979024e-02,  1.318873485752722167e-02,  1.217635128435554971e-02,  1.113508690419152958e-02,  1.006741157676500562e-02,  8.975857887848643843e-03,  7.863015238012549882e-03,  6.731523948359491305e-03,  5.584069730065691593e-03,  4.423379913182191175e-03,  3.252228984489823668e-03,  2.073516630281454706e-03,  8.916403608473938164e-04)'
	}
}

real scalar nextU()
{
	external real scalar legacyprng_flag, random_t
	if (legacyprng_flag) {
		// Legacy 32-bit LCG
		random_t = mod(64389.0 * random_t + 1.0, 4294967296.0)
		return(random_t / 4294967296.0)
	}
	else {
		// Native 64-bit Stata PRNG (respects set seed)
		return(runiform(1, 1))
	}
}

void jumble_s(real matrix s, real scalar m)
{
	external real scalar legacyprng_flag
	real scalar i, j, n
	real vector sj, p
	n = rows(s)
	if (legacyprng_flag) {
		for (i = 1; i <= m; i++) {
			j = floor(nextU() * n) + 1
			sj = s[i, .]
			s[i, .] = s[j, .]
			s[j, .] = sj
		}
	}
	else {
		p = order(runiform(n, 1), 1)
		s[1::m, .] = s[p[1::m], .]
	}
}

real vector raninds(real scalar n, real scalar capM)
{
	external real scalar legacyprng_flag
	real scalar i, j, j0
	real vector v, u, steps

	if (legacyprng_flag) {
		v = J(capM + 1, 1, 0)
		j = floor(n * nextU())
		for (i = 1; i <= capM + 1; i++) {
			v[i] = j + 1
			j = mod(j + 1 + floor(nextU() * (n - 1)), n)
		}
		return(v)
	}
	else {
		j0 = floor(runiform(1, 1) * n)
		u = floor(runiform(capM, 1) * (n - 1))
		steps = (1::capM) :+ runningsum(u)
		v = (j0 + 1) \ (mod(j0 :+ steps, n) :+ 1)
		return(v)
	}
}

void normalize_s(real matrix s)
{
	external real scalar latlongflag
	external real vector permfin
	real scalar n
	n = rows(s)
	if (latlongflag == 0) {
		s = s :- mean(s)
		permfin = order(rowsum(s:^2), 1)
		s = s[permfin, .]
	}
	else {
		s[., 2] = s[., 2] :- mean(s[., 2])
		s[., 2] = mod(s[., 2] :+ 180.0, 360.0) :- 180.0
		permfin = order(s, (1, 2))
		s = s[permfin, .]
	}
}

void LNsetWc0(real matrix s, real scalar avc0, real matrix W, real scalar c0, real scalar cmax)
{
	external real scalar m, capN, qmax, minavc
	struct sms vector ms
	real vector distvec
	real rowvector evals
	real colvector v
	real scalar i, j, n
	real matrix r, W0, Wx, Wall, evecs

	n = rows(s)
	ms = sms(capN)
	distvec = J(capN * m * (m - 1) / 2, 1, 0)

	r = s
	for (i = 1; i <= capN; i++) {
		jumble_s(r, m)
		ms[i].s = r[1::m, .]
		ms[i].mat = getdistmat(ms[i].s)
		distvec[(i - 1) * m * (m - 1) / 2 + 1::i * m * (m - 1) / 2] = lvech(ms[i].mat)
	}
	c0 = getc0_from_vec(distvec, avc0)
	cmax = getc0_from_vec(distvec, minavc)

	Wall = J(n, capN * qmax, 0)
	for (i = 1; i <= capN; i++) {
		symeigensystemselecti(demeanmat(exp(-c0 * ms[i].mat)), (1, qmax), W0, evals)
		Wx = exp(-c0 * getcrossdistmat(s, ms[i].s)) * W0
		Wx = Wx :- mean(Wx)
		Wx = Wx :/ sqrt(colsum(Wx:^2))
		for (j = 1; j <= qmax; j++) {
			Wall[., (j - 1) * capN + i] = Wx[., j]
		}
	}

	W = J(n, qmax, 0)
	for (i = 1; i <= qmax; i++) {
		Wx = Wall[., 1::capN * i]
		symeigensystemselecti(cross(Wx, Wx), (1, 1), evecs, evals)
		W[., i] = Wx * evecs
		W[., i] = W[., i] / sqrt(colsum(W[., i]:^2))
		Wall = Wall - W[., i] * cross(W[., i], Wall)
	}
	W = (J(n, 1, 1.0 / sqrt(n)), W)
}

struct mats vector LNgetOms(real matrix s, real scalar c0, real scalar cmax, real matrix W)
{
	external real scalar capM, cgridfac
	struct mats vector Oms
	real vector inds, dist, cd
	real matrix W1, W2, M, I_W
	real scalar n, nc, i, c

	n = rows(s)
	nc = getnc(c0, cmax)
	Oms = mats(nc)
	inds = raninds(n, capM)
	dist = getdistvec(s[inds[1::capM], .], s[inds[2::capM+1], .])

	W1 = W[inds[1::capM], .]
	W2 = W[inds[2::capM+1], .]
	I_W = I(cols(W))

	c = c0
	for (i = 1; i <= nc; i++) {
		cd = exp(-c * dist)
		// Symmetrized single transpose M + M' (halves 1M-row matrix cross products)
		M = cross(W1, cd, W2)
		Oms[i].mat = I_W + (0.5 * n * (n - 1) / capM) * (M + M')
		c = c * cgridfac
	}
	return(Oms)
}

void set_spatial_scpc(real scalar avc0, string scalar sel, string scalar latlong, string scalar method, real scalar threshold, real scalar cgridfac_arg, real scalar legacyprng_arg)
{
	external real scalar qmax, capM, m, cgridfac, capN, minavc, latlongflag, legacyprng_flag, random_t, c0, cmax
	external struct mats vector Omsfin, Sig_grid
	external real matrix Wfin, distmat, s
	external real vector permfin
	external real rowvector cvs_uncond
	external real scalar cvfin, condflag
	real scalar cv, n, nc, c, i_grid, is_exact
	string scalar slist
	real colvector slist_nums, slist_order
	real scalar slist_len, i
	struct mats vector Oms

	cgridfac = cgridfac_arg
	legacyprng_flag = legacyprng_arg
	random_t = 1.0
	cvs_uncond = J(1, 0, 0) // Reset unconditional CV cache
	Sig_grid = mats(0)

	// slist variable ordering validation
	stata("unab slist: s_*")
	slist = st_local("slist")
	slist = tokens(slist)'
	slist_nums = strtoreal(substr(slist, 3, .))
	slist_order = order(slist_nums, 1)
	slist = slist[slist_order]
	slist_len = length(slist)

	for (i = 1; i <= slist_len; i++) {
		if (substr(slist[i], 3, .) != strofreal(i)) {
			stata("disp as error" + char(34) + "s_* variables not continuously numbered starting from 1 (s_1, s_2, etc.)" + char(34))
			exit(999)
		}
	}

	s = st_data(., slist', sel)
	n = rows(s)
	if (sum(rowmissing(s) :== 0) < n) {
		stata("disp as error" + char(34) + "missing value(s) in s_* variable; aborting" + char(34))
		exit(999)
	}
	if (n < 5) {
		stata("disp as error" + char(34) + "too few locations found in variables s_1, s_2 etc; aborting" + char(34))
		exit(999)
	}

	stata("disp as text" + char(34) + "found " + strofreal(rows(s), "%6.0f") + " observations / clusters and " + strofreal(cols(s), "%3.0f") + "-dimensional locations in s_*" + char(34))
	latlongflag = (latlong != "")
	if (latlongflag == 0) {
		stata("disp as text" + char(34) + "using Euclidean norm to compute distance between locations stored in s_*" + char(34))
	}
	else {
		if (cols(s) != 2) {
			stata("disp as error" + char(34) + "with latlong option, there must only be s_1 and s_2 present" + char(34))
			exit(999)
		}
		else {
			stata("disp as text" + char(34) + "computing distances on surface of sphere treating s_1 as latitude and s_2 as longitude" + char(34))
		}
	}

	condflag = 0
	setGQxw(legacyprng_flag)
	minavc = 0.00001

	if (avc0 >= 0.05) qmax = 10
	else if (avc0 >= 0.01) qmax = 20
	else if (avc0 >= 0.005) qmax = 60
	else qmax = 120

	real matrix W

	while (1) {
		qmax = min((n - 1, qmax))
		is_exact = (method == "exact") | (method == "auto" & n < threshold)

		if (is_exact) {
			// Pure Mata Exact Computation Mode
			permfin = (1::n)
			distmat = getdistmat(s)
			c0 = getc0_from_mat(distmat, avc0)
			cmax = getc0_from_mat(distmat, minavc)
			W = getW(distmat, c0)

			// Precompute and cache Sig_grid once
			nc = getnc(c0, cmax)
			Sig_grid = mats(nc)
			c = c0
			for (i_grid = 1; i_grid <= nc; i_grid++) {
				Sig_grid[i_grid].mat = exp(-c * distmat)
				c = c * cgridfac
			}

			Oms = mats(nc)
			for (i_grid = 1; i_grid <= nc; i_grid++) {
				Oms[i_grid].mat = cross(W, Sig_grid[i_grid].mat) * W
			}
		}
		else {
			// Large-N Nyström approximation mode
			capN = 20
			capM = 1000000
			m = min((1000, n))
			normalize_s(s)
			LNsetWc0(s, avc0, W, c0, cmax)
			Oms = LNgetOms(s, c0, cmax, W)
		}

		setfinalW(Oms, W, cv)
		if (cols(W) - 1 < qmax | qmax == n - 1) break
		qmax = round(qmax + qmax / 2, 1)
	}

	stata("disp as text" + char(34) + "SCPC optimal q = " + strofreal(cols(W) - 1, "%3.0f") + " for maximal average pairwise correlation = " + strofreal(avc0, "%5.3f") + char(34))

	Wfin = W
	Omsfin = Oms
	cvfin = cv
}

void set_scpcstats(string scalar w, string scalar sel)
{
	real vector y
	real scalar tau, pval, SE, cv
	external real scalar condflag
	external struct mats vector Omsfin, Omsxfin
	external real matrix Wfin, Wx
	external real vector permfin
	external real scalar cvfin, cvxfin

	y = st_data(., w, sel)
	y = y[permfin]
	tau = gettau(y, Wfin)
	pval = maxrp(Omsfin, cols(Wfin) - 1, abs(tau) / sqrt(cols(Wfin) - 1))
	cv = cvfin

	if (condflag) {
		pval = max((pval, maxrp(Omsxfin, cols(Wx) - 1, abs(tau) / sqrt(cols(Wx) - 1))))
		cv = max((cv, cvxfin))
	}
	SE = norm(cross(Wfin[., 2::cols(Wfin)], y)) / (sqrt(cols(Wfin) - 1) * sqrt(rows(Wfin)))
	st_matrix("scpcstats", (mean(y), SE, tau, pval, mean(y) - cv * SE, mean(y) + cv * SE))
}

void set_ccspc()
{
	external struct mats vector Omsxfin, Sig_grid
	external real matrix s, Wx
	external real scalar c0, cmax, cvxfin
	real scalar i, use_cached

	use_cached = 0
	if (length(Sig_grid) > 0) {
		if (rows(Wx) == rows(Sig_grid[1].mat)) {
			use_cached = 1
		}
	}

	if (use_cached) {
		// Fast evaluation using cached Sig_grid
		Omsxfin = mats(length(Sig_grid))
		for (i = 1; i <= length(Sig_grid); i++) {
			Omsxfin[i].mat = cross(Wx, Sig_grid[i].mat) * Wx
		}
	}
	else {
		// Large-N evaluation
		Omsxfin = LNgetOms(s, c0, cmax, Wx)
	}
	cvxfin = getcv(Omsxfin, cols(Wx) - 1, 0.05)
}

void set_Wx_nocluster(string scalar Vname, string scalar sel)
{
	external real scalar condflag
	external real matrix Wfin, Wx
	external real vector permfin
	real colvector x, xs
	real matrix X, Y, R, invXX
	string rowvector clist_vars
	real vector r
	real scalar i

	condflag = 1
	x  = st_data(., "scpc_x", sel)[permfin]
	xs = st_data(., "scpc_xs", sel)[permfin]

	// Analytic OLS projection in Mata
	if (substr(strtrim(Vname), 1, 15) == "qui reg scpc_rx") {
		clist_vars = tokens(subinstr(strtrim(Vname), "qui reg scpc_rx", "", 1))
		if (cols(clist_vars) > 0) {
			X = (J(rows(Wfin), 1, 1), st_data(., clist_vars, sel)[permfin, .])
		}
		else {
			X = J(rows(Wfin), 1, 1)
		}
		invXX = invsym(cross(X, X))
		Y = Wfin[., 2::cols(Wfin)] :* xs
		R = Y - X * (invXX * cross(X, Y))
		Wx = (Wfin[., 1] :* xs :* x, R :* x)
	}
	else {
		// Fallback for non-standard models
		r = J(rows(Wfin), 1, 0)
		Wx = Wfin
		for (i = 1; i <= cols(Wfin); i++) {
			r[permfin] = Wfin[., i]
			stata("cap drop scpc_r")
			stata("qui gen scpc_r=.")
			st_store(., "scpc_r", sel, r)
			stata("cap drop scpc_rx")
			if (i == 1) {
				stata("qui gen scpc_rx=scpc_r*scpc_xs*scpc_x if e(sample)")
			}
			else {
				stata("qui gen scpc_rx=scpc_r*scpc_xs if e(sample)")
				stata("_estimates hold oreg, restore")
				stata(Vname)
				stata("cap drop scpc_r")
				stata("qui predict scpc_r, resid")
				stata("_estimates unhold oreg")
				stata("qui replace scpc_rx=scpc_r*scpc_x if e(sample)")
			}
			Wx[., i] = st_data(., "scpc_rx", sel)[permfin]
		}
	}
	set_ccspc()
}

void set_Wx_cluster(string scalar Vname, string scalar sel)
{
	external real scalar condflag
	external real matrix Wfin, Wx
	external real vector permfin
	real vector r
	real scalar i

	condflag = 1
	r = J(rows(Wfin), 1, 0)
	Wx = Wfin
	for (i = 1; i <= cols(Wfin); i++) {
		r[permfin] = Wfin[., i]
		stata("cap drop scpc_r")
		stata("qui gen scpc_r=.")
		stata("qui replace scpc_r=0 if e(sample)")
		st_store(., "scpc_r", sel, r)
		stata("cap drop scpc_rx")
		stata("qui by \`e(clustvar)': egen scpc_rx=total(scpc_r) if e(sample)")
		stata("cap drop scpc_r")
		if (i == 1) {
			stata("qui gen scpc_r=scpc_rx*scpc_xs if e(sample)")
		}
		else {
			stata("qui replace scpc_rx=scpc_rx*scpc_xs")
			stata("_estimates hold oreg, restore")
			stata(Vname)
			stata("qui predict scpc_r, resid")
			stata("_estimates unhold oreg")
		}
		stata("cap drop scpc_rx")
		stata("qui by \`e(clustvar)': egen scpc_rx=total(scpc_r*scpc_x) if e(sample)")
		Wx[., i] = st_data(., "scpc_rx", sel)[permfin]
	}
	set_ccspc()
}

void set_scpccvs()
{
	external struct mats vector Omsfin, Omsxfin
	external real matrix Wfin
	external real scalar condflag, cvxfin
	external real rowvector cvs_uncond
	real scalar q, j
	real matrix levels, cvs

	levels = (.32, .10, .05, .01)
	cvs = J(1, 4, 0)
	q = cols(Wfin) - 1

	// Cache unconditional CVs across covariates
	if (length(cvs_uncond) != 4) {
		cvs_uncond = J(1, 4, 0)
		for (j = 1; j <= 4; j++) {
			cvs_uncond[j] = getcv(Omsfin, q, levels[j])
		}
	}

	for (j = 1; j <= 4; j++) {
		cvs[j] = cvs_uncond[j]
		if (condflag) {
			if (j == 3) {
				// Level 0.05 was already computed in set_ccspc() as cvxfin
				cvs[j] = max((cvs[j], cvxfin))
			}
			else {
				cvs[j] = max((cvs[j], getcv(Omsxfin, q, levels[j])))
			}
		}
	}
	st_matrix("scpc_cvs_c", cvs)
}

end
