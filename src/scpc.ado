version 16.0

set matastrict on
mata 

real colvector scpc_euclidean_row_norm(real matrix differences)
// Scale before squaring to avoid overflowing or underflowing finite distances.
{
	real colvector scale, distance
	if (hasmissing(differences)) {
		errprintf("pairwise coordinate differences exceed the numeric range\n")
		exit(198)
	}
	scale = rowmax(abs(differences))
	scale = scale :+ (scale :== 0)
	distance = scale :* sqrt(rowsum((differences :/ scale):^2))
	if (hasmissing(distance)) {
		errprintf("pairwise distances exceed the numeric range\n")
		exit(198)
	}
	return(distance)
}

real matrix scpc_getdistmat(real matrix s)
// computes matrix of distances from locations
// if external scpc_latlongflag=0, distances computed as norm, otherwise latitude / longitude and haversine formula for sphere with radius 1/Pi
{	
	external real scalar scpc_latlongflag
	real matrix mat, d, rad
	real scalar n, i, c
	n=rows(s)
	mat=J(n,n,0)
	if(scpc_latlongflag==0){
		for(i=1;i<=n;i++){
			mat[.,i]=scpc_euclidean_row_norm(s :- s[i,.])
		}
	}
	else
		{
		c=pi()/180
		for(i=1;i<=n;i++){
			d=(.5*c)*(s:-s[i,.])
			rad=sin(d[.,1]):^2 + cos(c*s[i,1])*(cos(c*s[.,1]):*(sin(d[.,2]):^2))
			rad=rad :* (rad :> 0)
			rad=rad :* (rad :<= 1) :+ (rad :> 1)
			mat[.,i]=asin(sqrt(rad))/pi()
		}
		_diag(mat,0)
	}	
	return(mat)
}	

real vector scpc_lvech(real matrix mat)
// takes lower triangular part of matrix and puts it into vector
{
	real scalar i,c
	real vector v
	v=J(rows(mat)*(rows(mat)-1)/2,1,0)
	c=1
	for(i=1;i<=rows(mat)-1;i++){
		v[c::c+rows(mat)-i-1]=mat[i+1::rows(mat),i]	
		c=c+rows(mat)-i
	}	
	return(v)
}	

real matrix scpc_demeanmat(real matrix mat)
// demeans matrix
{
	mat=mat:-(rowsum(mat)/cols(mat))
	mat=mat:-(colsum(mat)/rows(mat))
	return(mat)
}

void scpc_avc_eval(real vector dist, real scalar c, real scalar target_avc,
                   real scalar near_one, real scalar res, real scalar slope)
{
	real scalar cutoff, n_d, n_big
	real vector d_act, e, m1

	n_d = length(dist)
	cutoff = (c > 0 ? 800 / c : .)
	d_act = select(dist, dist :< cutoff)
	n_big = n_d - length(d_act)
	if (near_one) {
		if (length(d_act) > 0) {
			m1 = expm1(-c * d_act)
			res = (1 - target_avc) + (sum(m1) - n_big) / n_d
			slope = -sum(d_act :* (1 :+ m1)) / n_d
		}
		else {
			res = (1 - target_avc) - n_big / n_d
			slope = 0
		}
	}
	else {
		if (length(d_act) > 0) {
			e = exp(-c * d_act)
			res = sum(e) / n_d - target_avc
			slope = -sum(d_act :* e) / n_d
		}
		else {
			res = -target_avc
			slope = 0
		}
	}
}

real scalar getc0fromavc(real vector dist,real scalar avc0)
// solves for c0 from vector of pairwise distances given avc0
// using co-location floor, scale-invariant initialization, and safeguarded Newton + bisection
{
	real scalar n_d, floor_corr, mean_pos, max_d, target_pos, near_one, raw_init
	real scalar lo, hi, c, f_lo, d_lo, f_hi, d_hi, f, slope, final_f
	real scalar residual_tol, accept_tol, newton, midpoint, next_c, iter
	real vector pos_dist

	n_d = length(dist)
	if (n_d == 0) {
		stata("disp as text " + char(34) + "at least two locations are required; aborting" + char(34))
		exit(999)
	}
	if (any(dist :< 0) | any(dist :>= .)) {
		stata("disp as text " + char(34) + "distances must be finite and nonnegative; aborting" + char(34))
		exit(999)
	}
	if (sum(dist :> 0) == 0) {
		stata("disp as text " + char(34) + "all coordinates are identical; aborting" + char(34))
		exit(999)
	}

	floor_corr = mean(dist :== 0)
	if (avc0 >= . | avc0 <= floor_corr | avc0 >= 1) {
		stata("disp as text " + char(34) + "average correlation is unattainable for these locations; aborting" + char(34))
		exit(999)
	}

	pos_dist = select(dist, dist :> 0)
	mean_pos = mean(pos_dist)
	max_d = max(pos_dist)
	target_pos = (avc0 - floor_corr) / (1 - floor_corr)
	near_one = (avc0 > 0.9999)

	raw_init = -log(target_pos) / mean_pos
	if (raw_init <= 0 | raw_init >= .) {
		raw_init = 1 / max_d
		if (raw_init <= 0 | raw_init >= .) raw_init = 1
	}

	lo = raw_init
	hi = raw_init
	f_lo = 0; d_lo = 0; f_hi = 0; d_hi = 0; f = 0; slope = 0; final_f = 0
	scpc_avc_eval(dist, lo, avc0, near_one, f_lo, d_lo)
	if (f_lo == 0) return(lo)

	if (f_lo < 0) {
		hi = lo
		for (iter = 1; iter <= 1024; iter++) {
			lo = 0.5 * lo
			if (lo == 0 | lo >= .) {
				stata("disp as text " + char(34) + "decay calibration lower bracket underflowed; aborting" + char(34))
				exit(999)
			}
			scpc_avc_eval(dist, lo, avc0, near_one, f_lo, d_lo)
			if (f_lo == 0) return(lo)
			if (f_lo > 0) break
			hi = lo
		}
		if (f_lo < 0) {
			stata("disp as text " + char(34) + "decay calibration lower bracket did not enclose root; aborting" + char(34))
			exit(999)
		}
	}
	else {
		f_hi = f_lo
		for (iter = 1; iter <= 1024; iter++) {
			lo = hi
			f_lo = f_hi
			hi = 2 * hi
			if (hi >= .) {
				stata("disp as text " + char(34) + "spatial decay calibration failed to bracket upper decay scale; aborting" + char(34))
				exit(999)
			}
			scpc_avc_eval(dist, hi, avc0, near_one, f_hi, d_hi)
			if (f_hi == 0) return(hi)
			if (f_hi < 0) break
			d_lo = d_hi
		}
		if (f_hi > 0) {
			stata("disp as text " + char(34) + "spatial decay calibration failed to bracket upper decay scale; aborting" + char(34))
			exit(999)
		}
	}

	residual_tol = 1e-12 * max((min((avc0 - floor_corr, 1 - avc0)), 1e-15))
	c = 0.5 * (lo + hi)
	if (d_lo < 0 & f_lo < .) {
		newton = lo - f_lo / d_lo
		if (newton > lo & newton < hi & newton < .) {
			c = newton
		}
		else if (lo > 0 & hi / lo > 4) {
			c = sqrt(lo * hi)
		}
	}
	else if (lo > 0 & hi / lo > 4) {
		c = sqrt(lo * hi)
	}
	if (c <= lo | c >= hi | c >= .) c = 0.5 * (lo + hi)

	for (iter = 1; iter <= 128; iter++) {
		scpc_avc_eval(dist, c, avc0, near_one, f, slope)
		if (abs(f) <= residual_tol) return(c)
		if (f > 0) lo = c
		else hi = c
		if ((hi - lo) <= 1e-13 * max((c, smallestdouble()))) {
			c = 0.5 * (lo + hi)
			break
		}

		midpoint = lo + 0.5 * (hi - lo)
		if (lo > 0 & hi / lo > 4) {
			next_c = sqrt(lo * hi)
			if (next_c > lo & next_c < hi & next_c < .) midpoint = next_c
		}

		next_c = midpoint
		if (slope < 0 & slope > -. & f < .) {
			newton = c - f / slope
			if (newton > lo & newton < hi & newton < .) next_c = newton
		}

		if (next_c == c) {
			c = midpoint
			if (c == lo | c == hi) break
		}
		else {
			c = next_c
		}
	}

	scpc_avc_eval(dist, c, avc0, near_one, final_f, slope)
	accept_tol = 1e-8 * max((min((avc0 - floor_corr, 1 - avc0)), 1e-15))
	if (abs(final_f) <= accept_tol | (hi - lo) <= 1e-12 * max((c, smallestdouble()))) {
		return(c)
	}
	stata("disp as text " + char(34) + "spatial decay calibration did not converge; aborting" + char(34))
	exit(999)
	return(.)
}

real matrix getW(real matrix distmat, real scalar c0)
// computes first scpc_qmax eigenvectors of demeaned Sigma(c0) and stores them in W, along with column vector of ones
// all columns of W are normalized to length 1
{
	external real scalar scpc_qmax
	real vector	evals
	real matrix W

	W=J(rows(distmat),scpc_qmax,1)
	evals=(1::scpc_qmax)
	symeigensystemselecti(scpc_demeanmat(exp(-c0*distmat)),(1,scpc_qmax),W,evals)
	W=(J(rows(distmat),1,1/sqrt(rows(distmat))),W)
	return(W)
}

struct scpc_mats{
// helper struct for vector of matrices Oms
		real matrix mat
}

real scalar scpc_rp_panel(real rowvector evals, real rowvector logs,
                          real scalar scale, real scalar is_comp,
                          real scalar lo, real scalar hi)
{
	external real matrix scpc_GQxw
	real scalar width, j, r, m, log_cos, s_log, xj, term
	real colvector theta, s_th, c_th, sin2, cos2, denom, prod_v, deficit, vals, extra, inv
	real rowvector x_row

	width = hi - lo
	theta = lo :+ width * scpc_GQxw[., 1]
	s_th = sin(theta)
	c_th = cos(theta)
	sin2 = s_th :* s_th
	cos2 = c_th :* c_th

	if (is_comp) {
		denom = cos2 + scale * sin2
		prod_v = J(rows(sin2), 1, 1)
		deficit = J(rows(sin2), 1, 0)
		m = length(evals)
		for (j = 1; j <= m; j++) {
			extra = (evals[j] / scale) * denom
			inv = 1 :/ (sin2 + extra)
			deficit = deficit + prod_v :* extra :* inv
			prod_v = prod_v :* sin2 :* inv
		}
		vals = deficit :/ ((1 :+ sqrt(prod_v)) :* denom)
	}
	else {
		prod_v = J(rows(cos2), 1, 1)
		m = length(evals)
		for (j = 1; j <= m; j++) {
			prod_v = prod_v :* (cos2 :/ (cos2 :+ evals[j]))
		}
		vals = sqrt(prod_v)
		for (r = 1; r <= rows(cos2); r++) {
			if (cos2[r] <= 0) {
				vals[r] = 0
			}
			else if (prod_v[r] < smallestdouble()) {
				log_cos = ln(cos2[r])
				x_row = logs :- log_cos
				s_log = 0
				for (j = 1; j <= m; j++) {
					xj = x_row[j]
					if (xj > 800) {
						term = xj
					}
					else if (xj < -800) {
						term = 0
					}
					else if (xj > 0) {
						term = xj + ln1p(exp(-xj))
					}
					else {
						term = ln1p(exp(xj))
					}
					s_log = s_log + term
				}
				vals[r] = (-0.5 * s_log < -800 ? 0 : exp(-0.5 * s_log))
			}
		}
	}
	return(width * sum(scpc_GQxw[., 2] :* vals))
}

real scalar scpc_rp_adaptive(real rowvector evals, real rowvector logs,
                             real scalar scale, real scalar is_comp,
                             real scalar lo, real scalar hi,
                             real scalar whole, real scalar tol, real scalar depth)
{
	real scalar mid, left, right, bound

	mid = 0.5 * (lo + hi)
	left = scpc_rp_panel(evals, logs, scale, is_comp, lo, mid)
	right = scpc_rp_panel(evals, logs, scale, is_comp, mid, hi)
	bound = max((tol, 1e-11 * abs(whole)))
	if (depth == 20 | abs(left + right - whole) <= bound) {
		return(left + right)
	}
	return(scpc_rp_adaptive(evals, logs, scale, is_comp, lo, mid, left, 0.5 * tol, depth + 1) +
	       scpc_rp_adaptive(evals, logs, scale, is_comp, mid, hi, right, 0.5 * tol, depth + 1))
}

// Mata's power operator rounds nontrivial powers of two via exp/log.
real scalar scpc_power_two(real scalar exponent)
{
    real scalar factor, remaining, value
    factor = exponent < 0 ? 0.5 : 2
    remaining = abs(exponent)
    value = 1
    while (remaining > 0) {
        if (mod(remaining,2)) value = value * factor
        remaining = floor(remaining/2)
        if (remaining > 0) factor = factor * factor
    }
    return(value)
}

// Dekker's error-free product for power-of-two scaled inputs.
real scalar scpc_product_tail(real scalar x, real scalar y, real scalar p)
{
    real scalar sx, sy, xh, yh, xl, yl
    sx = 134217729 * x
    sy = 134217729 * y
    xh = sx - (sx - x)
    yh = sy - (sy - y)
    xl = x - xh
    yl = y - yh
    return(xl*yl - (((p - xh*yh) - xl*yh) - xh*yl))
}

real scalar scpc_bivariate_rejection(real scalar a, real scalar b, real scalar d, real scalar cv)
{
	// Strict Gaussian quadrant probability for (X - cv Y)(X + cv Y) > 0.
	real scalar tol, magnitude, v0, v1, ea, ed, eb, an, dn, bn, ad, bb, det, t, td, ttd, spread, difference, angle
	tol = 64 * epsilon(1)
	magnitude = max(abs((a,b,d)))
	if (magnitude == 0) return(0)
	v0 = a/magnitude + tol
	v1 = d/magnitude + tol
	if (min((v0,v1)) < 0 | abs(b)/magnitude > sqrt(v0)*sqrt(v1)) {
		_error(3498, "SCPC covariance is not positive semidefinite")
	}
	a = max((0,a))
	d = max((0,d))
	if (min((a,d)) == 0) return(a > 0)
	// Evaluate atan2(2*c*sqrt(a*d-b^2), c^2*d-a) in scaled units.
	ea = max((-1022,min((1022,floor(ln(a)/ln(2))))))
	ed = max((-1022,min((1022,floor(ln(d)/ln(2))))))
	eb = floor((ea+ed)/2)
	an = a / scpc_power_two(ea)
	dn = (d / scpc_power_two(ed)) * (1 + ea + ed - 2*eb)
	bn = b / scpc_power_two(eb)
	ad = an*dn
	bb = bn*bn
	if (bb > 2*ad) det = 0
	else det = max((0, (ad-bb) + (scpc_product_tail(an,dn,ad) - scpc_product_tail(bn,bn,bb))))
	t = cv * scpc_power_two(eb-ea)
	if (t >= .) return(0)
	if (t > 2^500) {
		spread = 2 * sqrt(det) / t
		difference = dn
	}
	else {
		td = t*dn
		ttd = t*td
		spread = 2*t*sqrt(det)
		difference = (ttd-an) + (t*scpc_product_tail(t,dn,td) + scpc_product_tail(td,t,ttd))
	}
	if (spread == 0) return(difference < 0)
	if (difference == 0) return(0.5)
	if (abs(difference) > spread) angle = atan(spread / abs(difference)) / pi()
	else angle = 0.5 - atan(abs(difference) / spread) / pi()
	return(difference > 0 ? angle : 1 - angle)
}

real scalar scpc_one_factor_rejection(real matrix Om, real scalar cv)
{
    real matrix A
    real scalar scale, i, j, a, b, c, d, p, denominator
    if (Om[1,1] <= 0 | Om[1,1]*Om[2,2] != Om[2,1]*Om[2,1]) return(.)
    scale = floor(ln(max(diagonal(Om))) / ln(2))
    scale = scpc_power_two(max((-1022, min((1022, scale)))))
    A = Om / scale
    if (any(A*scale :!= Om) | any(Om :!= Om')) return(.)
    a = A[1,1]
    for (i=2; i<=rows(A); i++) {
        for (j=2; j<=i; j++) {
            b = A[i,j]
            c = A[i,1]
            d = A[j,1]
            p = a*b
            if (p != c*d) return(.)
            if (p == 0 & (b != 0 | (c != 0 & d != 0))) return(.)
            if (p != 0 & abs(p) < 2^(-918)) return(.)
            if (scpc_product_tail(a,b,p) != scpc_product_tail(c,d,p)) return(.)
        }
    }
    denominator = sum(diagonal(A[2..rows(A),2..rows(A)]))
    return(cv > 1 ? (a/cv)/cv > denominator : a > (denominator*cv)*cv)
}

real rowvector scpc_rejection_eigenvalues(real matrix Om, real scalar cv, real scalar tol)
{
	real scalar n, i, repair
	real colvector v, sd, divisor
	real matrix C, L, Q, W
	real rowvector vals, keep
	n = rows(Om)
	v = diagonal(Om)
	if (min(v) < -tol) _error(3498, "SCPC covariance is not positive semidefinite")
	sd = J(n, 1, sqrt(max((0, max(v[2..n])))))
	sd[1] = sqrt(max((0, v[1])))
	divisor = sd + (sd :== 0)
	C = (Om :/ divisor) :/ divisor'
	repair = 0
	for (i=1; i<=n; i++) {
		if (sd[i] == 0) {
			if (max(abs(Om[i,.])) > tol) repair = 1
			C[i,.] = J(1,n,0)
			C[.,i] = J(n,1,0)
		}
	}
	if (repair | hasmissing(C)) {
		if (max(v) <= 0) _error(3498, "SCPC covariance is not positive semidefinite")
		C = Om / max(v)
		sd = J(n,1,sqrt(max(v)))
	}
	L = cholesky(C)
	if (hasmissing(L)) {
		symeigensystem(C, Q, vals)
		if (min(vals) < -32*n*epsilon(1)) {
			if (max(v) <= 0) _error(3498, "SCPC covariance is not positive semidefinite")
			symeigensystem(Om / max(v), Q, vals)
			if (min(vals) < -32*n*epsilon(1)) _error(3498, "SCPC covariance is not positive semidefinite")
			sd = J(n,1,sqrt(max(v)))
		}
		keep = selectindex(vals :> 0)
		if (!length(keep)) return(0)
		L = Q[.,keep] :* sqrt(vals[keep])
	}
	if (cv > 1) sd[1] = sd[1] / cv
	else sd[2..n] = sd[2..n] * cv
	if (max(sd) == 0) return(0)
	W = L :* (sd / max(sd))
	return(symeigenvalues(cross(W[1,.], W[1,.]) - cross(W[2..n,.], W[2..n,.])))
}

real scalar getrp(real matrix Om, real scalar cv)
{
// computes rejection probability given Om and cv; cv here is cv/sqrt(q) in notation of paper
	external real matrix scpc_GQxw
	real rowvector vals, evals, logs
	real colvector active, omitted
	real scalar positive, scale, half_pi, whole, tol, prob, b, j

	if (rows(scpc_GQxw) == 0) setGQxw()
	if (cv >= .) _error(3498, "SCPC critical value must be finite")
	if (cv < 0) return(1)
	if (rows(Om) != cols(Om) | rows(Om) < 2 | hasmissing(Om)) {
		_error(3498, "SCPC covariance must be finite and square with at least two rows")
	}
	tol = 32 * rows(Om) * epsilon(1) * max((0, max(diagonal(Om))))
	if (max(abs(Om - Om')) > tol) _error(3498, "SCPC covariance is not symmetric")
	if (rows(Om) > 2) {
		active = selectindex(diagonal(Om)[2..rows(Om)] :!= 0) :+ 1
		if (length(active) <= 1) {
			j = (length(active) ? active[1] : 2)
			omitted = select((2::rows(Om)), (2::rows(Om)) :!= j)
			if (!any(Om[omitted,.]) & !any(Om[.,omitted])) {
				b = Om[1,j] + 0.5*(Om[j,1]-Om[1,j])
				return(scpc_bivariate_rejection(Om[1,1], b, Om[j,j], cv))
			}
		}
	}
	if (rows(Om) == 2) {
		b = Om[1,2] + 0.5 * (Om[2,1] - Om[1,2])
		return(scpc_bivariate_rejection(Om[1,1], b, Om[2,2], cv))
	}
	prob = scpc_one_factor_rejection(Om, cv)
	if (prob < .) return(prob)
	vals = scpc_rejection_eigenvalues(Om, cv, tol)
	positive = max(vals)
	if (positive <= 0 | positive >= .) return(0)
	evals = -select(vals, vals :< 0) / positive
	evals = select(evals, evals :> 0)
	if (length(evals) == 0) return(1)
	if (length(evals) == 1) return((2 / pi()) * atan(1 / sqrt(evals[1])))

	logs = ln(evals)
	scale = sum(evals)
	if (scale == 0 | scale >= .) return(1)
	half_pi = 0.5 * pi()
	if (scale < 0.0625) {
		whole = scpc_rp_panel(evals, logs, scale, 1, 0, half_pi)
		tol = (abs(whole) < 1e-8 ? 1e-10 * abs(whole) : 1e-11)
		prob = 1 - sqrt(scale) * (scpc_rp_adaptive(evals, logs, scale, 1, 0, half_pi, whole, tol, 0) * (2 / pi()))
	}
	else {
		whole = scpc_rp_panel(evals, logs, scale, 0, 0, half_pi)
		tol = (abs(whole) < 1e-8 ? 1e-10 * abs(whole) : 1e-11)
		prob = scpc_rp_adaptive(evals, logs, scale, 0, 0, half_pi, whole, tol, 0) * (2 / pi())
	}
	if (prob < 0) prob = 0
	if (prob > 1) prob = 1
	return(prob)
}

real scalar maxrp(struct scpc_mats vector Oms, real scalar q, real scalar cv)
{
// computes largest rejection probability given vector of Om matrices stored in Oms given q and cv
	real scalar i, val
	val=0
	for(i=1;i<=length(Oms);i++){
		val=max((val,getrp(Oms[i].mat[1::q+1,1::q+1],cv)))
	}
	return(val)
}

real scalar getcv(struct scpc_mats vector Oms, real scalar q, real scalar level)
// computes (two-sided) critical value from q and Oms of given level
{
	real scalar active, j, outer, iter, step, tol, cv0, cv1, mid, rp
	real vector qvec, rps
	real matrix kmat

	rps = J(length(Oms), 1, 0)
	active = 1			// compute cv1 first for Oms[1]; if that doesn't lead to size control, iterate
	step = 1 / sqrt(q)
	tol = 0.001 * step
	cv0 = invttail(q, level / 2) * step
	cv1 = cv0
	rp = 1

	for (outer = 1; outer <= length(Oms) + 5; outer++) {
		if (rp <= level) break
		cv1 = cv0
		for (iter = 0; iter < 1024; iter++) {
			if (getrp(Oms[active].mat[1::q+1, 1::q+1], cv1) > level) {
				cv0 = cv1
				cv1 = cv1 + (iter < 32 ? step : cv1)
				if (cv1 >= .) {
					stata("disp as text " + char(34) + "SCPC critical value could not be bracketed; aborting" + char(34))
					exit(999)
				}
			}
			else {
				break
			}
		}
		if (getrp(Oms[active].mat[1::q+1, 1::q+1], cv1) > level) {
			stata("disp as text " + char(34) + "SCPC critical value could not be bracketed; aborting" + char(34))
			exit(999)
		}

		while (cv1 - cv0 > tol) {
			mid = 0.5 * (cv1 + cv0)
			if (mid == cv0 | mid == cv1) break
			if (getrp(Oms[active].mat[1::q+1, 1::q+1], mid) > level) {
				cv0 = mid
			}
			else {
				cv1 = mid
			}
		}

		for (j = 1; j <= length(Oms); j++) {
			rps[j] = getrp(Oms[j].mat[1::q+1, 1::q+1], cv1)
		}
		maxindex(rps, 1, qvec, kmat)
		if (qvec[1] == active) {
			break
		}
		active = qvec[1]
		rp = rps[active]
		cv0 = cv1
	}
	return(cv1 * sqrt(q))
}

void setfinalW(struct scpc_mats vector Oms, real matrix W, real scalar cv)
// solves for optimal q and cv from Oms, stores results in W and cv
{
	external real scalar scpc_qmax
	real scalar q
	real vector cvs,lengths,qvec
	real matrix kmat
	
	lengths=(1::scpc_qmax)
	cvs=lengths
	for(q=1;q<=scpc_qmax;q++){
		cvs[q]=getcv(Oms,q,0.05)
		lengths[q]=cvs[q]*exp(lngamma(.5*(q+1))-lngamma(.5*q))/sqrt(q)
	} 
	minindex(lengths,1,qvec,kmat)
	q=qvec[1]
	cv=cvs[q]	
	W=W[.,1::q+1]
}

real rowvector scpc_decay_grid(real scalar c0, real scalar cmax)
{
	external real scalar scpc_cgridfac
	real scalar ratio, span, n, i, c
	real rowvector grid
	ratio = cmax / c0
	span = (ratio < . ? log(ratio) : log(cmax) - log(c0))
	n = max((2, ceil(span / log(scpc_cgridfac))))
	grid = J(1, n, .)
	c = c0
	for (i = 1; i <= n & c < .; i++) {
		grid[i] = c
		c = c * scpc_cgridfac
	}
	grid = grid[1..(i-1)]
	if (cmax < . & grid[cols(grid)] < cmax) grid = (grid, cmax)
	return(grid)
}

struct scpc_mats vector getOms(real matrix distmat, real scalar c0, real scalar cmax, real matrix W)
// computes vector of scpc_qmax x scpc_qmax Om(c) matrices from distmat, c0 and W
{
	struct scpc_mats vector Oms
	real scalar i
	real rowvector grid

	grid = scpc_decay_grid(c0, cmax)
	Oms=scpc_mats(cols(grid))
	for(i=1;i<=length(Oms);i++){	
		Oms[i].mat=cross(W,exp(-grid[i]*distmat))*W
	}
	return(Oms)
}

void setGQxw()
{
// sets 40-point Gauss-Legendre quadrature nodes and weights on [0, 1]
	external real matrix scpc_GQxw
	scpc_GQxw=(  0.0008811451447203744,  0.004636880650271513,  0.01137002500811285,  0.02104159039310416,  0.03359359586066174,  0.04895059651556283,  0.06702024839387027,  0.08769388458334415,  0.11084717428674029,  0.13634087240503645,  0.16402165769291022,  0.1937230551660099,  0.2252664374524359,  0.25846209915691065,  0.2931103978141975,  0.32900295458712076,  0.36592390749637316,  0.40365120964931445,  0.4419579646623724,  0.48061379124697456,  0.5193862087530254,  0.5580420353376276,  0.5963487903506856,  0.6340760925036268,  0.6709970454128793,  0.7068896021858024,  0.7415379008430893,  0.7747335625475641,  0.8062769448339902,  0.8359783423070898,  0.8636591275949635,  0.8891528257132597,  0.9123061154166558,  0.9329797516061298,  0.9510494034844372,  0.9664064041393383,  0.9789584096068958,  0.9886299749918872,  0.9953631193497285,  0.9991188548552796 \ 0.0022606385492667702,  0.005249142265577026,  0.008210529190953473,  0.011122924597083375,  0.013968503490011764,  0.016730097641273825,  0.019391083987236147,  0.0219354540928367,  0.024347903817536057,  0.026613923491968422,  0.02871988454969573,  0.030653121246464486,  0.0324020067283005,  0.0339560229076169,  0.03530582369564337,  0.036443291197901954,  0.03736158452898409,  0.03805518095031309,  0.03851990908212394,  0.038752973989212375,  0.038752973989212375,  0.03851990908212394,  0.03805518095031309,  0.03736158452898409,  0.036443291197901954,  0.03530582369564337,  0.0339560229076169,  0.0324020067283005,  0.030653121246464486,  0.02871988454969573,  0.026613923491968422,  0.024347903817536057,  0.0219354540928367,  0.019391083987236147,  0.016730097641273825,  0.013968503490011764,  0.011122924597083375,  0.008210529190953473,  0.005249142265577026,  0.0022606385492667702)'
}

/* 
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXX now routines for n>2000 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  */

real matrix getdistvec(real matrix s1, real matrix s2)
// computes vector of pariwise location distances
// if external scpc_latlongflag=0, distances computed as norm, otherwise latitude / longitude and haversine formula for sphere with radius 1/Pi
{
	external real scalar scpc_latlongflag
	real vector dist, rad
	real scalar c
	real matrix d
	if(scpc_latlongflag==0){
		dist=scpc_euclidean_row_norm(s1 - s2)
	}
	else
		{
		c=pi()/180
		d=(.5*c)*(s1-s2)
		rad=sin(d[.,1]):^2 + cos(c*s1[.,1]):*(cos(c*s2[.,1]):*(sin(d[.,2]):^2))
		rad=rad :* (rad :> 0)
		rad=rad :* (rad :<= 1) :+ (rad :> 1)
		dist=asin(sqrt(rad))/pi()
	}	
	return(dist)
}	

void normalize_s(real matrix s)
// normalizes locations to ensure invariance given probabilistic algorithm
// if locations are in lat/long format, just normalize longitudes (so invariance of rotations around polar axis)
{
	external real scalar scpc_latlongflag
	real matrix evecs
	real vector evals
	real scalar scale
	external real vector scpc_permfin
	
	if(scpc_latlongflag==0) {
		s=s:-(0.5 * colmin(s) + 0.5 * colmax(s))
		scale=max(abs(s))
		if(scale>0) s=s/scale
		s=s:-mean(s)
		symeigensystem(cross(s,s),evecs,evals)
		s=s*evecs		
		if(max(s[.,1])!=max(abs(s[.,1]))) s=-s
		scpc_permfin=order(s,(1::cols(s))')
		s=s[scpc_permfin,.]
		s=s:-colmin(s)
		if(max(s)>0) s=s/max(s)
	}
	else
	{
		s[.,2]=s[.,2]:-mean(s[.,2])
		s[.,2]=mod(s[.,2]:+180,360):-180
		scpc_permfin=order(s,(1,2))
		s=s[scpc_permfin,.]
	}
}

real scalar nextU()
// simple pseudo-random number generator, state (=seed) in global variable scpc_random_t
{
	external real scalar scpc_random_t
	scpc_random_t=mod(64389*scpc_random_t + 1,2^32)
	return(scpc_random_t/2^32)
}

void jumble_s(real matrix s, real scalar m)
// randomly jumbles first m rows of location matrix s without replacement
{
	real scalar i,j,n
	real vector sj
	n=rows(s)	
	for(i=1;i<=m;i++){
		j=i+floor(nextU()*(n-i+1))
		if (j != i) {
			sj=s[j,.]
			s[j,.]=s[i,.]
			s[i,.]=sj
		}
	}
}	

struct scpc_sms{
// helper struct for capN random sets of m locations: mat stores distmat, s stores locations
		real matrix mat
		real matrix s
}

void LNsetWc0(real matrix s,real scalar avc,real matrix W, real scalar c0, real scalar cmax)
{
// computes W and c0 when n is large
	external real scalar scpc_m,scpc_capN
	external real scalar scpc_qmax, scpc_minavc
	struct scpc_sms vector ms
	real vector distvec
	real rowvector evals
	real colvector v
	real scalar i,j,n,floor_corr,minavc_eff
	real matrix r,W0,scpc_Wx,Wall,evecs
	
	n=rows(s)
	ms=scpc_sms(scpc_capN)
	distvec=J(scpc_capN*scpc_m*(scpc_m-1)/2,1,0) // vector of all pairwise distances across scpc_capN random locations
	
	r=s	// ensure s does not get jumbled
	for(i=1;i<=scpc_capN;i++){
		jumble_s(r,scpc_m)
		ms[i].s=r[1::scpc_m,.]
		ms[i].mat=scpc_getdistmat(ms[i].s)
		distvec[(i-1)*scpc_m*(scpc_m-1)/2+1::i*scpc_m*(scpc_m-1)/2]=scpc_lvech(ms[i].mat)
	}
	floor_corr=mean(distvec :== 0)
	minavc_eff=floor_corr + (1 - floor_corr) * scpc_minavc
	c0=getc0fromavc(distvec,avc)
	cmax=getc0fromavc(distvec,minavc_eff)
	
	Wall=J(n,scpc_capN*scpc_qmax,0)		// matrix of all eigenvectors 1,...,scpc_qmax across all scpc_capN random locations
	for(i=1;i<=scpc_capN;i++){
		symeigensystemselecti(scpc_demeanmat(exp(-c0*ms[i].mat)),(1,scpc_qmax),W0,evals)
		scpc_Wx=J(n,scpc_qmax,0)		// m x scpc_qmax eigenvectors from ith random set of m locations
		for(j=1;j<=scpc_m;j++){
			v=exp(-c0*getdistvec(s, J(n, 1, 1)*ms[i].s[j,.]))
			scpc_Wx=scpc_Wx+v*W0[j,.]		// Nystrom extension
		}
		scpc_Wx=scpc_Wx:-mean(scpc_Wx)
		scpc_Wx=scpc_Wx:/sqrt(colsum(scpc_Wx:^2))
		for(j=1;j<=scpc_qmax;j++){
			Wall[.,(j-1)*scpc_capN+i]=scpc_Wx[.,j] // store such that first capN vectors correspond to eigenvectors corresponding to largest eigenvalue, etc
		}
	}	
	W=J(n,scpc_qmax,0)
	for(i=1;i<=scpc_qmax;i++){
		scpc_Wx=Wall[.,1::scpc_capN*i]
		symeigensystemselecti(cross(scpc_Wx,scpc_Wx),(1,1),evecs,evals)
		W[.,i]=scpc_Wx*evecs		// first (sample) principle component of first 1,...,q*capN columns of Wall
		W[.,i]=W[.,i]/sqrt(colsum(W[.,i]:^2))
		Wall=Wall-W[.,i]*cross(W[.,i],Wall) // project off new principle component
	}
	W=(J(n,1,1/sqrt(n)),W)
}

real vector raninds(real scalar n,real scalar capM)
{
// yields capM+1 vector of random indices in {1,...,n} such that each pair is distinct
 	real scalar i,j
	real vector v

	v=J(capM+1,1,0)
	j=floor(n*nextU())
	for(i=1;i<=capM+1;i++){
		v[i]=j+1
		j=mod(j+1+floor(nextU()*(n-1)),n)
	}
	return(v)
}

struct scpc_pair_sample {
    real colvector inds
    real colvector dist
}

struct scpc_pair_sample scalar scpc_draw_pairs(real matrix s)
{
    external real scalar scpc_capM
    struct scpc_pair_sample scalar pairs

    pairs.inds = raninds(rows(s), scpc_capM)
    pairs.dist = getdistvec(s[pairs.inds[1::scpc_capM], .], s[pairs.inds[2::scpc_capM+1], .])
    return(pairs)
}

struct scpc_mats vector LNgetOms(real scalar c0, real scalar cmax, real matrix W,
    struct scpc_pair_sample scalar pairs)
{
	struct scpc_mats vector Oms
	real scalar i, n, capM, covariance_scale
	real rowvector grid
	real vector cd, eigenvalues
	real matrix W1,W2,gram, cross_sample, normalized_covariance

	n=rows(W)
	capM=rows(pairs.dist)
	grid = scpc_decay_grid(c0, cmax)
	Oms=scpc_mats(cols(grid))
	W1=W[pairs.inds[1::capM],.]
	W2=W[pairs.inds[2::capM+1],.]
	gram=cross(W, W)
	for(i=1;i<=length(Oms);i++){
		cd=exp(-grid[i]*pairs.dist)
		cross_sample = cross(W1, cd, W2)
		Oms[i].mat = gram + .5 * (n * (n - 1) / capM) * (cross_sample + cross_sample')
		if (hasmissing(Oms[i].mat)) {
			errprintf("invalid sampled covariance: SCPC covariance grid is nonfinite; try exact SCPC\n")
			exit(198)
		}
		covariance_scale = max(abs(Oms[i].mat))
		if (covariance_scale > 0) {
			normalized_covariance = Oms[i].mat / covariance_scale
			eigenvalues = symeigenvalues(normalized_covariance)
			if (hasmissing(eigenvalues) | min(eigenvalues) < -32 * cols(W) * epsilon(1)) {
				errprintf("invalid sampled covariance: SCPC covariance grid is not positive semidefinite; try exact SCPC\n")
				exit(198)
			}
		}
	}
	return(Oms)
}


/* 
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXX driver routines XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  */

void set_method_opts(real scalar is_approx, real scalar seed_in)
{
	external real scalar scpc_approxflag, scpc_seedval
	scpc_approxflag = (is_approx != 0)
	scpc_seedval = seed_in
}

void setOmsWfin(real scalar avc0, string scalar sel, string scalar latlong)
{
	external real matrix scpc_geometry
	external struct scpc_pair_sample scalar scpc_cond_pairs
	external real scalar scpc_c0,scpc_cmax
	real scalar cv,n,floor_corr,minavc_eff,is_approx,seed_use,q_limit
	external real matrix scpc_distmat
	real matrix W
	external real scalar scpc_random_t, scpc_approxflag, scpc_seedval
	external real scalar scpc_latlongflag
	struct scpc_mats vector Oms
	external real scalar scpc_qmax, scpc_capM, scpc_m, scpc_cgridfac, scpc_capN, scpc_minavc

	external struct scpc_mats vector scpc_Omsfin
	external real matrix scpc_Wfin
	external real vector scpc_permfin
	external real scalar scpc_cvfin
	external real scalar scpc_condflag
	string scalar slist_raw
	string vector slist
	real vector slist_nums, slist_order, dvec
	real scalar slist_len, i
	
	// BUG FIX WITH COORDINATE VARIABLE ORDERING, thank you David Boll!
	
	stata("capture unab slist: s_*")  // create local macro slist containing a list of all s_* variables

	slist_raw = st_local("slist") // import to mata
	slist = tokens(slist_raw)' // tokenize list into vector
	slist_len = length(slist) // number of elements
	if (slist_len == 0) {
		stata("disp as text"+char(34)+"s_* variables not continuously numbered starting from 1 (s_1, s_2, etc.)"+char(34))
		exit(999)
	}
	slist_nums = strtoreal(substr(slist,3,.)) // extract numbers for sorting
	slist_order = order(slist_nums,1) // extract numerical order (alphabetical order would be incorrect for d>9)
	slist = slist[slist_order] // re-order elements of slist

	for (i=1; i<=slist_len; i++) {
		if (substr(slist[i],3,.) != strofreal(i)){ 	// check that i'th element of slist is equal to "s_i", prohibiting e.g. "s_2 s_3" or "s_1 s_3"
			stata("disp as text"+char(34)+"s_* variables not continuously numbered starting from 1 (s_1, s_2, etc.)"+char(34))
			exit(999)
		}
	}
	
	scpc_geometry=st_data(.,slist',sel)		// Import slist variables, now correctly ordered and continuously numbered
								// expects locations in s_1, s_2, s_3... etc. Dimension d is equal to the number of s_ variables 
	
	// ----------------------------------------------------------------------
	
	n=rows(scpc_geometry)
	if(sum(rowmissing(scpc_geometry) :== 0) < n){
		stata("disp as text"+char(34)+"missing value(s) in s_* variable; aborting"+char(34))
		exit(999)
	}
	if(n<5){
		stata("disp as text"+char(34)+"too few locations found in variables s_1, s_2 etc; aborting"+char(34))
		exit(999)
	}		
	stata("disp as text"+char(34)+"found "+strofreal(rows(scpc_geometry), "%6.0f")+" observations / clusters and "+strofreal(cols(scpc_geometry), "%3.0f")+"-dimensional locations in s_*"+char(34))
	scpc_latlongflag=latlong!=""
	if(scpc_latlongflag==0) {
		stata("disp as text"+char(34)+"using Euclidean norm to compute distance between locations stored in s_*")
	}
	else
	{
		if(cols(scpc_geometry)!=2) {
			stata("disp as text"+char(34)+"with latlong option, there must only be s_1 and s_2 present")
			exit(999)
		}
		else{
			stata("disp as text"+char(34)+"Computing distances on surface of sphere treating s_1 as latitude and s_2 as longitude")
		}
	}
	scpc_condflag=0
	scpc_cond_pairs=scpc_pair_sample()
	setGQxw()
	scpc_minavc=0.00001			// minimal avc value for which size control is checked
	scpc_cgridfac=1.2			// factor in c-grid for size control
	
	if(avc0>=0.05){
		scpc_qmax=10
	}
	else{
		if(avc0>=0.01){
		scpc_qmax=20
		}
		else{
			if(avc0>=0.005){
				scpc_qmax=60
			}
			else{
				scpc_qmax=120
			}
		}
	}
	is_approx = (length(scpc_approxflag) > 0 ? (scpc_approxflag != 0) : (n >= 4500))
	seed_use = (length(scpc_seedval) > 0 ? scpc_seedval : 1)
	if (is_approx == 0) {
		scpc_permfin = (1::n)
		scpc_distmat = scpc_getdistmat(scpc_geometry)
		dvec = scpc_lvech(scpc_distmat)
		floor_corr = mean(dvec :== 0)
		minavc_eff = floor_corr + (1 - floor_corr) * scpc_minavc
		scpc_c0 = getc0fromavc(dvec, avc0)
		scpc_cmax = getc0fromavc(dvec, minavc_eff)
		q_limit = n - 1
	}
	else {
		scpc_capN=20			// number of random sets of locations (of size scpc_m) used to approximate eigenvectors
		scpc_capM=1000000		// number of location pairs used to approximate Om(c) (in production: scpc_capM=1000000 or so)
		scpc_m=min((1000, n))	// size of set of locations for which eigenvector is computed, and then Nystrom extended from scpc_m to n
		normalize_s(scpc_geometry)
		q_limit = scpc_m - 1
	}
	while(1){
		scpc_qmax=min((q_limit,scpc_qmax))
		if(is_approx==0){
	// code for small n
			W=getW(scpc_distmat,scpc_c0)
			Oms=getOms(scpc_distmat,scpc_c0,scpc_cmax,W)
			}
		else{
			scpc_random_t=seed_use
			LNsetWc0(scpc_geometry,avc0,W,scpc_c0,scpc_cmax)
			Oms=LNgetOms(scpc_c0,scpc_cmax,W,scpc_draw_pairs(scpc_geometry))
		}
		setfinalW(Oms,W,cv)
		if(cols(W)-1<scpc_qmax | scpc_qmax==q_limit) break
		scpc_qmax=ceil(1.5*scpc_qmax)
	}
	stata("disp as text"+char(34)+"SCPC optimal q = "+strofreal(cols(W)-1, "%3.0f")+" for maximal average pairwise correlation = "+strofreal(avc0, "%5.3f")+`pavc'"")
	
	scpc_Wfin=W
	scpc_Omsfin=Oms
	scpc_cvfin=cv

}

	
void set_scpcstats(string scalar w, string scalar sel, string scalar outmat, real scalar beta)
{	
	real vector y
	real scalar tau, pval, SE, y_norm, proj_norm, scale, root_q
	external real scalar scpc_condflag
	external struct scpc_mats vector scpc_Omsfin, scpc_Omsxfin
	external real matrix scpc_Wfin, scpc_Wx
	external real vector scpc_permfin
	external real scalar scpc_cvfin,scpc_cvxfin
	real scalar cv
	

	y=st_data(.,w,sel)

	y=y[scpc_permfin]
	scale=max(abs(y))
	if (scale > 0 & scale < .) y = y / scale
	y = y :- mean(y)
	y_norm=norm(y)
	proj_norm=norm(cross(scpc_Wfin[.,2::cols(scpc_Wfin)],y))
	if(y_norm == 0 | y_norm >= . | proj_norm >= . | proj_norm <= 1e-12 * y_norm){
		stata("disp as text "+char(34)+"SCPC variance projection is degenerate; aborting"+char(34))
		exit(999)
	}
	root_q = sqrt(cols(scpc_Wfin) - 1)
	SE = (proj_norm / (root_q * sqrt(rows(scpc_Wfin)))) * scale
	if (SE == 0 | missing(SE)) {
		errprintf("SCPC standard error is outside the floating-point range\n")
		exit(999)
	}
	tau = beta / SE
	if (missing(tau)) {
		errprintf("SCPC t-statistic is outside the floating-point range\n")
		exit(999)
	}
	pval=maxrp(scpc_Omsfin,cols(scpc_Wfin)-1,abs(tau)/root_q)
	cv=scpc_cvfin
	if(scpc_condflag) {
		pval=max((pval,maxrp(scpc_Omsxfin,cols(scpc_Wx)-1,abs(tau)/root_q)))
		cv=max((cv,scpc_cvxfin))
	}
	st_matrix(outmat, (beta, SE, tau, pval, beta-cv*SE,beta+cv*SE ))	// coef, SE, tstat, p-value, 95% CI
}

void set_ccspc()
{
	external struct scpc_mats vector scpc_Omsxfin
	external real matrix scpc_geometry
	external struct scpc_pair_sample scalar scpc_cond_pairs
	external real scalar scpc_c0,scpc_cmax,scpc_cvxfin,scpc_approxflag
	external real matrix scpc_distmat, scpc_Wx
	real scalar is_approx, scale
	real matrix weights
	
	is_approx = (length(scpc_approxflag) > 0 ? (scpc_approxflag != 0) : (rows(scpc_Wx) >= 4500))
	// One common scale cancels from all rejection probabilities for this coefficient.
	scale = max(abs(scpc_Wx))
	weights = (scale > 0 ? scpc_Wx / scale : scpc_Wx)
	if(is_approx==0){
		scpc_Omsxfin=getOms(scpc_distmat,scpc_c0,scpc_cmax,weights)
		}
	else{
		// All coefficients use the first conditional pair chain after the unconditional grid.
		if (rows(scpc_cond_pairs.inds) == 0) scpc_cond_pairs = scpc_draw_pairs(scpc_geometry)
		scpc_Omsxfin=LNgetOms(scpc_c0,scpc_cmax,weights,scpc_cond_pairs)
	}
	scpc_cvxfin=getcv(scpc_Omsxfin,cols(scpc_Wx)-1,0.05)
}

void scpc_iv_rank_error()
{
	errprintf("IV residualizer is rank deficient\n")
	exit(999)
}

real matrix scpc_design_basis(real matrix X, real scalar full_rank)
{
    real matrix A, U, Vt
    real rowvector scale
    real colvector s
    real scalar rank

    if (cols(X) == 0) return(J(rows(X), 0, 0))
    scale = colmax(abs(X))
    A = X :/ (scale + (scale :== 0))
    if (rows(A) >= cols(A)) svd(A, U = ., s = ., Vt = .)
    else {
        svd(A', Vt = ., s = ., U = .)
        U = U'
    }
    rank = sum(s :> 1e-10 * s[1])
    if (full_rank & rank < cols(X)) {
        errprintf("OLS residualizer is rank deficient\n")
        exit(999)
    }
    if (rank == 0) return(J(rows(X), 0, 0))
    return(U[., 1..rank])
}

real matrix scpc_iv_projector(real matrix X, real matrix Z)
{
	real rowvector scale
	real matrix Qz, U, Vt, V
	real colvector s
	real scalar p

	p = cols(X)
	if (p == 0 | cols(Z) < p | rows(X) < p | rows(X) != rows(Z) |
		hasmissing(X) | hasmissing(Z)) scpc_iv_rank_error()

	// Normalize X in place so callers use the matching right factor.
	scale = colmax(abs(X))
	if (any(scale :== 0)) scpc_iv_rank_error()
	X = X :/ scale

	Qz = scpc_design_basis(Z, 0)
	if (cols(Qz) < p) scpc_iv_rank_error()

	// X is already equilibrated. A relative test on Qz'X alone can
	// accept a whole projected design made up of roundoff directions.
	svd(cross(Qz, X), U = ., s = ., Vt = .)
	if (hasmissing(s) | min(s) <= 1e-10 * max((1, max(s)))) scpc_iv_rank_error()
	V = Qz * ((U :/ s') * Vt)
	if (hasmissing(V)) scpc_iv_rank_error()
	return(V)
}

// Compact 2SLS coefficient matrices without losing factor/time-series names.
void scpc_select_iv_columns(string scalar b, string scalar bread, string scalar columns)
{
	real rowvector keep
	real matrix coefficients, covariance
	string matrix stripes

	keep = strtoreal(tokens(columns))
	coefficients = st_matrix(b)
	covariance = st_matrix(bread)
	stripes = st_matrixcolstripe(b)
	st_matrix(b, coefficients[., keep])
	st_matrixcolstripe(b, stripes[keep, .])
	st_matrix(bread, covariance[keep, keep])
	st_matrixrowstripe(bread, stripes[keep, .])
	st_matrixcolstripe(bread, stripes[keep, .])
}

// Q Q' for OLS and the transposed 2SLS hat matrix for IV.
void scpc_prepare_projection(string scalar x_vars, string scalar z_vars,
    string scalar sample_sel, real matrix left, real matrix right)
{
    real matrix X, Z
    X = st_data(., tokens(x_vars), sample_sel)
    if (z_vars == "") {
        left = scpc_design_basis(X, 1)
        right = left
    }
    else {
        Z = st_data(., tokens(z_vars), sample_sel)
        left = scpc_iv_projector(X, Z)
        right = X
    }
}

real rowvector scpc_column_norms(real matrix X)
{
    real rowvector scale
    scale = colmax(abs(X))
    return(scale :* sqrt(colsum((X :/ (scale + (scale :== 0))):^2)))
}

void scpc_conditional_weights(real matrix left, real matrix right,
    string scalar sample_sel, string scalar cluster_var,
    string scalar weight_before, string scalar weight_after)
{
    external real scalar scpc_condflag
    external real matrix scpc_Wfin, scpc_Wx
    external real vector scpc_permfin
    real matrix W, R, RX, residual, weighted, panels
    real colvector xb, xa, groups, first, singletons
    real scalar pre_max, q

    scpc_condflag = 1
    q = cols(scpc_Wfin) - 1
    xb = st_data(., weight_before, sample_sel)
    if (max(abs(xb)) > 0 & max(abs(xb)) < .) xb = xb / max(abs(xb))
    xa = st_data(., weight_after, sample_sel)
    W = J(rows(scpc_Wfin), cols(scpc_Wfin), 0)
    W[scpc_permfin, .] = scpc_Wfin
    if (cluster_var != "") {
        groups = st_data(., cluster_var, sample_sel)
        panels = panelsetup(groups, 1)
        singletons = select(panels[., 1], panels[., 1] :== panels[., 2])
        // A zero singleton uses the same sign as an unclustered zero.
        xa[singletons] = xa[singletons] :+ (xa[singletons] :== 0)
        R = W[groups, .]
        first = panelsum(R[., 1] :* xb :* xa, panels)
    }
    else {
        R = W
        first = R[., 1] :* xb :* xa
    }
    scpc_Wx = scpc_Wfin
    scpc_Wx[., 1] = first[scpc_permfin]
    if (q > 0) {
        RX = R[., 2..(q+1)] :* xb
        residual = RX - left * cross(right, RX)
        weighted = residual :* xa
        if (cluster_var != "") {
            weighted = panelsum(weighted, panels)
            pre_max = max(scpc_column_norms(W[., 2..(q+1)] :* (sqrt(rows(W)) * first)))
        }
        else pre_max = max(scpc_column_norms(RX))
        scpc_Wx[., 2..(q+1)] = weighted[scpc_permfin, .]
        if (pre_max == 0 | pre_max >= . | hasmissing(weighted) |
            max(scpc_column_norms(weighted)) <= 1e-12 * pre_max) {
            errprintf("conditional projection collapsed all non-constant W columns\n")
            exit(999)
        }
    }
    set_ccspc()
}

void set_scpccvs(string scalar outmat)
{
	external struct scpc_mats vector scpc_Omsfin, scpc_Omsxfin
	external real matrix scpc_Wfin
	external real scalar scpc_condflag, scpc_cvfin, scpc_cvxfin
	real scalar q
	real rowvector cvs

	q = cols(scpc_Wfin) - 1
	cvs = (getcv(scpc_Omsfin, q, .32), getcv(scpc_Omsfin, q, .10), scpc_cvfin, getcv(scpc_Omsfin, q, .01))
	if (scpc_condflag) {
		cvs = colmax(cvs \ (getcv(scpc_Omsxfin, q, .32), getcv(scpc_Omsxfin, q, .10), scpc_cvxfin, getcv(scpc_Omsxfin, q, .01)))
	}
	st_matrix(outmat, cvs)
}

void get_scpc_meta(string scalar c0name, string scalar cvname, string scalar qname)
{
	external real scalar scpc_c0, scpc_cvfin
	external real matrix scpc_Wfin
	st_numscalar(c0name, scpc_c0)
	st_numscalar(cvname, scpc_cvfin)
	st_numscalar(qname, cols(scpc_Wfin) - 1)
}

end

capture program drop scpc_setscores
program scpc_setscores
	syntax, scpc_sel(name) score_vars(string) bread(name) [IVColumns(string) RAWSCORES]
	tempvar e
	local slist ""
	quietly{
	if(e(cmd)=="regress" | e(cmd)=="logit" | e(cmd)=="probit"){
		predict double `e' if e(sample),sc
		local na : colnames e(V)
		if "`ivcolumns'" == "" {
			local n_b : word count `na'
			forvalues j = 1/`n_b' {
				local ivcolumns "`ivcolumns' `j'"
			}
		}
		local i = 1
		foreach j of local ivcolumns {
			local term : word `j' of `na'
			local score_i : word `i' of `score_vars'
			capture drop `score_i'
			if ("`term'" == "_cons") {
				gen double `score_i' = `e' if e(sample)
			}
			else {
				quietly fvrevar `term' if e(sample)
				gen double `score_i' = `e' * `r(varlist)' if e(sample)
			}
			local slist "`slist' `score_i'"
			local ++i
		}
	}
	else if(e(cmd)=="areg"){
		predict double `e' if e(sample),sc
		local na : colnames e(V)
		if "`ivcolumns'" == "" {
			local n_b : word count `na'
			forvalues j = 1/`n_b' {
				local ivcolumns "`ivcolumns' `j'"
			}
		}
		// Materialize factor/interaction/lag terms in original dataset order before sorting by absvar.
		local fv_vars ""
		foreach j of local ivcolumns {
			local term : word `j' of `na'
			if ("`term'" == "_cons") {
				local fv_vars "`fv_vars' _cons"
			}
			else {
				quietly fvrevar `term' if e(sample)
				local fv_vars "`fv_vars' `r(varlist)'"
			}
		}
		tempvar initial_order
		gen long `initial_order' = _n
		sort `e(absvar)'
		local n_act : word count `ivcolumns'
		forvalues act = 1/`n_act' {
			local fv : word `act' of `fv_vars'
			local score_i : word `act' of `score_vars'
			capture drop `score_i'
			if ("`fv'" == "_cons") {
				gen double `score_i' = `e' if e(sample)
			}
			else {
				tempvar x
				by `e(absvar)': egen double `x' = mean(`fv') if e(sample)
				sum `fv' if e(sample)
				gen double `score_i' = `e' * (`fv' - `x' + `r(mean)') if e(sample)
				drop `x'
			}
			local slist "`slist' `score_i'"
		}
		sort `initial_order'
	}
	else if(e(cmd)=="ivregress"){
		tempvar stata_scr
		predict double `stata_scr'* if e(sample), scores
		predict double `e' if e(sample), residuals

		local instd `"`e(instd)'"'
		local instd : subinstr local instd "bn." ".", all
		local na : colnames e(b)
		if "`ivcolumns'" == "" {
			local n_b : word count `na'
			forvalues j = 1/`n_b' {
				local ivcolumns "`ivcolumns' `j'"
			}
		}
		local i = 1
		foreach j of local ivcolumns {
			local term : word `j' of `na'
			local score_i : word `i' of `score_vars'
			// e(b) can mark an active factor level with bn. while e(instd) uses .
			local key : subinstr local term "bn." ".", all
			local score_pos : list posof "`key'" in instd
			capture drop `score_i'
			if `score_pos' > 0 {
				gen double `score_i' = `stata_scr'`score_pos' if e(sample)
			}
			else if "`term'" == "_cons" {
				gen double `score_i' = `e' if e(sample)
			}
			else {
				quietly fvrevar `term' if e(sample)
				gen double `score_i' = `e' * `r(varlist)' if e(sample)
			}
			local slist "`slist' `score_i'"
			local ++i
		}
	}
	else {
		di "scpc_setscores does not support " e(cmd)
		exit(999)
	}
	
	matrix colnames `bread' = `slist'
	gen byte `scpc_sel'=0
	if( "`e(clustvar)'"=="") qui replace `scpc_sel'=1 if e(sample)
	else{
		tempvar score_sort_sample
		gen byte `score_sort_sample' = !e(sample)
		sort `e(clustvar)' `score_sort_sample'
		if "`rawscores'" == "" scpc_sum_scores `score_vars'
		by `e(clustvar)': replace `scpc_sel' = 1 if _n==1 & e(sample)
	}
	}
end


capture program drop scpc_sum_scores
program scpc_sum_scores
    syntax varlist(numeric)
    tempvar total
    foreach score of local varlist {
        quietly by `e(clustvar)': egen double `total' = total(`score') if e(sample)
        quietly replace `score' = `total' if e(sample)
        drop `total'
    }
end

capture program drop scpc
program scpc, eclass sortpreserve
	version 16.0
	syntax , ///
		[ ///
			cvs	///
			k(real 10) ///
			avc(real -1) ///
			latlong		///
			uncond   ///
			exact ///
			approx ///
			seed(real 1) ///
			backend(string) ///
			nthreads(integer -1) ///
			CAPabilities ///
		]
		//	syntax [if] [in], [
//		avc(real 1.0)
//		]
	
	// The SPUR wrapper uses this option to check ado-level compatibility.
	// Return without reading or replacing the active estimation results.
	if "`capabilities'" != "" exit

	if `"`e(wtype)'"' != "" | `"`e(wexp)'"' != "" {
		di as error "SCPC does not support weighted estimation"
		exit 198
	}

	if "`e(cmd)'" == "ivregress" & "`e(estimator)'" != "2sls" & "`uncond'" == "" {
		di as error "Conditional IV SCPC requires ivregress 2sls"
		exit 198
	}

	if "`exact'" != "" & "`approx'" != "" {
		di as error "cannot specify both exact and approx options"
		exit 198
	}
	if `seed' < 0 | `seed' >= 4294967296 | `seed' != floor(`seed') {
		di as error "seed() must be an integer in [0, 2^32)"
		exit 198
	}
	local seed_arg = trim("`: display %15.0f `seed''")

	local req_backend = lower(trim(`"`backend'"'))
	if "`req_backend'" == "" {
		local req_backend = lower(trim(`"$SPATIAL_BACKEND"'))
	}
	if "`req_backend'" == "" {
		local env_backend : environment SPATIAL_BACKEND
		local req_backend = lower(trim(`"`env_backend'"'))
	}
	if "`req_backend'" == "" {
		local req_backend "auto"
	}
	if "`req_backend'" != "auto" & "`req_backend'" != "rust" & "`req_backend'" != "fallback" {
		di as error "invalid backend() option: `req_backend' (expected auto, rust, or fallback)"
		exit 198
	}
	if `nthreads' != -1 & `nthreads' < 1 {
		di as error "nthreads() must be a positive integer (or -1 for default)"
		exit 198
	}
	if missing(`k') | `k' < 1 | `k' != floor(`k') {
		di as error "k() must be a positive integer"
		exit 198
	}

	tempname b_input bread_input scpc_Bread
	matrix `b_input' = e(b)
	matrix `bread_input' = e(V_modelbased)
	local iv_score_options ""
	local n_coefficients = colsof(`b_input')
	if rowsof(`bread_input') != `n_coefficients' | colsof(`bread_input') != `n_coefficients' {
		di as error "e(b) and e(V_modelbased) dimensions do not match"
		exit 198
	}
	if inlist("`e(cmd)'", "regress", "areg", "logit", "probit", "ivregress") {
		tempname omitted
		quietly _ms_omit_info `b_input'
		matrix `omitted' = r(omit)
		local active_columns ""
		forvalues j = 1/`n_coefficients' {
			if `omitted'[1, `j'] == 0 & !(`b_input'[1, `j'] == 0 & `bread_input'[`j', `j'] == 0) {
				local active_columns "`active_columns' `j'"
			}
		}
		if trim("`active_columns'") == "" {
			di as error "SCPC requires at least one estimable coefficient"
			exit 198
		}
		if `: word count `active_columns'' < `n_coefficients' {
			mata scpc_select_iv_columns("`b_input'", "`bread_input'", "`active_columns'")
		}
		local iv_score_options ivcolumns("`active_columns'")
	}
	matrix `scpc_Bread' = `bread_input'
	if( `scpc_Bread'[1,1]==.) {
		disp as text "e(V_modelbased) missing; aborting"
		exit(999)
	}	
	if `avc'==-1 {
		local avc 0.03	// set default
	}
	if missing(`avc') | `avc' <= 0.001 | `avc' >= 0.99 {
		di as error "avc() must lie strictly between 0.001 and 0.99"
		exit 198
	}

	// Build numerically sorted s_1..s_d coordinate list (avoids unab alphabetical/dataset ordering bugs)
	local slist ""
	local d 0
	local slist_ok 1
	capture unab raw_slist : s_*
	if _rc != 0 {
		local slist_ok 0
	}
	else {
		local d : word count `raw_slist'
		forvalues s_idx = 1/`d' {
			local in_list : list posof "s_`s_idx'" in raw_slist
			if `in_list' == 0 {
				local slist_ok 0
				continue, break
			}
			capture confirm numeric variable s_`s_idx', exact
			if _rc != 0 {
				local slist_ok 0
				continue, break
			}
			local slist "`slist' s_`s_idx'"
		}
	}

	if `slist_ok' == 0 {
		disp as text "s_* variables not continuously numbered starting from 1 (s_1, s_2, etc.)"
		exit 999
	}
	foreach svar of local slist {
		qui count if missing(`svar') & e(sample)
		if r(N) > 0 {
			disp as text "missing value(s) in s_* variable; aborting"
			exit 999
		}
	}
	if "`latlong'" != "" {
		if `d' != 2 {
			di as error "latlong option requires exactly s_1 and s_2"
			exit 198
		}
		qui count if e(sample) & !inrange(s_1, -90, 90)
		if r(N) > 0 {
			di as error "latitude s_1 must be in [-90, 90] in the estimation sample"
			exit 198
		}
		qui count if e(sample) & !inrange(s_2, -180, 180)
		if r(N) > 0 {
			di as error "longitude s_2 must be in [-180, 180] in the estimation sample"
			exit 198
		}
	}

	local has_cluster = ("`e(clustvar)'" != "")
	if `has_cluster' == 1 {
		tempvar cl_tag
		qui egen `cl_tag' = tag(`e(clustvar)') if e(sample)
		qui count if `cl_tag' == 1 & e(sample)
		local n_sites = r(N)
		foreach svar of local slist {
			tempvar unique_coords_tag unique_coords_num
			qui egen `unique_coords_tag' = tag(`e(clustvar)' `svar') if e(sample)
			qui egen `unique_coords_num' = total(`unique_coords_tag') if e(sample), by(`e(clustvar)')
			qui count if e(sample) & `unique_coords_num' != 1
			if r(N) > 0 {
				di as error "coordinates not constant within clusters in the estimation sample"
				exit 999
			}
		}
	}
	else {
		qui count if e(sample)
		local n_sites = r(N)
	}

	local approx_mode = ("`approx'" != "" | ("`exact'" == "" & `n_sites' >= 4500))
	local method_used "exact"
	if `approx_mode' == 1 {
		local method_used "approx"
	}

	// Check Rust plugin availability and exact SCPC eligibility
	local use_rust 0
	local is_iv 0
	local n_endog 0
	local n_z_excl 0
	local z_excl_vars ""
	local backend_reason ""
	if "`req_backend'" != "fallback" {
		local plugin_ok 0
		local plugin_reason "compatible scpc_rust.plugin (v7) could not be loaded. Install the matching plugin"
		capture plugin call _scpc_rust_plugin, ping_v7
		if _rc == 0 {
			local plugin_ok 1
		}
		else {
			capture plugin call _scpc_rust_plugin, ping
			if _rc != 0 {
				capture program _scpc_rust_plugin, plugin using("scpc_rust.plugin")
				if _rc != 0 {
					capture findfile "scpc_rust.plugin"
					if _rc == 0 {
						capture program _scpc_rust_plugin, plugin using("`r(fn)'")
					}
				}
			}
			capture plugin call _scpc_rust_plugin, ping_v7
			if _rc == 0 {
				local plugin_ok 1
			}
			else {
				capture plugin call _scpc_rust_plugin, ping
				if _rc == 0 {
					local plugin_reason "loaded scpc_rust.plugin is incompatible with the required API (v7). Install the matching plugin and restart Stata"
				}
			}
		}

		if `plugin_ok' == 1 {
			capture plugin call _scpc_rust_plugin, require_spec 0.3.8
			if _rc != 0 {
				local plugin_ok 0
				local plugin_reason "loaded scpc_rust.plugin needs methods spec 0.3.8 or newer. Install the matching plugin and restart Stata"
			}
		}

		local rust_eligible 1
		local rust_reason ""
		if `plugin_ok' == 0 {
			local rust_eligible 0
			local rust_reason "`plugin_reason'"
		}
		else if "`e(cmd)'" == "regress" {
			local is_iv 0
		}
		else if "`e(cmd)'" == "areg" {
			if "`uncond'" == "" {
				local rust_eligible 0
				local rust_reason "Rust backend requires uncond option after areg"
			}
		}
		else if "`e(cmd)'" == "ivregress" {
			if "`e(estimator)'" != "2sls" {
				local rust_eligible 0
				local rust_reason "Rust backend only supports ivregress 2sls (got e(estimator)==`e(estimator)')"
			}
			else {
				local is_iv 1
			}
		}
		else {
			local rust_eligible 0
			local rust_reason "Rust backend only supports regress, ivregress 2sls, and uncond areg (got `e(cmd)')"
		}

		if `rust_eligible' == 1 & (`avc' <= 0.001 | `avc' >= 0.99) {
			local rust_eligible 0
			local rust_reason "Rust backend requires avc strictly between 0.001 and 0.99"
		}

		if `rust_eligible' == 1 {
			if `n_sites' < 5 {
				local rust_eligible 0
				local rust_reason "Rust backend requires at least 5 observations/clusters (got `n_sites')"
			}
		}

		if `rust_eligible' == 1 {
			tempname b_full bread_full
			matrix `b_full' = `b_input'
			matrix `bread_full' = `bread_input'
			local p = colsof(`b_full')
			if rowsof(`bread_full') != `p' | colsof(`bread_full') != `p' {
				local rust_eligible 0
				local rust_reason "e(b) and e(V_modelbased) dimensions do not match"
			}
			else {
				local b_names : colnames `b_full'
				local has_cons 0
				local xvars ""
				forvalues j = 1/`p' {
					local colname : word `j' of `b_names'
					if `b_full'[1, `j'] == . | `bread_full'[`j', `j'] == . {
						local rust_eligible 0
						local rust_reason "missing entries in e(b) or e(V_modelbased)"
						continue, break
					}
					if `bread_full'[`j', `j'] <= 0 {
						local rust_eligible 0
						local rust_reason "collinear or dropped regressor '`colname'' detected"
						continue, break
					}
					if "`colname'" == "_cons" {
						if `j' != `p' {
							local rust_eligible 0
							local rust_reason "_cons must be the final column of e(b)"
							continue, break
						}
						local has_cons 1
					}
					else if "`e(cmd)'" == "regress" {
						capture quietly fvrevar `colname' if e(sample)
						if _rc != 0 {
							local rust_eligible 0
							local rust_reason "could not materialize regressor '`colname''"
							continue, break
						}
						local xvars "`xvars' `r(varlist)'"
					}
					else if "`uncond'" == "" {
						if strpos("`colname'", ".") > 0 | strpos("`colname'", "#") > 0 {
							local rust_eligible 0
							local rust_reason "factor, interaction, or time-series term (`colname') is not supported in Rust backend"
							continue, break
						}
						capture confirm numeric variable `colname', exact
						if _rc != 0 {
							local rust_eligible 0
							local rust_reason "regressor '`colname'' is not a numeric variable"
							continue, break
						}
						local xvars "`xvars' `colname'"
					}
				}
			}
		}

		if `rust_eligible' == 1 & `is_iv' == 1 & "`uncond'" == "" {
			local instd_list `"`e(instd)'"'
			local instd_list : list instd_list & xvars
			local instrument_terms `"`e(insts)'"'
			local insts_list ""
			foreach term of local instrument_terms {
				quietly _ms_parse_parts `term'
				if r(omit) == 0 local insts_list "`insts_list' `term'"
			}
			local n_endog : word count `instd_list'
			local nx_count : word count `xvars'
			if `n_endog' < 1 | `n_endog' > `nx_count' {
				local rust_eligible 0
				local rust_reason "ivregress e(instd) does not align with e(b) regressors"
			}
			else {
				forvalues j = 1/`n_endog' {
					local xw : word `j' of `xvars'
					local iw : word `j' of `instd_list'
					if "`xw'" != "`iw'" {
						local rust_eligible 0
						local rust_reason "ivregress e(b) ordering differs from e(instd)"
						continue, break
					}
				}
			}
			if `rust_eligible' == 1 {
				local x_exog_vars ""
				local exog_start = `n_endog' + 1
				forvalues j = `exog_start'/`nx_count' {
					local xw : word `j' of `xvars'
					local x_exog_vars "`x_exog_vars' `xw'"
				}
				local missing_exog : list x_exog_vars - insts_list
				if trim("`missing_exog'") != "" {
					local rust_eligible 0
					local rust_reason "ivregress e(insts) does not contain all exogenous regressors"
				}
			}
			if `rust_eligible' == 1 {
				local z_excl_vars : list insts_list - x_exog_vars
				local z_excl_vars : list z_excl_vars - instd_list
				foreach zname of local z_excl_vars {
					if strpos("`zname'", ".") > 0 | strpos("`zname'", "#") > 0 {
						local rust_eligible 0
						local rust_reason "factor, interaction, or time-series instrument (`zname') is not supported in Rust backend"
						continue, break
					}
					capture confirm numeric variable `zname', exact
					if _rc != 0 {
						local rust_eligible 0
						local rust_reason "instrument '`zname'' is not a numeric variable"
						continue, break
					}
				}
				local n_z_excl : word count `z_excl_vars'
				if `n_z_excl' < `n_endog' {
					local rust_eligible 0
					local rust_reason "insufficient excluded instruments in ivregress"
				}
			}
		}
		else if `is_iv' == 1 {
			local n_endog 0
			local n_z_excl 0
		}

		if `rust_eligible' == 1 {
			local use_rust 1
		}
		else if "`req_backend'" == "rust" {
			di as error "`rust_reason'"
			exit 198
		}
		else {
			local backend_reason "`rust_reason'"
		}
	}
	else {
		local backend_reason "requested fallback"
	}

	// Both paths use the same scores, cluster representatives, and row order.
	local n_scores = colsof(`scpc_Bread')
	local score_vars ""
	forvalues j = 1/`n_scores' {
		tempvar score_`j'
		local score_vars "`score_vars' `score_`j''"
	}
	tempvar scpc_sel
	scpc_setscores, scpc_sel(`scpc_sel') score_vars("`score_vars'") bread(`scpc_Bread') `iv_score_options' rawscores

	local tied_val 0
	if `use_rust' == 1 {
		local k_use = min(`p', int(`k'))
		tempname b_sub c0_scal cv_scal q_scal tied_scal scpctab scpc_cvs
		matrix `b_sub' = `b_full'[1..., 1..`k_use']
		matrix `scpctab' = `b_sub'' * (1, 1, 1, 1, 1, 1)
		local cvs_arg "_"
		if "`cvs'" == "cvs" {
			matrix `scpc_cvs' = `b_sub'' * (1, 1, 1, 1)
			local cvs_arg "`scpc_cvs'"
		}

		local latlong_flag = ("`latlong'" != "")
		local uncond_flag = ("`uncond'" != "")
		local cvs_flag = ("`cvs'" != "")
		local plugin_vars "`slist' `score_vars'"
		if `has_cluster' == 1 {
			tempvar cl_idx
			qui gen double `cl_idx' = sum(`scpc_sel') if e(sample)
			local plugin_vars "`plugin_vars' `cl_idx'"
			local plugin_cond "e(sample)"
		}
		else {
			local plugin_cond "`scpc_sel' == 1"
		}
		if `uncond_flag' == 0 {
			local plugin_vars "`plugin_vars' `xvars'"
			if `n_endog' > 0 & `n_z_excl' > 0 {
				local plugin_vars "`plugin_vars' `z_excl_vars'"
			}
		}

		local avc_arg = trim("`: display %21.16g `avc''")
		if `approx_mode' == 1 {
			plugin call _scpc_rust_plugin `plugin_vars' if `plugin_cond', ///
				model_approx_scpc `d' `p' `k_use' `has_cons' "`avc_arg'" `latlong_flag' `uncond_flag' `cvs_flag' "`nthreads'" ///
				"`b_full'" "`bread_full'" "`scpctab'" "`cvs_arg'" "`c0_scal'" "`cv_scal'" "`q_scal'" "`tied_scal'" `has_cluster' `n_endog' `n_z_excl' "`seed_arg'"
		}
		else {
			plugin call _scpc_rust_plugin `plugin_vars' if `plugin_cond', ///
				model_exact_scpc `d' `p' `k_use' `has_cons' "`avc_arg'" `latlong_flag' `uncond_flag' `cvs_flag' "`nthreads'" ///
				"`b_full'" "`bread_full'" "`scpctab'" "`cvs_arg'" "`c0_scal'" "`cv_scal'" "`q_scal'" "`tied_scal'" `has_cluster' `n_endog' `n_z_excl'
		}

		if "`req_backend'" == "auto" & scalar(`tied_scal') == 1 {
			local use_rust 0
			local tied_val 1
			local backend_reason "tied eigenspace"
			capture matrix drop `b_sub' `scpctab' `scpc_cvs'
			capture scalar drop `c0_scal' `cv_scal' `q_scal' `tied_scal'
		}
		else if "`req_backend'" == "rust" & `approx_mode' == 1 & scalar(`tied_scal') == 1 {
			di as error "tied eigenspace in Nystrom approximation"
			exit 198
		}
		else {
			local k `k_use'
			di as text "found " %6.0f `n_sites' " observations / clusters and " %3.0f `d' "-dimensional locations in s_*"
			if `latlong_flag' == 0 {
				di as text "using Euclidean norm to compute distance between locations stored in s_*"
			}
			else {
				di as text "Computing distances on surface of sphere treating s_1 as latitude and s_2 as longitude"
			}
			di as text "SCPC optimal q = " %3.0f scalar(`q_scal') " for maximal average pairwise correlation = " %5.3f `avc'

			local c0_val = scalar(`c0_scal')
			local cv_val = scalar(`cv_scal')
			local q_val = scalar(`q_scal')
			local tied_val = scalar(`tied_scal')
			local backend_used "rust"
		}
	}
	if `use_rust' == 0 {
		if `has_cluster' scpc_sum_scores `score_vars'
		mata set_method_opts(`approx_mode', `seed')
		mata setOmsWfin(`avc',"`scpc_sel'","`latlong'")
		
		qui sum `scpc_sel', detail
		tempname neff b_mat scpc_Bread0 scpctab scpc_cvs coef_mat stats_row cvs_row oreg_hold
		tempvar x_var xs_var
		matrix `neff'=r(sum)
		matrix `b_mat'=`b_input'
		matrix `scpc_Bread0'=`bread_input'
		local clist 
		local weight_before `xs_var'
		local weight_after `x_var'
		local na : colnames `b_mat'
		if(e(cmd)=="regress" & "`uncond'"==""){
			// Section 5.1 uses Xi_s' M_V Xi R for the OLS denominator weights.
			local weight_before `x_var'
			local weight_after `xs_var'
			tempvar ess one_reg cl_idx_fb cluster_order_reg
			quietly generate byte `ess' = e(sample)
			quietly generate double `one_reg' = 1 if `ess'
			local timevar : char _dta[_TStvar]
			if `has_cluster' & "`timevar'" != "" {
				local panelvar : char _dta[_TSpanel]
				quietly generate long `cluster_order_reg' = _n
				sort `panelvar' `timevar'
			}
			tempname omitted_ols
			quietly _ms_omit_info `b_mat'
			matrix `omitted_ols' = r(omit)
			local iv_u_vars ""
			local iv_z_vars ""
			local j = 0
			foreach term of local na {
				local ++j
				if "`term'" == "_cons" local rv "`one_reg'"
				else {
					quietly fvrevar `term' if `ess'
					local rv "`r(varlist)'"
				}
				local clist "`clist' `rv'"
				if `omitted_ols'[1,`j'] == 0 local iv_u_vars "`iv_u_vars' `rv'"
			}
			matrix colnames `scpc_Bread0' = `clist'
			if `has_cluster' {
				if "`timevar'" != "" sort `e(clustvar)' `cluster_order_reg'
				quietly gen long `cl_idx_fb' = sum(`scpc_sel') if `ess'
			}
		}
		else if(e(cmd)=="regress"){
		}	
		else if(e(cmd)=="ivregress" & "`uncond'"==""){
			// Section 5.3 uses Xi_s' M' Xi R where M' = I - P_Z X (X' P_Z X)^-1 X'.
			local weight_before `x_var'
			local weight_after `xs_var'
			local elist `e(insts)'
			local instd `e(instd)'
			local instd : subinstr local instd "bn." ".", all

			tempvar ess one_iv cl_idx_fb cluster_order
			qui gen byte `ess' = e(sample)
			qui gen double `one_iv' = 1 if `ess'
			local timevar : char _dta[_TStvar]
			if `has_cluster' & "`timevar'" != "" {
				// fvrevar needs time order to materialize lags; preserve the cluster order.
				local panelvar : char _dta[_TSpanel]
				quietly generate long `cluster_order' = _n
				sort `panelvar' `timevar'
			}
			quietly fvrevar `elist' if `ess'
			local iv_z_vars "`r(varlist)'"
			local has_cons_iv : list posof "_cons" in na
			if `has_cons_iv' > 0 {
				local iv_z_vars "`iv_z_vars' `one_iv'"
			}

			_estimates hold `oreg_hold', restore
			local clist
			local iv_u_vars
			foreach term of local na {
				if "`term'" == "_cons" {
					local rv "`one_iv'"
				}
				else {
					quietly fvrevar `term' if `ess'
					local rv "`r(varlist)'"
				}
				local iv_u_vars "`iv_u_vars' `rv'"
				local key : subinstr local term "bn." ".", all
				local is_endog : list posof "`key'" in instd
				if `is_endog' > 0 {
					tempvar fitted
					qui regress `rv' `iv_z_vars' if `ess', noconstant
					qui predict double `fitted' if `ess', xb
					local clist "`clist' `fitted'"
				}
				else {
					local clist "`clist' `rv'"
				}
			}
			matrix colnames `scpc_Bread0' = `clist'
			_estimates unhold `oreg_hold'
			if `has_cluster' {
				if "`timevar'" != "" sort `e(clustvar)' `cluster_order'
				qui gen double `cl_idx_fb' = sum(`scpc_sel') if `ess'
			}
		}
		else if("`uncond'"==""){
			disp as text "conditional critical values only supported for ivregress 2sls and regress, use option uncond for only unconditionally valid critical values; aborting"
			exit(999)
		}
		
		if "`uncond'" == "" {
			tempname proj_left proj_right
			mata scpc_prepare_projection("`iv_u_vars'", "`iv_z_vars'", "`ess'", `proj_left'=., `proj_right'=.)
		}

		local k=min(colsof(`b_mat'),int(`k'))
		matrix `b_mat'=`b_mat'[1..., 1..`k']
		local na : colnames `b_mat'
		matrix `scpctab' =`b_mat''*(1,1,1,1,1,1)
		if("`cvs'"=="cvs"){
			matrix `scpc_cvs' =`b_mat''*(1,1,1,1)
		}

		tokenize `na'
		
		// Loop through covariates
		local i = 1	
		while "``i''" != "" {
			tempvar w
			matrix `coef_mat' = `neff'*`scpc_Bread'[`i', 1...]
			qui matrix score double `w' = `coef_mat' if `scpc_sel' == 1		
			if("`uncond'"==""){	
				capture drop `x_var' 
				capture drop `xs_var'			
				matrix `coef_mat' = `neff'*`scpc_Bread0'[`i', 1...]	
				qui matrix score double `x_var' = `coef_mat' if e(sample)						
				local group_id ""
				if !`has_cluster' {
					quietly gen double `xs_var' = cond(missing(`x_var'), ., cond(`x_var' < 0, -1, 1)) if e(sample)
				}
				else {
					local group_id "`cl_idx_fb'"
					tempvar abs_score score_scale scaled_score
					quietly generate double `abs_score' = abs(`x_var') if e(sample)
					quietly by `e(clustvar)': egen double `score_scale' = max(`abs_score') if e(sample)
					quietly generate double `scaled_score' = cond(`score_scale' == 0, 0, `x_var' / `score_scale') if e(sample)
					quietly by `e(clustvar)': egen double `xs_var' = total(`scaled_score'^2) if e(sample)
					quietly replace `xs_var' = `scaled_score' / sqrt(`xs_var') if e(sample) & `xs_var' > 0
					drop `abs_score' `score_scale' `scaled_score'
				}
				mata scpc_conditional_weights(`proj_left', `proj_right', "`ess'", "`group_id'", "`weight_before'", "`weight_after'")
			}
			mata set_scpcstats("`w'","`scpc_sel'","`stats_row'",st_matrix("`b_mat'")[1,`i'])
			matrix `scpctab'[`i', 1] = `stats_row'[1, 1..6]
			if("`cvs'"=="cvs"){
				mata set_scpccvs("`cvs_row'")
				matrix `scpc_cvs'[`i', 1] = `cvs_row'[1, 1..4]
			}
			local ++i
		}

		if "`uncond'" == "" {
			mata: mata drop `proj_left' `proj_right'
			mata: scpc_cond_pairs = scpc_pair_sample()
		}
		tempname c0_scal cv_scal q_scal
		mata get_scpc_meta("`c0_scal'", "`cv_scal'", "`q_scal'")
		local c0_val = scalar(`c0_scal')
		local cv_val = scalar(`cv_scal')
		local q_val = scalar(`q_scal')
		local backend_used "fallback"
	}

	matrix colnames `scpctab' = "Coef" "Std_Err" "  t  " "P>|t|" "95% Conf" "Interval"
	local n = rowsof(`scpctab')
	local ands = `n'*"&"
	local rs &-`ands'
	local pavc : di %5.3f `avc'
	local tit "SCPC Inference for first `k' coefficients"
	matlist `scpctab', border(all) title(`tit') cspec(o2& %12s | %9.0g o2 & %9.0g o2 &o1 %5.2f o1& o2 %6.3f o1 & o2 %9.0g & o1 %9.0g o2&) rspec(`rs')
	if("`cvs'"=="cvs"){
//		mata set_scpccvs()
//		matrix rownames scpc_cvs ="Two-Sided" "One-Sided"
		matrix colnames `scpc_cvs' ="32%" "10%" "5%" "1%" 
		local n = rowsof(`scpc_cvs')
		local ands = `n'*"&"
		local rs &-`ands'
		matlist `scpc_cvs', border(all) title("Two-Sided Critical values of SCPC t-tests")  cspec(o2& %12s | %6.3f o2 & %6.3f o2 & %6.3f o2 & %6.3f o2 &) rspec(`rs')
		ereturn matrix scpccvs = `scpc_cvs'
	}
	else {
		mata: st_matrix("e(scpccvs)", J(0, 0, .))
	}
 	// Return results
	ereturn matrix scpcstats = `scpctab'
	ereturn scalar q = `q_val'
	ereturn scalar c0 = `c0_val'
	ereturn scalar cv = `cv_val'
	ereturn scalar tied = `tied_val'
	ereturn scalar seed = `seed'
	ereturn local method "`method_used'"
	ereturn local backend "`backend_used'"
	ereturn local backend_reason "`backend_reason'"
end


capture program drop analyticV
program analyticV, sortpreserve
// compares stata covariance matrix with covariance matrix implied by scores and sandwich formula
// should return values close to one or something is wrong
    tempname scpc_Bread
    matrix `scpc_Bread'=e(V_modelbased)
    if( `scpc_Bread'[1,1]==.) {
        disp as text "e(V_modelbased) missing; aborting"
        exit(999)
    } 
	if( "`e(clustvar)'"!="") {
		tempvar es
		gen `es'=0 if e(sample)
		sort `e(clustvar)' `es'
	}

    local n_scores = colsof(`scpc_Bread')
    local score_vars ""
    forvalues j = 1/`n_scores' {
        tempvar score_`j'
        local score_vars "`score_vars' `score_`j''"
    }
    tempvar scpc_sel
    scpc_setscores, scpc_sel(`scpc_sel') score_vars("`score_vars'") bread(`scpc_Bread')
    qui corr `score_vars' if `scpc_sel' == 1, cov
    tempname V
//matlist `scpc_Bread'
//mat `V'= r(C)
//matlist `V'
    mat `V'= r(N)*`scpc_Bread'*r(C)*`scpc_Bread''
	di "number of terms in score:", r(N)
	mata Vx=st_matrix("`V'")
	mata V=st_matrix("e(V)")
	di "ratio of STATA covariance matrix and standard Sandwich formula for first 10 elements"
	mata Vx=Vx:/V
	mata Vx[1::min((10,cols(Vx))),1::min((10,cols(Vx)))]
end
