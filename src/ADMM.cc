#include <iostream>
#include <cmath>
#include <cfloat>

#include "ADMM.h"
#include "util.h"

/*** public methods ***/
ADMM::ADMM (double alpha, double lambda, double mu)
{
	__init__ ();
	_mu_ = mu;
	set_params (alpha, lambda);
}

ADMM::ADMM (double alpha, double lambda, double mu, double nu, mm_real *lower)
{
	__init__ ();
	_mu_ = mu;

	_nu_ = nu;
	_lower_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, lower->m, lower->n, lower->nnz);
	mm_real_memcpy (_lower_, lower);
	_apply_lower_bound_ = true;

	set_params (alpha, lambda);
}

// set alpha and lambda
void
ADMM::set_params (double alpha, double lambda)
{
	_alpha_ = alpha;
	_lambda_ = lambda;
}

// set simultaneous equations to be solved
void
ADMM::simeq (mm_real *f, mm_real *X, bool normalize)
{
	if (f->m != X->m) throw std::runtime_error ("size of g and X does not match.");

	_size1_ = X->m;
	_size2_ = X->n;

	_f_ = f;
	_X_ = X;

	if (normalize) {
		if (_w_) mm_real_free (_w_);
		_w_ = _normalize_ (_X_);
	}
	// initialize primal and dual variables, and set to 0
	_initialize_ ();

	if (_apply_lower_bound_) {
		for (size_t j = 0; j < _size2_; j++) _lower_->data[j] *= _w_->data[j];
	}
}

// return zeta
mm_real *
ADMM::get_zeta ()
{
	if (!_zeta_) return NULL;

	mm_real	*zeta = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _zeta_->m, 1, _zeta_->m);
	mm_real_memcpy (zeta, _zeta_);
	if (_w_ == NULL) {
		for (size_t j = 0; j < zeta->m; j++) zeta->data[j] /= _w_->data[j];
	}

	return zeta;
}

// start ADMM iteration
size_t
ADMM::start (const double tol, const size_t maxiter)
{
	size_t	k;
	for (k = 0; k < maxiter; k++) {

		_one_cycle_ ();

		_residual_ = _eval_residuals_ ();
#ifdef VERBOS
		if (k % 100 == 0) {
			fprintf (stderr, "residual[%ld] = %.4e / %.4e\n", k, _residual_, tol);
			check_mem (NULL);
		}
#endif
		if (_residual_ < tol) break;
	}
	return k;
}

// restart ADMM iteration
size_t
ADMM::restart (const double tol, const size_t maxiter)
{
	_initialize_ ();
	return start (tol, maxiter);
}

// recover the input data
mm_real *
ADMM::recover ()
{
	mm_real	*f = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _size1_, 1, _size1_);
	mm_real_x_dot_yk (false, 1., _X_, _zeta_, 0, 0., f);
	return f;
}

// read _Ci_ from file
void
ADMM::fread_Cinv (FILE *fp)
{
	_Ci_ = mm_real_fread_binary (fp);
	if (_Ci_ == NULL) throw std::runtime_error ("fread CXinv failed");
}

// write _Ci_ to file
void
ADMM::fwrite_Cinv (FILE *fp)
{
	if (_Ci_) mm_real_fwrite_binary (fp, _Ci_);
}

/*** protected methods ***/
// initialize zeta, s, and u
void
ADMM::_initialize_ ()
{
	if (_size1_ == 0 || _size2_ == 0) throw std::runtime_error ("size not specified");

	_residual_ = 0.;

	// allocate and initialize
	if (_zeta_) mm_real_free (_zeta_); 
	_zeta_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _size2_, 1, _size2_);
	mm_real_set_all (_zeta_, 0.);

	if (_zeta_prev_) mm_real_free (_zeta_prev_);
	_zeta_prev_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _size2_, 1, _size2_);
	mm_real_set_all (_zeta_prev_, 0.);

	if (_s_) mm_real_free (_s_);
	_s_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _size2_, 1, _size2_);
	mm_real_set_all (_s_, 0.);

	if (_u_) mm_real_free (_u_);
	_u_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _size2_, 1, _size2_);
	mm_real_set_all (_u_, 0.);

	// lower bound onstraint
	if (_apply_lower_bound_) {
		if (_t_) mm_real_free (_t_);
		_t_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _size2_, 1, _size2_);
		mm_real_set_all (_t_, 0.);

		if (_v_) mm_real_free (_v_);
		_v_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _size2_, 1, _size2_);
		mm_real_set_all (_v_, 0.);
	}

}

// compute Ci and CYi:
// Ci = (X.T * X + (mu + nu) * I)^-1
void
ADMM::_calc_Ci_ ()
{
	if (_Ci_ != NULL) mm_real_free (_Ci_);
	_Ci_ = _Cinv_SMW_ (_mu_ + _nu_, _X_);
}

// update b and by:
// b = X.T * f + mu * (s + u) + nu * (t + v)
void
ADMM::_update_b_ ()
{
	// b = X.T * f + mu * (s + u)
	if (_b_ == NULL) _b_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _X_->n, 1, _X_->n);
	if (_c_ == NULL) {
		_c_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _X_->n, 1, _X_->n);
		mm_real_x_dot_yk (true, 1., _X_, _f_, 0, 0., _c_);
	}
	for (size_t i = 0; i < _b_->m; i++) _b_->data[i] = _c_->data[i] + _mu_ * (_s_->data[i] + _u_->data[i]);
	if (_apply_lower_bound_) {
		for (size_t i = 0; i < _b_->m; i++) _b_->data[i] += _nu_ * (_t_->data[i] + _v_->data[i]);
	}
}

// update zeta:
// zeta = (I - X.T * Ci * X) * b / (mu + nu) 
// Ci = (X * X.T + (mu + nu) * I)^-1
void
ADMM::_update_zeta_ ()
{
	if (_Ci_ == NULL) _calc_Ci_ ();
	_update_b_ ();

	if (_zeta_) {
		mm_real_memcpy (_zeta_prev_, _zeta_);
		mm_real_free (_zeta_);
	}
	// zeta = (b - X.T * Ci * X * b) / (mu + nu)
	_zeta_ = _inv_SMW_ (_mu_ + _nu_, _X_, _Ci_, _b_);
}

// update s:
// s = C1 * S(q, lambda / mu), q = zeta - v
void
ADMM::_update_s_ ()
{
	double	ck = _mu_ / (_mu_ + (1. - _alpha_) * _lambda_);

#pragma omp parallel for
	for (size_t j = 0; j < _size2_; j++) {
		double	qj = _zeta_->data[j] - _u_->data[j];
		double	cj = _soft_threshold_ (qj, _alpha_ * _lambda_ / _mu_);
		_s_->data[j] = ck * cj; 
	}
}

// update t:
// t = max (lower, zeta - v)
void
ADMM::_update_t_ ()
{
#pragma omp parallel for
	for (size_t j = 0; j < _size2_; j++) {
		double	lj = _lower_->data[j];
		double	qj = _zeta_->data[j] - _v_->data[j];
		_t_->data[j] = (lj < qj) ? qj : lj;
	}
}

// update u: u = u + mu * (s - zeta)
void
ADMM::_update_u_ ()
{
	for (size_t j = 0; j < _size2_; j++)
		_u_->data[j] += _mu_ * (_s_->data[j] - _zeta_->data[j]);
}

// update v: v = v + nu * (t - zeta)
void
ADMM::_update_v_ ()
{
	for (size_t j = 0; j < _size2_; j++)
		_v_->data[j] += _nu_ * (_t_->data[j] - _zeta_->data[j]);
}

// perform ADMM iteration at once
void
ADMM::_one_cycle_ ()
{
	_update_zeta_ ();
	_update_s_ ();
	_update_u_ ();
	if (_apply_lower_bound_) {
		_update_t_ ();
		_update_v_ ();
	}
}

double
ADMM::_eval_residuals_ ()
{
	mm_real	*dt = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _size2_, 1, _size2_);
	mm_real	*dz = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _size2_, 1, _size2_);

#pragma omp parallel for
	for (size_t j = 0; j < _size2_; j++) {
		// dt = s - zeta
		dt->data[j] = _s_->data[j] - _zeta_->data[j];
		// dz = mu * (zeta - zeta_prev)
		dz->data[j] = _mu_ * (_zeta_->data[j] - _zeta_prev_->data[j]);
	}

	double	dr1 = mm_real_xj_nrm2 (dt, 0) / sqrt ((double) dt->m);
	mm_real_free (dt);
	double	dr2 = mm_real_xj_nrm2 (dz, 0) / sqrt ((double) dz->m);
	mm_real_free (dz);

	return (dr1 >= dr2) ? dr1 : dr2;
}

// soft threshold: S(x, l) = sign(x) * (|x| - l)+
double
ADMM::_soft_threshold_ (double gamma, double lambda)
{
	double	sign = (gamma >= 0.) ? 1. : -1.;
	double	ci = fabs (gamma) - lambda;
	return sign * ((ci >= 0) ? ci : 0.);
}

// normalize matrix *K and store || kj || in vector *w
mm_real *
ADMM::_normalize_ (mm_real *K)
{
	mm_real	*w = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, K->n, 1, K->n);
#pragma omp parallel for
	for (size_t j = 0; j < K->n; j++) {
		double	wj = mm_real_xj_nrm2 (K, j);
		mm_real_xj_scale (K, j, 1. / wj);
		w->data[j] = wj;
	}
	return w;
}

// compute inverse using dpotrf and dpotri (cholesky)
void
ADMM::_cholinv_ (mm_real *C)
{
	char	uplo = (C->symm == MM_REAL_SYMMETRIC_UPPER) ? 'U' : 'L';
	int		info;
	dpotrf_ (&uplo, (int *) &C->m, C->data, (int *) &C->m, &info);
	dpotri_ (&uplo, (int *) &C->m, C->data, (int *) &C->m, &info);
}

// compute inverse using dgetrf and dgetri (LU)
void
ADMM::_LUinv_ (mm_real *C)
{
	int		info;
	size_t	m = C->m;
	size_t	n = C->n;
	size_t	min_mn = (m < n) ? m : n; 
	int	*ipiv = new int [min_mn];
	dgetrf_ ((int *) &m, (int *) &n, C->data, (int *) &m, ipiv, &info);

	double	w;
	int		lwork = -1;
	dgetri_ ((int *) &n, C->data, (int *) &m, ipiv, &w, &lwork, &info);
	lwork = (int) w;
	double	*work = new double [lwork];
	dgetri_ ((int *) &n, C->data, (int *) &m, ipiv, work, &lwork, &info);

	delete [] ipiv;
	delete [] work;
}

// compute (K * K.T + mu * I)^-1:
mm_real *
ADMM::_Cinv_SMW_ (double mu, mm_real *K)
{
	size_t	m = K->m;
	size_t	n = K->n;
	// Ci = (K * K.T + mu * I)^-1
	mm_real	*Ci = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, m, m * m);
	mm_real_x_dot_y (false, true, 1., K, K, 0., Ci); // Ci = K * K.T
	for (size_t i = 0; i < Ci->m; i++) Ci->data[i + i * Ci->m] += mu; // Ci = K * K.T + mu * I

#ifdef USE_LUINV
	_LUinv_ (Ci);
#else
	mm_real_general_to_symmetric ('U', Ci);
	_cholinv_ (Ci);
#endif
	return Ci;
}

// compute (I - K.T * (K * K.T + mu * I)^-1 * K) * b / mu
mm_real *
ADMM::_inv_SMW_ (double mu, mm_real *K, mm_real *Ci, mm_real *b)
{
	size_t	m = K->m;
	size_t	n = K->n;

	mm_real	*zeta = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, n, 1, n);

	// y1 = K * b
	mm_real	*y1 = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, 1, m);
	mm_real_x_dot_yk (false, 1., K, b, 0, 0., y1);

	// y2 = Ci * K * b
	mm_real	*y2 = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, 1, m);
	mm_real_x_dot_yk (false, 1., Ci, y1, 0, 0., y2);
	mm_real_free (y1);

	// zeta = - K.T * Ci * K * b
	mm_real_x_dot_yk (true, -1., K, y2, 0, 0., zeta);
	mm_real_free (y2);

	// zeta = b - K.T * Ci * K * b
	mm_real_axjpy (1., b, 0, zeta);

	// zeta = (b - K.T * Ci * K * b) / mu
	if (fabs (_mu_ - 1.) > DBL_EPSILON) mm_real_xj_scale (zeta, 0, 1. / mu);

	return zeta;
}

/*** private methods ***/
void
ADMM::__init__ ()
{
	_size1_ = 0;
	_size2_ = 0;

	_f_ = NULL;
	_X_ = NULL;
	_w_ = NULL;

	_mu_ = 0.;
	_nu_ = 0.;
	_apply_lower_bound_ = false;
	_lower_ = NULL;

	_zeta_ = NULL;
	_zeta_prev_ = NULL;

	_s_ = NULL;
	_u_ = NULL;

	_t_ = NULL;
	_v_ = NULL;

	_c_ = NULL;
	_b_ = NULL;
	_Ci_ = NULL;
}


