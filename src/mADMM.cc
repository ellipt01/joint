#include <iostream>
#include <cmath>
#include <cfloat>

#include "mADMM.h"
#include "util.h"

/*** public methods ***/
mADMM::mADMM (double alpha, double lambda, double mu, double nu, mm_real *lower)
{
	__init__ ();
	_mu_ = mu;

	if (nu > 0.) {
		_nu_ = nu;
		_lower_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, lower->m, lower->n, lower->nnz);
		mm_real_memcpy (_lower_, lower);
		_apply_lower_bound_ = true;
	}

	set_params (alpha, lambda);
}

void
mADMM::set_params (double alpha, double lambda)
{
	_alpha_ = alpha;
	_lambda_ = lambda;
}

// set simultaneous equations to be solved
void
mADMM::simeq (mm_real *f, mm_real *g, mm_real *X, mm_real *Y, bool normalize)
{
	if (f->m != g->m) throw std::runtime_error ("size of g and f does not match.");
	if (X->m != Y->m || X->n != Y->n) throw std::runtime_error ("size of X and Y does not match.");
	if (f->m != X->m) throw std::runtime_error ("size of g and X does not match.");
	if (g->m != Y->m) throw std::runtime_error ("size of f and Y does not match.");

	_size1_ = X->m;
	_size2_ = X->n;

	_f_ = f;
	_g_ = g;

	_X_ = X;
	_Y_ = Y;

	if (normalize) {
		_wx_ = _normalize_ (_X_);
		_wy_ = _normalize_ (_Y_);
	}
	// initialize zeta (beta, rho), s, and u, and set to 0
	_initialize_ ();

	if (_apply_lower_bound_) {
		for (size_t j = 0; j < _size2_; j++) {
			_lower_->data[j] *= _wx_->data[j];
			_lower_->data[j + _size2_] *= _wy_->data[j];
		}
	}
}


// return beta
mm_real *
mADMM::get_beta ()
{
	if (!_beta_) return NULL;

	mm_real	*beta = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _beta_->m, 1, _beta_->m);
	mm_real_memcpy (beta, _beta_);

	if (_wx_ != NULL) {
		for (size_t j = 0; j < beta->m; j++) beta->data[j] /= _wx_->data[j];
	}

	return beta;
}

// return rho
mm_real *
mADMM::get_rho ()
{
	if (!_rho_) return NULL;

	mm_real	*rho = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _rho_->m, 1, _rho_->m);
	mm_real_memcpy (rho, _rho_);
	if (_wy_ != NULL) {
		for (size_t j = 0; j < rho->m; j++) rho->data[j] /= _wy_->data[j];
	}

	return rho;
}

// start ADMM iteration
size_t
mADMM::start (const double tol, const size_t maxiter)
{
	size_t	k;
	for (k = 0; k < maxiter; k++) {

		_one_cycle_ ();

		_residual_ = _eval_residuals_ ();
#ifdef VERBOS
		if (k % 100 == 0) {
			fprintf (stderr, "residual[%ld] = %.4e / %.4e\n", k, _residual_, tol);
			check_mem ("used memory is ");
		}
#endif
		if (_residual_ < tol) break;
	}
	return k;
}

// restart ADMM iteration
size_t
mADMM::restart (const double tol, const size_t maxiter)
{
	_initialize_ ();
	return start (tol, maxiter);
}

// recover the input data
void
mADMM::recover (mm_real *f, mm_real *g)
{
	if (f->m != _size1_ || g->m != _size1_) throw std::runtime_error ("size of mm_real invalid");
	
	mm_real_x_dot_yk (false, 1., _X_, _beta_, 0, 0., f);
	mm_real_x_dot_yk (false, 1., _Y_, _rho_,  0, 0., g);

	return;
}

// read _Ci_ from file
void
mADMM::fread_Cinv (FILE *fp)
{
	_CXi_ = mm_real_fread_binary (fp);
	if (_CXi_ == NULL) throw std::runtime_error ("fread CXinv failed");
	_CYi_ = mm_real_fread_binary (fp);
	if (_CYi_ == NULL) throw std::runtime_error ("fread CYinv failed");
}

// write _Ci_ to file
void
mADMM::fwrite_Cinv (FILE *fp)
{
	if (_CXi_) mm_real_fwrite_binary (fp, _CXi_);
	if (_CYi_) mm_real_fwrite_binary (fp, _CYi_);
}

/*** protected methods ***/
// initialize beta, rho, s, and u
void
mADMM::_initialize_ ()
{
	if (_size1_ == 0 || _size2_ == 0) throw std::runtime_error ("size not specified");

	_residual_ = 0.;

	// allocate and initialize
	if (_beta_) mm_real_free (_beta_);
	_beta_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _size2_, 1, _size2_);
	mm_real_set_all (_beta_, 0.);

	if (_rho_) mm_real_free (_rho_);
	_rho_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _size2_, 1, _size2_);
	mm_real_set_all (_rho_, 0.);

	if (_beta_prev_) mm_real_free (_beta_prev_);
	_beta_prev_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _size2_, 1, _size2_);
	mm_real_set_all (_beta_prev_, 0.);

	if (_rho_prev_) mm_real_free (_rho_prev_);
	_rho_prev_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _size2_, 1, _size2_);
	mm_real_set_all (_rho_prev_, 0.);

	if (_s_) mm_real_free (_s_);
	_s_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, 2 * _size2_, 1, 2 * _size2_);
	mm_real_set_all (_s_, 0.);

	if (_u_) mm_real_free (_u_); 
	_u_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, 2 * _size2_, 1, 2 * _size2_);
	mm_real_set_all (_u_, 0.);

	if (_apply_lower_bound_) {
		if (_t_) mm_real_free (_t_);
		_t_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, 2 * _size2_, 1, 2 * _size2_);
		mm_real_set_all (_t_, 0.);

		if (_v_) mm_real_free (_v_); 
		_v_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, 2 * _size2_, 1, 2 * _size2_);
		mm_real_set_all (_v_, 0.);
	}
}

// compute CXi and CYi:
// CXi = (X.T * X + mu * I)^-1, CYi = (Y.T * Y + mu * I)^-1
void
mADMM::_calc_Ci_ ()
{
	if (_CXi_) mm_real_free (_CXi_);
	_CXi_ = _Cinv_SMW_ (_mu_ + _nu_, _X_);
	if (_CYi_) mm_real_free (_CYi_);
	_CYi_ = _Cinv_SMW_ (_mu_ + _nu_, _Y_);
}

// update bx and by:
// bx = X.T * f + mu * (t + v)[1:M], by = Y.T * g + mu * (t + v)[M:2M]
void
mADMM::_update_bx_ ()
{
	// bx = X.T * f + mu * (s + u)[1:M]
	if (_bx_ == NULL) _bx_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _X_->n, 1, _X_->n);
	if (_cx_ == NULL) {
		_cx_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _X_->n, 1, _X_->n);
		mm_real_x_dot_yk (true, 1., _X_, _f_, 0, 0., _cx_);
	}
	for (size_t i = 0; i < _bx_->m; i++) _bx_->data[i] = _cx_->data[i] + _mu_ * (_s_->data[i] + _u_->data[i]);
	if (_apply_lower_bound_) {
		for (size_t i = 0; i < _bx_->m; i++) _bx_->data[i] += _nu_ * (_t_->data[i] + _v_->data[i]);
	}
}

void
mADMM::_update_by_ ()
{
	// by = Y.T * g + mu * (s + u)[M:2*M]
	if (_by_ == NULL) _by_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _Y_->n, 1, _Y_->n);
	if (_cy_ == NULL) {
		_cy_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _Y_->n, 1, _Y_->n);
		mm_real_x_dot_yk (true, 1., _Y_, _g_, 0, 0., _cy_);
	}
	for (size_t i = 0; i < _by_->m; i++) _by_->data[i] = _cy_->data[i] + _mu_ * (_s_->data[i + _size2_] + _u_->data[i + _size2_]);
	if (_apply_lower_bound_) {
		for (size_t i = 0; i < _by_->m; i++) _by_->data[i] += _nu_ * (_t_->data[i + _size2_] + _v_->data[i + _size2_]);
	}
}

void
mADMM::_update_b_ ()
{
	_update_bx_ ();
	_update_by_ ();
}

// update zeta:
// rho  = (I - X.T * CXi * X) * bx / (mu + nu) 
// beta = (I - Y.T * CYi * Y) * by / (mu + nu)
// CXi = (X * X.T + mu * I)^-1, CYi = (Y * Y.T + (mu + nu) * I)^-1
void
mADMM::_update_zeta_ ()
{
	if (_CXi_ == NULL || _CYi_ == NULL) _calc_Ci_ ();
	_update_b_ ();

	if (_beta_) {
		mm_real_memcpy (_beta_prev_, _beta_);
		mm_real_free (_beta_);
	}
	if (_rho_) {
		mm_real_memcpy (_rho_prev_, _rho_);
		mm_real_free (_rho_);
	}

	// beta = (bx - X.T * CXi * X * bx) / (mu + nu)
	_beta_ = _inv_SMW_ (_mu_ + _nu_, _X_, _CXi_, _bx_);
	// rho  = (by - Y.T * CYi * Y * by) / (mu + nu)
	_rho_  = _inv_SMW_ (_mu_ + _nu_, _Y_, _CYi_, _by_);
}

// update s:
// s[gk] = C1 * (1 + C2 / (mu * || q[gk] ||))_+ * q[gk], q = zeta - v
void
mADMM::_update_s_ ()
{
	double	ck = _mu_ / (_mu_ + (1. - _alpha_) * _lambda_);

#pragma omp parallel for
	for (size_t j = 0; j < _size2_; j++) {

		double	q1 = _beta_->data[j] - _u_->data[j];
		double	q2 = _rho_->data[j] - _u_->data[j + _size2_];
		double	snrm = sqrt (pow (q1, 2.) + pow (q2, 2.));

		//double	cj = 1. - _alpha_ * _lambda_ / (_mu_ * snrm); 
		
		double	l = _alpha_ * _lambda_ / (_mu_ * snrm);
		double	cj = _soft_threshold_ (1., l);
		_s_->data[j] = ck * q1 * cj;
		_s_->data[j + _size2_] = ck * q2 * cj;
	}
}

// update t:
// t = max (lower, zeta - v)
void
mADMM::_update_t_ ()
{
#pragma omp parallel for
	for (size_t j = 0; j < _size2_; j++) {
		double	lbj = _lower_->data[j];
		double	qbj = _beta_->data[j] - _v_->data[j];
		_t_->data[j] = (lbj < qbj) ? qbj : lbj;
		double	lrj = _lower_->data[j + _size2_];
		double	qrj = _rho_->data[j] - _v_->data[j + _size2_];
		_t_->data[j + _size2_] = (lrj < qrj) ? qrj : lrj;
	}
}

// update u: u = u + mu * (s - zeta)
void
mADMM::_update_u_ ()
{
	for (size_t j = 0; j < _size2_; j++) {
		_u_->data[j] += _mu_ * (_s_->data[j] - _beta_->data[j]);
		_u_->data[j + _size2_] += _mu_ * (_s_->data[j + _size2_] - _rho_->data[j]); 
	}
}

// update u: u = u + mu * (s - zeta)
void
mADMM::_update_v_ ()
{
	for (size_t j = 0; j < _size2_; j++) {
		_v_->data[j] += _nu_ * (_t_->data[j] - _beta_->data[j]);
		_v_->data[j + _size2_] += _nu_ * (_t_->data[j + _size2_] - _rho_->data[j]); 
	}
}

// perform ADMM iteration at once
void
mADMM::_one_cycle_ ()
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
mADMM::_eval_residuals_ ()
{
	mm_real	*dt = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, 2 * _size2_, 1, 2 * _size2_);
	mm_real	*dz = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, 2 * _size2_, 1, 2 * _size2_);

#pragma omp parallel for
	for (size_t j = 0; j < _size2_; j++) {
		// dt = s - zeta
		dt->data[j] = _s_->data[j] - _beta_->data[j];
		dt->data[j + _size2_] = _s_->data[j + _size2_] - _rho_->data[j];
		// dz = mu * (zeta - zeta_prev)
		dz->data[j] = _mu_ * (_beta_->data[j] - _beta_prev_->data[j]);
		dz->data[j + _size2_] = _mu_ * (_rho_->data[j] - _rho_prev_->data[j]);
	}

	double	dr1 = mm_real_xj_nrm2 (dt, 0) / sqrt ((double) dt->m);
	mm_real_free (dt);
	double	dr2 = mm_real_xj_nrm2 (dz, 0) / sqrt ((double) dz->m);
	mm_real_free (dz);

	return (dr1 >= dr2) ? dr1 : dr2;
}

// soft threshold (overwrited): S(x, l) = (x - l)+
double
mADMM::_soft_threshold_ (double gamma, double lambda)
{
	double	ci = gamma - lambda;
	return (ci >= 0.) ? ci : 0.;
}

/*** private methods ***/
void
mADMM::__init__ ()
{
	_size1_ = 0;
	_size2_ = 0;

	_f_ = NULL;
	_g_ = NULL;

	_X_ = NULL;
	_Y_ = NULL;

	_wx_ = NULL;
	_wy_ = NULL;

	_mu_ = 0.;
	_nu_ = 0.;
	_apply_lower_bound_ = false;
	_lower_ = NULL;

	_beta_ = NULL;
	_rho_ = NULL;

	_beta_prev_ = NULL;
	_rho_prev_ = NULL;

	_s_ = NULL;
	_u_ = NULL;

	_t_ = NULL;
	_v_ = NULL;

	_cx_ = NULL;
	_cy_ = NULL;

	_bx_ = NULL;
	_by_ = NULL;

	_CXi_ = NULL;
	_CYi_ = NULL;

}

