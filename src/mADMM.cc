#include <iostream>
#include <cmath>
#include <cfloat>

#include "mgcal.h"
#include "ADMM.h"
#include "mADMM.h"

#ifdef __cplusplus
	extern "C" {
#endif // __cplusplus

static char		trans = 'T';
static char		notrans = 'N';
static size_t	ione = 1;
static double	dzero = 0.;
static double	done = 1.;
static double	dmone = -1.;

double	dnrm2_ (size_t *n, double *x, size_t *inc);
void	dgemv_(char *trans, size_t *m, size_t *n, double *alpha, double *A, size_t *ldA,
			   double *x, size_t *incx, double *beta , double *y, size_t *incy);

/*** public methods ***/
mADMM::mADMM (double lambda1, double lambda2, double mu)
{
	__init__ ();
	_mu_ = mu;
	set_params (lambda1, lambda2);
}

// set simultaneous equations to be solved
void
mADMM::simeq (size_t size1_mag, size_t size1_grv, size_t size2,
			  double *f, double *g, double *X, double *Y, bool normalize, double nu, double *lower)
{
	_size1_mag_ = size1_mag;
	_size1_grv_ = size1_grv;
	_size2_ = size2;

	_m_ = _size1_mag_ + _size1_grv_;
	_n_ = 2 * _size2_;

	_f_ = f;
	_g_ = g;

	_X_ = X;
	_Y_ = Y;

	if (normalize) {
		if (_wx_) delete [] _wx_;
		_wx_ = _normalize_ (_size1_mag_, _size2_, _X_);
		if (_wy_) delete [] _wy_;
		_wy_ = _normalize_ (_size1_grv_, _size2_, _Y_);
	}

	if (nu > DBL_EPSILON && lower != NULL) {
		_nu_ = nu;
		_lower_ = new double [_n_];
		// copy lower bounds and apply weight
		for (size_t j = 0; j < _size2_; j++) {
			_lower_[j] = lower[j] * _wx_[j];
			_lower_[j + _size2_] = lower[j + _size2_] * _wy_[j];
		}
		_apply_lower_bound_ = true;
	}
}


// return beta
double *
mADMM::get_beta ()
{
	if (!_beta_) return NULL;

	double	*beta = new double [_size2_];
	for (size_t j = 0; j < _size2_; j++) beta[j] = _beta_[j];

	if (_wx_ != NULL) {
		for (size_t j = 0; j < _size2_; j++) beta[j] /= _wx_[j];
	}

	return beta;
}

// return rho
double *
mADMM::get_rho ()
{
	if (!_rho_) return NULL;

	double	*rho = new double [_size2_];
	for (size_t j = 0; j < _size2_; j++) rho[j] = _rho_[j];

	if (_wy_ != NULL) {
		for (size_t j = 0; j < _size2_; j++) rho[j] /= _wy_[j];
	}

	return rho;
}

// start ADMM iteration
size_t
mADMM::start (const double tol, const size_t maxiter, bool verbos)
{
	// initialize zeta (beta, rho), s, and u, and set to 0
	_initialize_ ();

	// tol: tolerance, maxiter: maximum number of ADMM iterations
	size_t	k;
	for (k = 0; k < maxiter; k++) {

		_one_cycle_ ();

		_residual_ = _eval_residuals_ ();
		if (verbos && k % 100 == 0)
			fprintf (stderr, "residual[%ld] = %.4e / %.4e\n", k, _residual_, tol);
		if (_residual_ < tol) break;
	}
	return k;
}

// recover the input data
void
mADMM::recover (double *f, double *g)
{
	dgemv_ (&notrans, &_size1_mag_, &_size2_, &done, _X_, &_size1_mag_, _beta_, &ione, &dzero, f, &ione);
	dgemv_ (&notrans, &_size1_grv_, &_size2_, &done, _Y_, &_size1_grv_, _rho_,  &ione, &dzero, g, &ione);
	return;
}

/*** protected methods ***/
// initialize beta, rho, s, and u by padding 0,
// also t and v are initialized when bound constraint is applied
void
mADMM::_initialize_ ()
{
	if (_size1_mag_ == 0 || _size1_grv_ == 0) throw std::runtime_error ("number of data not specified");
	if (_size2_ == 0) throw std::runtime_error ("number of grid cells not specified");

	_residual_ = 0.;

	// allocate and initialize
	// magnetization
	if (_beta_) delete [] _beta_;
	_beta_ = new double [_size2_];
	for (size_t i = 0; i < _size2_; i++) _beta_[i] = 0.;

	// density
	if (_rho_) delete [] _rho_;
	_rho_ = new double [_size2_];
	for (size_t i = 0; i < _size2_; i++) _rho_[i] = 0.;

	// backup
	if (_beta_prev_) delete [] _beta_prev_;
	_beta_prev_ = new double [_size2_];
	for (size_t i = 0; i < _size2_; i++) _beta_prev_[i] = 0.;

	// density
	if (_rho_prev_) delete [] _rho_prev_;
	_rho_prev_ = new double [_size2_];
	for (size_t i = 0; i < _size2_; i++) _rho_prev_[i] = 0.;

	// slack vector introduced to separate penalty
	if (_s_) delete [] _s_;
	_s_ = new double [_n_];
	for (size_t i = 0; i < _n_; i++) _s_[i] = 0.;

	// Lagrange dual
	if (_u_) delete [] _u_;
	_u_ = new double [_n_];
	for (size_t i = 0; i < _n_; i++) _u_[i] = 0.;

	// lower bound onstraint
	if (_apply_lower_bound_) {
		// slack vector for bound constraint
		if (_t_) delete [] _t_;
		_t_ = new double [_n_];
		for (size_t i = 0; i < _n_; i++) _t_[i] = 0.;
		// Lagrange dual
		if (_v_) delete [] _v_;
		_v_ = new double [_n_];
		for (size_t i = 0; i < _n_; i++) _v_[i] = 0.;
	}
}

// update bx and by:
// bx = X.T * f + mu * (t + v)[1:M], by = Y.T * g + mu * (t + v)[M:2M]
void
mADMM::_update_bx_ ()
{
	// bx = X.T * f + mu * (s + u)[1:M]
	if (_bx_ == NULL) _bx_ = new double [_size2_];
	if (_cx_ == NULL) {
		_cx_ = new double [_size2_];
		dgemv_ (&trans, &_size1_mag_, &_size2_, &done, _X_, &_size1_mag_, _f_, &ione, &dzero, _cx_, &ione);
	}
	for (size_t i = 0; i < _size2_; i++) _bx_[i] = _cx_[i] + _mu_ * (_s_[i] + _u_[i]);
	if (_apply_lower_bound_) {
		for (size_t i = 0; i < _size2_; i++) _bx_[i] += _nu_ * (_t_[i] + _v_[i]);
	}
}

void
mADMM::_update_by_ ()
{
	// by = Y.T * g + mu * (s + u)[M:2*M]
	if (_by_ == NULL) _by_ = new double [_size2_];
	if (_cy_ == NULL) {
		_cy_ = new double [_size2_];
		dgemv_ (&trans, &_size1_grv_, &_size2_, &done, _Y_, &_size1_grv_, _g_, &ione, &dzero, _cy_, &ione);
	}
	for (size_t i = 0; i < _size2_; i++) _by_[i] = _cy_[i] + _mu_ * (_s_[i + _size2_] + _u_[i + _size2_]);
	if (_apply_lower_bound_) {
		for (size_t i = 0; i < _size2_; i++) _by_[i] += _nu_ * (_t_[i + _size2_] + _v_[i + _size2_]);
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

	if (_beta_) {
		for (size_t i = 0; i < _size2_; i++) _beta_prev_[i] = _beta_[i];
		delete [] _beta_;
	}
	if (_rho_) {
		for (size_t i = 0; i < _size2_; i++) _rho_prev_[i] = _rho_[i];
		delete [] _rho_;
	}
	// beta = (bx - X.T * CXi * X * bx) / (mu + nu)
	_beta_ = _eval_beta_SMW_ (_mu_ + _nu_, _size1_mag_, _size2_, _X_, _CXi_, _bx_);
	// rho  = (by - Y.T * CYi * Y * by) / (mu + nu)
	_rho_  = _eval_rho_SMW_ (_mu_ + _nu_, _size1_grv_, _size2_, _Y_, _CYi_, _by_);
}

// update s:
// s[gk] = C1 * (1 + lambda1 / (mu * || q[gk] ||))_+ * q[gk], q = zeta - v
void
mADMM::_update_s_ ()
{
	double	ck = _mu_ / (_mu_ + _lambda2_);

#pragma omp parallel for
	for (size_t j = 0; j < _size2_; j++) {

		double	q1 = _beta_[j] - _u_[j];
		double	q2 = _rho_[j] - _u_[j + _size2_];
		double	qnrm = sqrt (pow (q1, 2.) + pow (q2, 2.));
		double	l = _lambda1_ / (_mu_ * qnrm);
		double	cj = _soft_threshold_ (1., l);
		_s_[j] = ck * q1 * cj;
		_s_[j + _size2_] = ck * q2 * cj;
	}
}

// update t:
// t = max (lower, zeta - v)
void
mADMM::_update_t_ ()
{
#pragma omp parallel for
	for (size_t j = 0; j < _size2_; j++) {
		double	lbj = _lower_[j];
		double	qbj = _beta_[j] - _v_[j];
		_t_[j] = (lbj < qbj) ? qbj : lbj;
		double	lrj = _lower_[j + _size2_];
		double	qrj = _rho_[j] - _v_[j + _size2_];
		_t_[j + _size2_] = (lrj < qrj) ? qrj : lrj;
	}
}

// update u: u = u + mu * (s - zeta)
void
mADMM::_update_u_ ()
{
#pragma omp parallel for
	for (size_t j = 0; j < _size2_; j++) {
		_u_[j] += _mu_ * (_s_[j] - _beta_[j]);
		_u_[j + _size2_] += _mu_ * (_s_[j + _size2_] - _rho_[j]); 
	}
}

// update u: u = u + mu * (s - zeta)
void
mADMM::_update_v_ ()
{
#pragma omp parallel for
	for (size_t j = 0; j < _size2_; j++) {
		_v_[j] += _nu_ * (_t_[j] - _beta_[j]);
		_v_[j + _size2_] += _nu_ * (_t_[j + _size2_] - _rho_[j]); 
	}
}

// perform ADMM iteration at once
void
mADMM::_one_cycle_ ()
{
	_update_b_ ();
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
	double	*dt = new double [_n_];
	double	*dz = new double [_n_];

#pragma omp parallel for
	for (size_t j = 0; j < _size2_; j++) {
		// dt = s - zeta
		dt[j] = _s_[j] - _beta_[j];
		dt[j + _size2_] = _s_[j + _size2_] - _rho_[j];
		// dz = mu * (zeta - zeta_prev)
		dz[j] = _mu_ * (_beta_[j] - _beta_prev_[j]);
		dz[j + _size2_] = _mu_ * (_rho_[j] - _rho_prev_[j]);
	}

	double	dr1 = dnrm2_ (&_n_, dt, &ione) / sqrt ((double) _n_);
	delete [] dt;
	double	dr2 = dnrm2_ (&_n_, dz, &ione) / sqrt ((double) _n_);
	delete [] dz;

	return (dr1 >= dr2) ? dr1 : dr2;
}

// compute CXi and CYi:
// CXi = (X.T * X + (mu + nu) * I)^-1, CYi = (Y.T * Y + (mu + nu) * I)^-1
void
mADMM::_calc_Ci_ ()
{
	if (_CXi_) delete [] _CXi_;
	_CXi_ = _Cinv_SMW_ (_mu_ + _nu_, _size1_mag_, _size2_, _X_);
	if (_CYi_) delete [] _CYi_;
	_CYi_ = _Cinv_SMW_ (_mu_ + _nu_, _size1_grv_, _size2_, _Y_);
}

// compute beta = (I - X.T * CXi * X) * bx / coef
// interface of _eval_zeta_SMW_ ()
double *
mADMM::_eval_beta_SMW_ (double coef, size_t m, size_t n, double *X, double *CXi, double *bx)
{
	return _eval_zeta_SMW_ (coef, m, n, X, CXi, bx);
}

// compute rho = (I - Y.T * CYi * Y) * by / coef
// interface of _eval_zeta_SMW_ ()
double *
mADMM::_eval_rho_SMW_ (double coef, size_t m, size_t n, double *Y, double *CYi, double *by)
{
	return _eval_zeta_SMW_ (coef, m, n, Y, CYi, by);
}

// soft threshold (overwrited): S(x, l) = max(x - l, 0)
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

	_size1_mag_ = 0;
	_size1_grv_ = 0;

	_f_ = NULL;
	_g_ = NULL;

	_X_ = NULL;
	_Y_ = NULL;

	_wx_ = NULL;
	_wy_ = NULL;

	_beta_ = NULL;
	_rho_ = NULL;

	_beta_prev_ = NULL;
	_rho_prev_ = NULL;

	_cx_ = NULL;
	_cy_ = NULL;

	_bx_ = NULL;
	_by_ = NULL;

	_CXi_ = NULL;
	_CYi_ = NULL;

}

#ifdef __cplusplus
	}
#endif // __cplusplus

