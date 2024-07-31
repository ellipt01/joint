#include <iostream>
#include <cmath>
#include <cfloat>

#include "mgcal.h"
#include "ADMM.h"

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
void	dscal_ (size_t *n, double *scale, double *x, size_t *inc);
void	daxpy_ (size_t *n, double *alpha, double *x, size_t *incx, double *y, size_t *incy);

int		dpotrf_ (char *uplo, size_t *n, double *a, size_t *lda, long *nfo);
int		dpotri_ (char *uplo, size_t *n, double *a, size_t *lda, long *nfo);
int		dgetrf_ (size_t *m, size_t *n, double *a, size_t *lda, long *ipiv, long *info);
int		dgetri_ (size_t *n, double *a, size_t *m, long *ipiv, double *work, size_t *lwork, long *info);

void	dgemv_ (char *trans, size_t *m, size_t *n, double *alpha, double *A, size_t *ldA,
				double *x, size_t *incx, double *beta , double *y, size_t *incy);
void	dsymv_ (char *uplo, size_t *m, double *alpha, double *A, size_t *ldA,
				double *x, size_t *incx, double *beta, double *y, size_t *incy);
void	dgemm_ (char *transa, char *transb, size_t *m, size_t *n, size_t *k, double *alpha, double *A, size_t *ldA,
				double *B, size_t *ldB, double *beta , double *C, size_t *ldC);


/*** public methods ***/
ADMM::ADMM (double lambda1, double lambda2, double mu)
{
	__init__ ();
	_mu_ = mu;
	set_params (lambda1, lambda2);
}

// set lambda1 and lambda2
void
ADMM::set_params (double lambda1, double lambda2)
{
	_lambda1_ = lambda1;
	_lambda2_ = lambda2;
}

// set simultaneous equations to be solved
void
ADMM::simeq (size_t size1, size_t size2, double *f, double *X, bool normalize, double nu, double *lower)
{
	_size1_ = size1;
	_size2_ = size2;

	_f_ = f;
	_X_ = X;

	if (normalize) {
		if (_w_) delete [] _w_;
		_w_ = _normalize_ (_size1_, _size2_, _X_);
	}

	_nu_ = nu;
	if (_nu_ > DBL_EPSILON) {
		_lower_ = new double [_size2_];
		// copy lower bounds and apply weight
		for (size_t j = 0; j < _size2_; j++) _lower_[j] = lower[j] * _w_[j];
		_apply_lower_bound_ = true;
	}
}

// return zeta
double *
ADMM::get_zeta ()
{
	if (!_zeta_) return NULL;

	double	*zeta = new double [_size2_];
	for (size_t j = 0; j < _size2_; j++) zeta[j] = _zeta_[j];
	if (_w_ != NULL) {
		for (size_t j = 0; j < _size2_; j++) zeta[j] /= _w_[j];
	}
	return zeta;
}

// start ADMM iteration
size_t
ADMM::start (const double tol, const size_t maxiter, bool verbos)
{
	// initialize variables
	_initialize_ ();

	// tol: tolerance, maxiter: maximum number of ADMM iterations
	size_t	k = 0;
	while (k < maxiter) {

		_one_cycle_ ();

		_residual_ = _eval_residuals_ ();
		if (verbos && k % 100 == 0)
			fprintf (stderr, "residual[%ld] = %.4e / %.4e\n", k, _residual_, tol);
		if (_residual_ < tol) break;

		k++;
	}
	return k;
}

// recover the input data
double *
ADMM::recover ()
{
	double	*f = new double [_size1_];
	// f = X * zeta
	dgemv_ (&notrans, &_size1_, &_size2_, &done, _X_, &_size1_, _zeta_, &ione, &dzero, f, &ione);
	return f;
}

/*** protected methods ***/
// initialize zeta, s, and u by padding 0,
// also t and v are initialized when bound constraint is applied
void
ADMM::_initialize_ ()
{
	if (_size1_ == 0 || _size2_ == 0) throw std::runtime_error ("size not specified");

	_residual_ = 0.;

	// allocate and initialize
	// model vector
	if (_zeta_) delete [] _zeta_; 
	_zeta_ = new double [_size2_];
	for (size_t i = 0; i < _size2_; i++) _zeta_[i] = 0.;

	// backup
	if (_zeta_prev_) delete [] _zeta_prev_;
	_zeta_prev_ = new double [_size2_];
	for (size_t i = 0; i < _size2_; i++) _zeta_prev_[i] = 0.;

	// slack vector introduced to separate penalty
	if (_s_) delete [] _s_;
	_s_ = new double [_size2_];
	for (size_t i = 0; i < _size2_; i++) _s_[i] = 0.;

	// Lagrange dual
	if (_u_) delete [] _u_;
	_u_ = new double [_size2_];
	for (size_t i = 0; i < _size2_; i++) _u_[i] = 0.;

	// lower bound onstraint
	if (_apply_lower_bound_) {
		// slack vector for bound constraint
		if (_t_) delete [] _t_;
		_t_ = new double [_size2_];
		for (size_t i = 0; i < _size2_; i++) _t_[i] = 0.;
		// Lagrange dual
		if (_v_) delete [] _v_;
		_v_ = new double [_size2_];
		for (size_t i = 0; i < _size2_; i++) _v_[i] = 0.;
	}

}

// compute Ci and CYi:
// Ci = (X.T * X + (mu + nu) * I)^-1
void
ADMM::_calc_Ci_ ()
{
	if (_Ci_ != NULL) delete [] _Ci_;
	_Ci_ = _Cinv_SMW_ (_mu_ + _nu_, _size1_, _size2_, _X_);
}

// update b and by:
// b = X.T * f + mu * (s + u) + nu * (t + v)
void
ADMM::_update_b_ ()
{
	// b = X.T * f + mu * (s + u)
	if (_b_ == NULL) _b_ = new double [_size2_];
	if (_c_ == NULL) {
		_c_ = new double [_size2_];
		dgemv_ (&trans, &_size1_, &_size2_, &done, _X_, &_size1_, _f_, &ione, &dzero, _c_, &ione);
	}
	for (size_t i = 0; i < _size2_; i++) _b_[i] = _c_[i] + _mu_ * (_s_[i] + _u_[i]);
	if (_apply_lower_bound_) {
		for (size_t i = 0; i < _size2_; i++) _b_[i] += _nu_ * (_t_[i] + _v_[i]);
	}
}

// update zeta:
// zeta = (I - X.T * Ci * X) * b / (mu + nu) 
// Ci = (X * X.T + (mu + nu) * I)^-1
void
ADMM::_update_zeta_ ()
{
	if (_Ci_ == NULL) _calc_Ci_ ();
	if (_zeta_) {
		for (size_t i = 0; i < _size2_; i++) _zeta_prev_[i] = _zeta_[i];
		delete [] _zeta_;
	}
	// zeta = (b - X.T * Ci * X * b) / (mu + nu)
	_zeta_ = _eval_zeta_SMW_ (_mu_ + _nu_, _size1_, _size2_, _X_, _Ci_, _b_);
}

// update s:
// s = C1 * S(q, lambda1 / mu), q = zeta - v
void
ADMM::_update_s_ ()
{
	double	ck = _mu_ / (_mu_ + _lambda2_);

#pragma omp parallel for
	for (size_t j = 0; j < _size2_; j++) {
		double	qj = _zeta_[j] - _u_[j];
		double	cj = _soft_threshold_ (qj, _lambda1_ / _mu_);
		_s_[j] = ck * cj; 
	}
}

// update t: lower boubd constraint
// t = max (lower, zeta - v)
void
ADMM::_update_t_ ()
{
#pragma omp parallel for
	for (size_t j = 0; j < _size2_; j++) {
		double	lj = _lower_[j];
		double	qj = _zeta_[j] - _v_[j];
		_t_[j] = (lj < qj) ? qj : lj;
	}
}

// update u: u = u + mu * (s - zeta)
void
ADMM::_update_u_ ()
{
	for (size_t j = 0; j < _size2_; j++) _u_[j] += _mu_ * (_s_[j] - _zeta_[j]);
}

// update v: v = v + nu * (t - zeta)
void
ADMM::_update_v_ ()
{
	for (size_t j = 0; j < _size2_; j++) _v_[j] += _nu_ * (_t_[j] - _zeta_[j]);
}

// perform ADMM iteration at once
void
ADMM::_one_cycle_ ()
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

/*
	evaluate current residual dr:
	dr = MAX (||s - zeta|| / sqrt (M), mu * ||zeta - zeta_prev|| / sqrt (M))
 */
double
ADMM::_eval_residuals_ ()
{
	double	*dt = new double [_size2_];
	double	*dz = new double [_size2_];

#pragma omp parallel for
	for (size_t j = 0; j < _size2_; j++) {
		// dt = s - zeta
		dt[j] = _s_[j] - _zeta_[j];
		// dz = mu * (zeta - zeta_prev)
		dz[j] = _mu_ * (_zeta_[j] - _zeta_prev_[j]);
	}

	double	dr1 = dnrm2_ (&_size2_, dt, &ione) / sqrt ((double) _size2_);
	delete [] dt;
	double	dr2 = dnrm2_ (&_size2_, dz, &ione) / sqrt ((double) _size2_);
	delete [] dz;

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
double *
ADMM::_normalize_ (size_t m, size_t n, double *K)
{
	double	*w = new double [n];
#pragma omp parallel for
	for (size_t j = 0; j < n; j++) {
		double	wj = dnrm2_ (&m, K + j * m, &ione);
		double	wjinv = 1. / wj;
		dscal_ (&m, &wjinv, K + j * m, &ione);
		w[j] = wj;
	}
	return w;
}

// compute inverse using dpotrf and dpotri (cholesky)
void
ADMM::_cholinv_ (char uplo, size_t n, double *C)
{
	long	info;
	dpotrf_ (&uplo, &n, C, &n, &info);
	dpotri_ (&uplo, &n, C, &n, &info);
}

// compute inverse using dgetrf and dgetri (LU)
void
ADMM::_LUinv_ (size_t m, size_t n, double *C)
{
	long	info;
	size_t	min_mn = (m < n) ? m : n; 
	long	*ipiv = new long [min_mn];
	dgetrf_ (&m, &n, C, &m, ipiv, &info);

	double	w;
	size_t	lwork = -1;
	dgetri_ (&n, C, &m, ipiv, &w, &lwork, &info);
	lwork = (size_t) w;
	double	*work = new double [lwork];
	dgetri_ (&n, C, &m, ipiv, work, &lwork, &info);

	delete [] ipiv;
	delete [] work;
}

// compute inverse (K * K.T + coef * I)^-1 for the data-space inversion
// using following Sherman Morrison Woodbury (SMW) formula:
// (c * I + K.T * K)^-1 = (I - K.T * (c * I + K * K.T)^-1 * K) / c
// This method computes Ci = (c * I + K * K.T)^-1
double *
ADMM::_Cinv_SMW_ (double coef, size_t m, size_t n, double *K)
{
	// Ci = (K * K.T + coef * I)^-1
	double	*Ci = new double [m * m];
	dgemm_ (&notrans, &trans, &m, &m, &n, &done, K, &m, K, &m, &dzero, Ci, &m);
	for (size_t i = 0; i < m; i++) Ci[i + i * m] += coef; // Ci = K * K.T + coef * I
#ifdef USE_LUINV
	_LUinv_ (Ci);
#else
	_cholinv_ ('U', m, Ci);
#endif
	return Ci;
}

// compute zeta using SMW formula
// via compute zeta = (I - K.T * (K * K.T + coef * I)^-1 * K) * b / coef
double *
ADMM::_eval_zeta_SMW_ (double coef, size_t m, size_t n, double *K, double *Ci, double *b)
{
	double *zeta = new double [n];

	// y1 = K * b
	double	*y1 = new double [m];
	dgemv_ (&notrans, &m, &n, &done, K, &m, b, &ione, &dzero, y1, &ione);

	// y2 = Ci * K * b
	double	*y2 = new double [m];
	char	uplo = 'U';
	dsymv_ (&uplo, &m, &done, Ci, &m, y1, &ione, &dzero, y2, &ione);
	delete [] y1;

	// zeta = - K.T * Ci * K * b
	dgemv_ (&trans, &m, &n, &dmone, K, &m, y2, &ione, &dzero, zeta, &ione);
	delete [] y2;

	// zeta = b - K.T * Ci * K * b
	daxpy_ (&n, &done, b, &ione, zeta, &ione);

	// zeta = (b - K.T * Ci * K * b) / coef
	if (fabs (coef - 1.) > DBL_EPSILON) {
		double	scale = 1. / coef;
		dscal_ (&n, &scale, zeta, &ione);
	}
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

#ifdef __cplusplus
	}
#endif // __cplusplus

