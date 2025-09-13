#include <iostream>
#include <cmath>
#include <cfloat>

#include "mgcal.h"
#include "ADMM.h"

#ifdef __cplusplus
	extern "C" {
#endif // __cplusplus

static char	trans = 'T';
static char	notrans = 'N';
static size_t	ione = 1;
static double	dzero = 0.;
static double	done = 1.;
static double	dmone = -1.;

double	dnrm2_ (size_t *n, double *x, size_t *inc);
void		dscal_ (size_t *n, double *scale, double *x, size_t *inc);
void		daxpy_ (size_t *n, double *alpha, double *x, size_t *incx, double *y, size_t *incy);

int		dpotrf_ (char *uplo, size_t *n, double *a, size_t *lda, long *nfo);
int		dpotri_ (char *uplo, size_t *n, double *a, size_t *lda, long *nfo);
int		dgetrf_ (size_t *m, size_t *n, double *a, size_t *lda, long *ipiv, long *info);
int		dgetri_ (size_t *n, double *a, size_t *m, long *ipiv, double *work, size_t *lwork, long *info);

void		dgemv_ (char *trans, size_t *m, size_t *n, double *alpha, double *A, size_t *ldA,
			  double *x, size_t *incx, double *beta , double *y, size_t *incy);
void		dsymv_ (char *uplo, size_t *m, double *alpha, double *A, size_t *ldA,
			  double *x, size_t *incx, double *beta, double *y, size_t *incy);
void		dgemm_ (char *transa, char *transb, size_t *m, size_t *n, size_t *k, double *alpha, double *A, size_t *ldA,
			  double *B, size_t *ldB, double *beta , double *C, size_t *ldC);

#ifdef __cplusplus
	}
#endif // __cplusplus


/*** public methods ***/
// Constructor
ADMM::ADMM (double lambda1, double lambda2, double mu)
{
	mu_ = mu;
	set_regularization_parameters (lambda1, lambda2);
}

// Destructor
ADMM::~ADMM ()
{
	delete [] lower_;
	delete [] w_;
	delete [] zeta_;
	delete [] zeta_prev_;
	delete [] s_;
	delete [] u_;
	delete [] t_;
	delete [] v_;
	delete [] c_;
	delete [] b_;
	delete [] Ci_;
}

// Sets the regularization parameters.
void
ADMM::set_regularization_parameters (double lambda1, double lambda2)
{
	lambda1_ = lambda1;
	lambda2_ = lambda2;
}

// Sets up the problem with the given data and parameters.
void
ADMM::setup_problem (size_t size1, size_t size2, double *f, double *X, bool normalize, double nu, double *lower)
{
	size1_ = size1;
	size2_ = size2;

	f_ = f;
	X_ = X;

	if (normalize) {
		if (w_) delete [] w_;
		w_ = normalize_matrix (size1_, size2_, X_);
	}

	nu_ = nu;
	if (nu_ > DBL_EPSILON) {
		lower_ = new double [size2_];
		// Copies lower bounds and applies weights.
		for (size_t j = 0; j < size2_; j++) lower_[j] = lower[j] * w_[j];
		apply_lower_bound_ = true;
	}
}

// Retrieves the model vector (zeta).
double *
ADMM::get_model_vector ()
{
	if (!zeta_) return NULL;

	double	*zeta = new double [size2_];
	for (size_t j = 0; j < size2_; j++) zeta[j] = zeta_[j];
	if (w_ != NULL) {
		for (size_t j = 0; j < size2_; j++) zeta[j] /= w_[j];
	}
	return zeta;
}

// Solves the problem by running ADMM iterations until convergence.
size_t
ADMM::solve (const double tol, const size_t maxiter, bool verbos)
{
	// Initializes the ADMM variables.
	initialize_variables ();

	// tol: tolerance, maxiter: maximum number of ADMM iterations
	size_t	k = 0;
	while (k < maxiter) {

		iterate ();

		residual_ = eval_residuals ();
		if (verbos && k % 100 == 0)
			fprintf (stderr, "residual[%ld] = %.4e / %.4e\n", k, residual_, tol);
		if (residual_ < tol) break;

		k++;
	}
	return k;
}

// Recovers the input data vector `f` from the solution `zeta`.
double *
ADMM::recover_data ()
{
	double	*f = new double [size1_];
	// f = X * zeta
	dgemv_ (&notrans, &size1_, &size2_, &done, X_, &size1_, zeta_, &ione, &dzero, f, &ione);
	return f;
}

/*** protected methods ***/
// Initializes all ADMM variables (zeta, s, u) to zero.
// Also initializes t and v if a lower bound constraint is applied.
void
ADMM::initialize_variables ()
{
	if (size1_ == 0 || size2_ == 0) throw std::runtime_error ("size not specified");

	residual_ = 0.;

	// Allocates and initializes memory for all vectors.
	// Model vector
	if (zeta_) delete [] zeta_;
	zeta_ = new double [size2_];
	for (size_t i = 0; i < size2_; i++) zeta_[i] = 0.;

	// Backup for the previous model vector
	if (zeta_prev_) delete [] zeta_prev_;
	zeta_prev_ = new double [size2_];
	for (size_t i = 0; i < size2_; i++) zeta_prev_[i] = 0.;

	// Slack variable for regularization penalty
	if (s_) delete [] s_;
	s_ = new double [size2_];
	for (size_t i = 0; i < size2_; i++) s_[i] = 0.;

	// Lagrange dual for regularization penalty
	if (u_) delete [] u_;
	u_ = new double [size2_];
	for (size_t i = 0; i < size2_; i++) u_[i] = 0.;

	// Lower bound constraint variables
	if (apply_lower_bound_) {
		// Slack variable for the bound constraint
		if (t_) delete [] t_;
		t_ = new double [size2_];
		for (size_t i = 0; i < size2_; i++) t_[i] = 0.;
		// Lagrange dual for the bound constraint
		if (v_) delete [] v_;
		v_ = new double [size2_];
		for (size_t i = 0; i < size2_; i++) v_[i] = 0.;
	}

}

// Computes the intermediate matrix `Ci` for the zeta update step.
// Ci = (X.T * X + (mu + nu) * I)^-1
void
ADMM::compute_Ci ()
{
	if (Ci_ != NULL) delete [] Ci_;
	Ci_ = compute_Cinv_for_SMW (mu_ + nu_, size1_, size2_, X_);
}

// Updates the intermediate vector `b` for the zeta update step.
// b = X.T * f + mu * (s + u) + nu * (t + v)
void
ADMM::update_b ()
{
	// b = X.T * f + mu * (s + u)
	if (b_ == NULL) b_ = new double [size2_];
	if (c_ == NULL) {
		c_ = new double [size2_];
		dgemv_ (&trans, &size1_, &size2_, &done, X_, &size1_, f_, &ione, &dzero, c_, &ione);
	}
	for (size_t i = 0; i < size2_; i++) b_[i] = c_[i] + mu_ * (s_[i] + u_[i]);
	if (apply_lower_bound_) {
		for (size_t i = 0; i < size2_; i++) b_[i] += nu_ * (t_[i] + v_[i]);
	}
}

// Updates the model vector `zeta`.
// zeta = (I - X.T * Ci * X) * b / (mu + nu)
// Ci = (X * X.T + (mu + nu) * I)^-1
void
ADMM::update_zeta ()
{
	if (Ci_ == NULL) compute_Ci ();
	if (zeta_) {
		for (size_t i = 0; i < size2_; i++) zeta_prev_[i] = zeta_[i];
		delete [] zeta_;
	}
	// zeta = (b - X.T * Ci * X * b) / (mu + nu)
	zeta_ = eval_zeta_using_SMW (mu_ + nu_, size1_, size2_, X_, Ci_, b_);
}

// Updates the slack variable `s`.
// s = C1 * S(q, lambda1 / mu), where q = zeta - u and S is the soft-thresholding operator.
void
ADMM::update_s ()
{
	double	ck = mu_ / (mu_ + lambda2_);

#pragma omp parallel for
	for (size_t j = 0; j < size2_; j++) {
		double	qj = zeta_[j] - u_[j];
		double	cj = soft_threshold (qj, lambda1_ / mu_);
		s_[j] = ck * cj;
	}
}

// Updates the slack variable `t` for the lower bound constraint.
// t = max(lower, zeta - v)
void
ADMM::update_t ()
{
#pragma omp parallel for
	for (size_t j = 0; j < size2_; j++) {
		double	lj = lower_[j];
		double	qj = zeta_[j] - v_[j];
		t_[j] = (lj < qj) ? qj : lj;
	}
}

// Updates the dual variable `u`.
// u = u + mu * (s - zeta)
void
ADMM::update_u ()
{
	for (size_t j = 0; j < size2_; j++) u_[j] += mu_ * (s_[j] - zeta_[j]);
}

// Updates the dual variable `v`.
// v = v + nu * (t - zeta)
void
ADMM::update_v ()
{
	for (size_t j = 0; j < size2_; j++) v_[j] += nu_ * (t_[j] - zeta_[j]);
}

// Performs a single full cycle of ADMM iterations.
void
ADMM::iterate ()
{
	update_b ();
	update_zeta ();
	update_s ();
	update_u ();
	if (apply_lower_bound_) {
		update_t ();
		update_v ();
	}
}

/*
	Evaluates the current residual for convergence check.
	dr = MAX (||s - zeta|| / sqrt (M), mu * ||zeta - zeta_prev|| / sqrt (M))
*/
double
ADMM::eval_residuals ()
{
	double	*dt = new double [size2_];
	double	*dz = new double [size2_];

#pragma omp parallel for
	for (size_t j = 0; j < size2_; j++) {
		// dt = s - zeta
		dt[j] = s_[j] - zeta_[j];
		// dz = mu * (zeta - zeta_prev)
		dz[j] = mu_ * (zeta_[j] - zeta_prev_[j]);
	}

	double	dr1 = dnrm2_ (&size2_, dt, &ione) / sqrt ((double) size2_);
	delete [] dt;
	double	dr2 = dnrm2_ (&size2_, dz, &ione) / sqrt ((double) size2_);
	delete [] dz;

	return (dr1 >= dr2) ? dr1 : dr2;
}

// Applies the soft-thresholding operator: S(x, l) = sign(x) * (|x| - l)+
double
ADMM::soft_threshold (double gamma, double lambda)
{
	double	sign = (gamma >= 0.) ? 1. : -1.;
	double	ci = fabs (gamma) - lambda;
	return sign * ((ci >= 0) ? ci : 0.);
}

// Normalizes the columns of a matrix `K` and returns a vector `w`
// containing the Euclidean norm of each original column.
double *
ADMM::normalize_matrix (size_t m, size_t n, double *K)
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

// Computes the inverse of a matrix `C` using Cholesky decomposition.
// Calls the LAPACK routines dpotrf and dpotri.
void
ADMM::cholinv (char uplo, size_t n, double *C)
{
	long	info;
	dpotrf_ (&uplo, &n, C, &n, &info);
	dpotri_ (&uplo, &n, C, &n, &info);
}

// Computes the inverse of a matrix `C` using LU decomposition.
// Calls the LAPACK routines dgetrf and dgetri.
void
ADMM::LUinv (size_t m, size_t n, double *C)
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

// Computes `inv(K * K.T + coef * I)` for use in the Sherman-Morrison-Woodbury formula.
// The computed matrix is the `C_i` (intermediate inverse) in the formula.
double *
ADMM::compute_Cinv_for_SMW (double coef, size_t m, size_t n, double *K)
{
	// Ci = (K * K.T / coef + I)^-1
	double	*Ci = new double [m * m];
	double	alpha = 1. / coef;
	dgemm_ (&notrans, &trans, &m, &m, &n, &alpha, K, &m, K, &m, &dzero, Ci, &m);
	for (size_t i = 0; i < m; i++) Ci[i + i * m] += 1.; // Ci = K * K.T / coef + I
#ifdef USE_LUINV
	LUinv (Ci);
#else
	cholinv ('U', m, Ci);
#endif
	return Ci;
}

// Computes the `zeta` update step using the Sherman-Morrison-Woodbury formula.
// The calculation is: zeta = [ I - (1 / coef) * K.T * Ci * K ] * b / coef,
// where Ci = (K * K.T / coef + I)^-1.
double *
ADMM::eval_zeta_using_SMW (double coef, size_t m, size_t n, double *K, double *Ci, double *b)
{
	if (coef <= 0.) throw std::runtime_error ("coef must be > 0.");

	double	*zeta = new double [n];

	// y1 = K * b
	double	*y1 = new double [m];
	dgemv_ (&notrans, &m, &n, &done, K, &m, b, &ione, &dzero, y1, &ione);

	// y2 = Ci * K * b
	double	*y2 = new double [m];
	char	uplo = 'U';
	dsymv_ (&uplo, &m, &done, Ci, &m, y1, &ione, &dzero, y2, &ione);
	delete [] y1;

	// zeta = - (1. / coef) * K.T * Ci * K * b
	double	scale = - 1. / coef;
	dgemv_ (&trans, &m, &n, &scale, K, &m, y2, &ione, &dzero, zeta, &ione);
	delete [] y2;

	// zeta = b - (1. / coef) * K.T * Ci * K * b
	daxpy_ (&n, &done, b, &ione, zeta, &ione);

	// zeta = [ b - (1. / coef) * K.T * Ci * K * b] / coef
	if (coef > 1.) {
		scale = 1. / coef;
		dscal_ (&n, &scale, zeta, &ione);
	}
	return zeta;
}
