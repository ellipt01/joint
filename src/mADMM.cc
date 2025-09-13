#include <iostream>
#include <cmath>
#include <cfloat>

#include "mgcal.h"
#include "ADMM.h"
#include "mADMM.h"

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
void		dgemv_ (char *trans, size_t *m, size_t *n, double *alpha, double *A, size_t *ldA,
			  double *x, size_t *incx, double *beta, double *y, size_t *incy);

#ifdef __cplusplus
	}
#endif // __cplusplus

/*** public methods ***/
// Contructor
mADMM::mADMM (double lambda1, double lambda2, double mu)
{
	mu_ = mu;
	set_regularization_parameters (lambda1, lambda2);
}

// Destructor
// Destructor
mADMM::~mADMM ()
{
	delete [] wx_;
	delete [] wy_;
	delete [] beta_;
	delete [] rho_;
	delete [] beta_prev_;
	delete [] rho_prev_;
	delete [] lower_;
	delete [] cx_;
	delete [] cy_;
	delete [] bx_;
	delete [] by_;
	delete [] CXi_;
	delete [] CYi_;
}

// Sets up the optimization problem.
void
mADMM::setup_problem (size_t size1_mag, size_t size1_grv, size_t size2,
				double *f, double *g, double *X, double *Y, bool normalize, double nu, double *lower)
{
	if (size1_mag == 0 || size1_grv == 0)
		throw std::runtime_error ("Problem dimension error: The number of data points cannot be zero.");
	if (size2 == 0)
		throw std::runtime_error ("Model dimension error: The number of grid cells cannot be zero.");

	size1_mag_ = size1_mag;
	size1_grv_ = size1_grv;
	size2_ = size2;

	m_ = size1_mag_ + size1_grv_;
	n_ = 2 * size2_;

	f_ = f;
	g_ = g;

	X_ = X;
	Y_ = Y;

	if (normalize) {
		if (wx_) delete [] wx_;
		wx_ = normalize_matrix (size1_mag_, size2_, X_);
		if (wy_) delete [] wy_;
		wy_ = normalize_matrix (size1_grv_, size2_, Y_);
	}

	if (nu > DBL_EPSILON && lower != NULL) {
		nu_ = nu;
		lower_ = new double [n_];
		// Copy lower bounds and apply normalization weights.
		for (size_t j = 0; j < size2_; j++) {
			lower_[j] = lower[j] * wx_[j];
			lower_[j + size2_] = lower[j + size2_] * wy_[j];
		}
		apply_lower_bound_ = true;
	}
}


// Returns the solution vector beta, un-normalizing if necessary.
double *
mADMM::get_magnetization ()
{
	if (!beta_) return NULL;

	double	*beta = new double [size2_];
	for (size_t j = 0; j < size2_; j++) beta[j] = beta_[j];

	if (wx_ != NULL) {
		for (size_t j = 0; j < size2_; j++) beta[j] /= wx_[j];
	}

	return beta;
}

// Returns the solution vector rho, un-normalizing if necessary.
double *
mADMM::get_density ()
{
	if (!rho_) return NULL;

	double	*rho = new double [size2_];
	for (size_t j = 0; j < size2_; j++) rho[j] = rho_[j];

	if (wy_ != NULL) {
		for (size_t j = 0; j < size2_; j++) rho[j] /= wy_[j];
	}

	return rho;
}

// Solves the problem by running ADMM iterations.
size_t
mADMM::solve (const double tol, const size_t maxiter, bool verbos)
{
	// Initialize ADMM variables zeta (beta, rho), s, and u to zero.
	initialize_variables ();

	// tol: convergence tolerance, maxiter: maximum number of iterations
	size_t	k;
	for (k = 0; k < maxiter; k++) {

		iterate ();

		residual_ = eval_residuals ();
		if (verbos && k % 100 == 0)
			fprintf (stderr, "residual[%ld] = %.4e / %.4e\n", k, residual_, tol);
		if (residual_ < tol) break;
	}
	return k;
}

// Reconstructs the data vectors f and g using the final solution.
void
mADMM::recover_data (double *f, double *g)
{
	dgemv_ (&notrans, &size1_mag_, &size2_, &done, X_, &size1_mag_, beta_, &ione, &dzero, f, &ione);
	dgemv_ (&notrans, &size1_grv_, &size2_, &done, Y_, &size1_grv_, rho_,  &ione, &dzero, g, &ione);
	return;
}

/*** protected methods ***/
// Allocates and initializes ADMM variables.
// If bound constraints are active, t and v are also initialized.
void
mADMM::initialize_variables ()
{
	residual_ = 0.;

	// Allocate and initialize variables.
	// Primary variable for magnetization
	if (beta_) delete [] beta_;
	beta_ = new double [size2_];
	for (size_t i = 0; i < size2_; i++) beta_[i] = 0.;

	// Primary variable for density
	if (rho_) delete [] rho_;
	rho_ = new double [size2_];
	for (size_t i = 0; i < size2_; i++) rho_[i] = 0.;

	// Store previous iteration's beta for residual calculation
	if (beta_prev_) delete [] beta_prev_;
	beta_prev_ = new double [size2_];
	for (size_t i = 0; i < size2_; i++) beta_prev_[i] = 0.;

	// Store previous iteration's rho for residual calculation
	if (rho_prev_) delete [] rho_prev_;
	rho_prev_ = new double [size2_];
	for (size_t i = 0; i < size2_; i++) rho_prev_[i] = 0.;

	// Slack variable for the separable penalty term
	if (s_) delete [] s_;
	s_ = new double [n_];
	for (size_t i = 0; i < n_; i++) s_[i] = 0.;

	// Lagrangian dual variable for the constraint zeta = s
	if (u_) delete [] u_;
	u_ = new double [n_];
	for (size_t i = 0; i < n_; i++) u_[i] = 0.;

	// Variables for the lower bound constraint
	if (apply_lower_bound_) {
		// Slack variable for the bound constraint
		if (t_) delete [] t_;
		t_ = new double [n_];
		for (size_t i = 0; i < n_; i++) t_[i] = 0.;
		// Lagrangian dual variable for the constraint zeta = t
		if (v_) delete [] v_;
		v_ = new double [n_];
		for (size_t i = 0; i < n_; i++) v_[i] = 0.;
	}
}

// Updates the right-hand side vectors bx and by for the zeta-update step.
// bx = X^T * f + mu * (s + u)[1:M] + nu * (t + v)[1:M]
void
mADMM::update_bx ()
{
	// bx = X^T * f + mu * (s + u)[1:M]
	if (bx_ == NULL) bx_ = new double [size2_];
	if (cx_ == NULL) {
		cx_ = new double [size2_];
		dgemv_ (&trans, &size1_mag_, &size2_, &done, X_, &size1_mag_, f_, &ione, &dzero, cx_, &ione);
	}
	for (size_t i = 0; i < size2_; i++) bx_[i] = cx_[i] + mu_ * (s_[i] + u_[i]);
	if (apply_lower_bound_) {
		for (size_t i = 0; i < size2_; i++) bx_[i] += nu_ * (t_[i] + v_[i]);
	}
}

void
mADMM::update_by ()
{
	// by = Y^T * g + mu * (s + u)[M+1:2M]
	if (by_ == NULL) by_ = new double [size2_];
	if (cy_ == NULL) {
		cy_ = new double [size2_];
		dgemv_ (&trans, &size1_grv_, &size2_, &done, Y_, &size1_grv_, g_, &ione, &dzero, cy_, &ione);
	}
	for (size_t i = 0; i < size2_; i++) by_[i] = cy_[i] + mu_ * (s_[i + size2_] + u_[i + size2_]);
	if (apply_lower_bound_) {
		for (size_t i = 0; i < size2_; i++) by_[i] += nu_ * (t_[i + size2_] + v_[i + size2_]);
	}
}

void
mADMM::update_b ()
{
	update_bx ();
	update_by ();
}

// Updates zeta (beta and rho) using the Sherman-Morrison-Woodbury formula.
// beta = (I - X^T * CXi * X) * bx / (mu + nu)
// rho  = (I - Y^T * CYi * Y) * by / (mu + nu)
// where CXi = (X * X^T + (mu + nu) * I)^-1
// and   CYi = (Y * Y^T + (mu + nu) * I)^-1
void
mADMM::update_zeta ()
{
	if (CXi_ == NULL || CYi_ == NULL) compute_Ci ();

	if (beta_) {
		for (size_t i = 0; i < size2_; i++) beta_prev_[i] = beta_[i];
		delete [] beta_;
	}
	if (rho_) {
		for (size_t i = 0; i < size2_; i++) rho_prev_[i] = rho_[i];
		delete [] rho_;
	}
	// beta = (bx - X^T * CXi * X * bx) / (mu + nu)
	beta_ = eval_beta_using_SMW (mu_ + nu_, size1_mag_, size2_, X_, CXi_, bx_);
	// rho  = (by - Y^T * CYi * Y * by) / (mu + nu)
	rho_  = eval_rho_using_SMW (mu_ + nu_, size1_grv_, size2_, Y_, CYi_, by_);
}

// Updates the slack variable s via group soft-thresholding.
// s_k = c * max(0, 1 - lambda1 / (mu * ||q_k||)) * q_k, where q = zeta - u
void
mADMM::update_s ()
{
	double	ck = mu_ / (mu_ + lambda2_);

#pragma omp parallel for
	for (size_t j = 0; j < size2_; j++) {

		double	q1 = beta_[j] - u_[j];
		double	q2 = rho_[j] - u_[j + size2_];
		double	qnrm = sqrt (pow (q1, 2.) + pow (q2, 2.));
		double	l = lambda1_ / (mu_ * qnrm);
		double	cj = soft_threshold (1., l);
		s_[j] = ck * q1 * cj;
		s_[j + size2_] = ck * q2 * cj;
	}
}

// Updates the slack variable t by projecting onto the lower bounds.
// t = max(lower, zeta - v)
void
mADMM::update_t ()
{
#pragma omp parallel for
	for (size_t j = 0; j < size2_; j++) {
		double	lbj = lower_[j];
		double	qbj = beta_[j] - v_[j];
		t_[j] = (lbj < qbj) ? qbj : lbj;
		double	lrj = lower_[j + size2_];
		double	qrj = rho_[j] - v_[j + size2_];
		t_[j + size2_] = (lrj < qrj) ? qrj : lrj;
	}
}

// Updates the dual variable u.
void
mADMM::update_u ()
{
#pragma omp parallel for
	for (size_t j = 0; j < size2_; j++) {
		u_[j] += mu_ * (s_[j] - beta_[j]);
		u_[j + size2_] += mu_ * (s_[j + size2_] - rho_[j]);
	}
}

// Updates the dual variable v.
void
mADMM::update_v ()
{
#pragma omp parallel for
	for (size_t j = 0; j < size2_; j++) {
		v_[j] += nu_ * (t_[j] - beta_[j]);
		v_[j + size2_] += nu_ * (t_[j + size2_] - rho_[j]);
	}
}

// Performs one full ADMM iteration.
void
mADMM::iterate ()
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

double
mADMM::eval_residuals ()
{
	double	*dt = new double [n_];
	double	*dz = new double [n_];

#pragma omp parallel for
	for (size_t j = 0; j < size2_; j++) {
		// Primal residual term: dt = s - zeta
		dt[j] = s_[j] - beta_[j];
		dt[j + size2_] = s_[j + size2_] - rho_[j];
		// Dual residual term: dz = mu * (zeta - zeta_prev)
		dz[j] = mu_ * (beta_[j] - beta_prev_[j]);
		dz[j + size2_] = mu_ * (rho_[j] - rho_prev_[j]);
	}

	double	dr1 = dnrm2_ (&n_, dt, &ione) / sqrt ((double) n_);
	delete [] dt;
	double	dr2 = dnrm2_ (&n_, dz, &ione) / sqrt ((double) n_);
	delete [] dz;

	return (dr1 >= dr2) ? dr1 : dr2;
}

// Computes the inverse matrices required for the SMW formula.
// CXi = (X*X^T + (mu + nu)*I)^-1, CYi = (Y*Y^T + (mu + nu)*I)^-1
void
mADMM::compute_Ci ()
{
	if (CXi_) delete [] CXi_;
	CXi_ = compute_Cinv_for_SMW (mu_ + nu_, size1_mag_, size2_, X_);
	if (CYi_) delete [] CYi_;
	CYi_ = compute_Cinv_for_SMW (mu_ + nu_, size1_grv_, size2_, Y_);
}

// Computes beta using the SMW identity.
// This is a wrapper for the generic eval_zeta_using_SMW function.
double *
mADMM::eval_beta_using_SMW (double coef, size_t m, size_t n, double *X, double *CXi, double *bx)
{
	return eval_zeta_using_SMW (coef, m, n, X, CXi, bx);
}

// Computes rho using the SMW identity.
// This is a wrapper for the generic eval_zeta_using_SMW function.
double *
mADMM::eval_rho_using_SMW (double coef, size_t m, size_t n, double *Y, double *CYi, double *by)
{
	return eval_zeta_using_SMW (coef, m, n, Y, CYi, by);
}

// Soft-thresholding operator: S(x, l) = max(x - l, 0).
double
mADMM::soft_threshold (double gamma, double lambda)
{
	double	ci = gamma - lambda;
	return (ci >= 0.) ? ci : 0.;
}
