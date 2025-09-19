#ifndef _MADMM_H_
#define _MADMM_H_

enum class FieldType {
	Magnetic,
	Gravity
};

enum class ModelType {
	Magnetic,
	Gravity
};

/*
 * mADMM: A subclass of ADMM that extends the algorithm for joint magnetic and gravity data inversion.
 *
 * This class is designed to solve a linear problem with an L2-norm and Group Lasso combined regularization:
 *
 * min (1/2) ||b - Z * zeta||^2 + lambda1 * sum_j=0^M ||zeta_gj|| + lambda2 ||zeta||^2/2,
 *
 * where b = [f;g] (magnetic and gravity data), and Z = [X;Y] (magnetic and gravity kernel matrices).
 * zeta = [beta; rho] where beta is the magnetization and rho is the density model vector.
 * The second term is a Group Lasso penalty, with M being the number of grid cells.
 *
 * The mADMM algorithm solves the equivalent constrained optimization problem:
 *
 * min (1/2) ||b - Z * zeta||^2 + lambda1 * sum_j=0^M ||s_gj|| + lambda2 ||s||^2/2  s.t. s = zeta
 *
 * The augmented Lagrangian function for this problem is:
 *
 * L = (1/2) ||f - X * zeta||^2 + lambda1 * |s| + lambda2 ||s||^2/2 + (mu / 2) * ||s - zeta + u||^2,
 *
 * where u is the Lagrange dual vector and mu is a penalty parameter.
 *
 * mADMM iteratively minimizes this function by the following coordinate descent steps until convergence:
 *
 * 1. Update zeta: zeta^{k+1} = (X.T * X + mu * I)^-1 * {X.T * f + mu * (s^k - u^k)}
 * 2. Update s: s^{k+1}_gj = C * max(0, 1 + lambda1 / (mu * || q^k_gj ||)) * q^k_gj,
 * where q^k = zeta^{k+1} - u^k, and C = mu / (mu + lambda2).
 * 3. Update u (dual variable): u^{k+1} = u^k + mu * (s^{k+1} - zeta^{l+1}),
 *
 * The steps above are implemented by the `update_b()`, `update_zeta()`, `update_s()`, and `update_u()` functions.
 *
 * This class also allows applying lower bound constraints (zeta_min <= zeta_i)
 * by adding the following steps to the iteration cycle:
 * 4. Update t: t_i^{k+1} = max(zeta_min, zeta_i^k - v_i^k)
 * 5. Update v: v^{k+1} = v^k + nu * (t^{k+1} - zeta^{k+1})
 *
 * where t is a slack vector, v is a Lagrange dual vector, and nu is a penalty parameter for the lower bound constraint.
 * Steps 4 and 5 are implemented by `update_t()` and `update_v()`.
 * A single cycle of updates is performed by `iterate()`.
 *
 * The entire process is handled by the `solve()` function.
 */
class mADMM : public ADMM {

	size_t	size1_mag_ = 0;	// Number of magnetic data points (rows of f)
	size_t	size1_grv_ = 0;	// Number of gravity data points (rows of g)

	size_t	m_ = 0;		// Total number of data points (m = dim(f) + dim(g))
	size_t	n_= 0;		// Number of model variables (n = 2 * number of grid cells)

	// input data
	double	*f_ = NULL;		// Magnetic anomaly data vector
	double	*g_ = NULL;		// Gravity anomaly data vector

	// transform matrix
	double	*X_ = NULL;		// Magnetic kernel matrix
	double	*Y_ = NULL;		// Gravity kernel matrix

	// weighting
	double	*wx_ = NULL;	// Weighting vector for the magnetic kernel matrix
	double	*wy_ = NULL;	// Weighting vector for the gravity kernel matrix

	double	*beta_ = NULL;	// Magnetization model vector
	double	*rho_ = NULL;	// Density model vector

	double	*beta_prev_ = NULL;	// Backup of magnetization from previous iteration
	double	*rho_prev_ = NULL;	// Backup of density from previous iteration

	double	*cx_ = NULL;	// Pre-computed vector: X.T * f
	double	*cy_ = NULL;	// Pre-computed vector: Y.T * g
	double	*bx_ = NULL;	// Intermediate vector for beta update: cx + mu * (s + u)[1:M] + nu * (t + v)[1:M]
	double	*by_ = NULL;	// Intermediate vector for rho update: cy + mu * (s + u)[M:2M] + nu * (t + v)[M:2M]
	double	*CXi_ = NULL;	// Pre-computed inverse matrix: inv(X * X.T + (mu + nu) * I)
	double	*CYi_ = NULL;	// Pre-computed inverse matrix: inv(Y * Y.T + (mu + nu) * I)

	double	residual_ = 0.;

	// Temporal vectors used to update the beta and rho based on the SMW formula.
	double	*tmp1_mag_ = NULL;
	double	*tmp2_mag_ = NULL;
	double	*tmp1_grv_ = NULL;
	double	*tmp2_grv_ = NULL;

public:
	mADMM () { }
	mADMM (double lambda1, double lambda2, double mu);
	~mADMM ();

	double	*getModel (ModelType type);

	double	*getDepthWeights (ModelType type);

	void		setupProblem (size_t size1_mag, size_t size1_grv, size_t size2,
			 		  double *f, double *g, double *X, double *Y, bool normalize, double nu, double *lower);

	// Solves the problem by running ADMM iterations until convergence.
	size_t	solve (const double tol, const size_t maxiter, bool verbos = false);

	double	getResidual () const { return residual_; }

	double	*recoverData (FieldType type);

protected:
	void		initialize_variables (); // Initializes all ADMM variables to zero.

	// ADMM iteration steps
	void		update_bx (); // Updates the intermediate vector `bx` for the beta update.
	void		update_by (); // Updates the intermediate vector `by` for the rho update.
	void		update_b ();  // Calls update_bx() and update_by().
	void		update_zeta (); // Updates the combined model vector `zeta = [beta; rho]`.
	void		update_s (); // Updates the slack vector for the regularization penalty.
	void		update_t (); // Updates the slack vector for the lower bound constraint.
	void		update_u (); // Updates the Lagrange dual variable for the regularization penalty.
	void		update_v (); // Updates the Lagrange dual variable for the lower bound constraint.

	void		iterate ();	// Performs a single full cycle of ADMM updates.

	double	*get_magnetization (); // Get the magnetization model vector
	double	*get_density (); // Get the density model vector

	double	eval_residuals (); // Evaluates the maximum of primal and dual residuals for convergence check.

	void		compute_Ci ();  // Computes the intermediate matrices CXi and CYi.
	// Computes the beta update step using the Sherman-Morrison-Woodbury formula.
	double	*update_beta_using_SMW (double coef, size_t m, size_t n, double *X, double *CXi, double *bx, double *y1, double *y2);
	// Computes the rho update step using the Sherman-Morrison-Woodbury formula.
	double	*update_rho_using_SMW (double coef, size_t m, size_t n, double *Y, double *CYi, double *by, double *y1, double *y2);

	// Applies the soft-thresholding operator for the L2 + Group Lasso penalty.
	double	soft_threshold (double gamma, double lambda);

private:
};

#endif
