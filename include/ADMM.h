#ifndef ADMM_H_
#define ADMM_H_

/*
 * ADMM: An implementation of the Alternating Direction Method of Multipliers (ADMM) algorithm.
 *
 * This class solves a linear problem with L1-L2 norm regularization, defined as:
 * min (1/2) ||f - X * zeta||^2 + lambda1 * |zeta| + lambda2 * ||zeta||^2 / 2
 *
 * The ADMM algorithm transforms this into a constrained optimization problem:
 * min (1/2) ||f - X * zeta||^2 + lambda1 * |s| + lambda2 * ||s||^2 / 2  s.t. s = zeta
 *
 * It then minimizes the augmented Lagrangian function:
 * L = (1/2) ||f - X * zeta||^2 + lambda1 * |s| + lambda2 * ||s||^2 / 2 + (mu / 2) * ||s - zeta + u||^2
 *
 * The algorithm iteratively updates the variables (zeta, s, u) until convergence.
 *
 * Optional lower bound constraints (zeta_min <= zeta_i) can be applied by adding
 * two extra update steps for slack vector 't' and Lagrange dual 'v'.
 */
class ADMM {
protected:
	size_t	size1_ = 0;	// Number of data points (rows of X)
	size_t	size2_ = 0;	// Number of model variables (columns of X)

	double	lambda1_ = 0.;	// Regularization parameter for the L1-norm penalty
	double	lambda2_ = 0.;	// Regularization parameter for the L2-norm penalty

	double	mu_ = 0.;		// Penalty parameter for the regularization term
	
	double	nu_ = 0.;		// Penalty parameter for the lower bound constraint
	double	*lower_ = NULL;	// Array of lower bound values for each variable
	bool		apply_lower_bound_ = false;	// True if lower bound constraint is active

	// input data
	double	*f_ = NULL;	// Observed data vector

	// transform matrix
	double	*X_ = NULL;	// Kernel matrix (transformation matrix)

	// weighting
	double	*w_ = NULL;	// Reciprocal of sensitivity weighting, used for normalization

	double	*zeta_ = NULL;	// The primary model vector being solved for

	double	*zeta_prev_ = NULL;	// Backup of the model vector from the previous iteration

	/*** regularization penalty ***/
	double	*s_ = NULL;	// Slack variable for the L1-L2 penalty term
	double	*u_ = NULL;	// Lagrange dual variable for the L1-L2 penalty term

	/**** lower bound constraint ***/
	double	*t_ = NULL;	// Slack variable for the lower bound constraint
	double	*v_ = NULL;	// Lagrange dual variable for the lower bound constraint

	double	*c_ = NULL;	// Pre-computed vector: X.T * f
	double	*b_ = NULL;	// Pre-computed vector: c + mu * (s - u)
	double	*Ci_ = NULL;	// Pre-computed inverse matrix: inv(X * X.T + mu * I)

	double	residual_ = 0.;

public:
	ADMM () { }
	ADMM (double lambda1, double lambda2, double mu);
	~ADMM ();

	// Set regulatization parameters
	void		setLambdas (double lambda1, double lambda2);
	void		setupProblem (size_t size1, size_t size2, double *f, double *X, bool normalize, double nu = 0., double *lower = NULL);

	// Get model vector zeta
	double	*getModel ();

	// Get depth weightings
	double	*getDepthWeights () { return w_; } // Returns the depth weighting vector

	// Solve problem by running ADMM iteration until model converged
	size_t	solve (const double tol, const size_t maxiter, bool verbos = false);

	// Get residual
	double	getResidual () const { return residual_; }

	// Recover the input data
	double	*recoverData ();

protected:
	// Initializes zeta, s, and u to default values
	void		initialize_variables ();

	// Computes and stores the pre-factorized inverse matrix Ci
	void		compute_Ci ();

	// ADMM iterations
	void		update_b (); // Updates the intermediate vector 'b' for the zeta update step
	void		update_zeta ();
	void		update_s ();
	void		update_t ();

	// update duals
	void		update_u ();
	void		update_v ();

	// Performs a single full cycle of the ADMM updates
	void		iterate ();

	// Evaluates the primal and dual residuals for convergence check
	double	eval_residuals ();

	// soft threshold for L2 + L1
	double	soft_threshold (double gamma, double lambda);

	// normalize matrix K and return weighting
	double	*normalize_matrix (size_t m, size_t n, double *K);
	// compute inv(C) using cholesky decomposition
	void		cholinv (char uplo, size_t n, double *C);
	// compute inv(C) using LU decomposition
	void		LUinv (size_t m, size_t n, double *C);
	// Computes inv(X * X.T + coef * I) for the Sherman-Morrison-Woodbury formula.
	double	*compute_Cinv_for_SMW (double coef, size_t m, size_t n, double *K);
	// Computes the zeta update step (I - X.T * Ci * X) * b / coef using the Sherman-Morrison-Woodbury formula.
	double	*eval_zeta_using_SMW (double coef, size_t m, size_t n, double *K, double *Ci, double *b);

private:
};

#endif
