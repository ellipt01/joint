#ifndef _ADMM_H_
#define _ADMM_H_

/***
	ADMM: class implements the Alternative Direction Method of Multiplier (ADMM) algorithm

	This class is designed for solving a linear problem with L1-L2 norm regularization:
	
	min (1/2) ||f - X * zeta||^2 + lambda1 * |zeta| + lambda2 ||zeta||^2/2

	where ||a|| indicates the Euclidean norm, and |a| is the absolute norm of a vector a.

	Instead to minimize this objective function, ADMM searchs an optimal solution
	of the following constrained optimization problem:
	
	min (1/2) ||f - X * zeta||^2 + lambda1 * |s| + lambda2 ||s||^2/2 s.t. s = zeta,
	
	where s is a slack vector.
	The augmented Lagrange function of this problem is
	
	L = (1/2) ||f - X * zeta||^2 + lambda1 * |s| + lambda2 ||s||^2/2
		+ (mu / 2) * ||s - zeta + u||^2,

	where u is the Lagrange dual vector, and mu is a penalty parameter.
	
	ADMM searhcs a model zeta which minimizes this function by the following
	coordinate descent method, e.g. repeat the following cycle until converged:
	
	1. update zeta: zeta^{k+1} = (X.T * X + mu * I)^-1 * {X.T * f + mu * (s^k - u^k)}
	2. update s: s^{k+1} = (mu / (mu + lambda2)) * S(zeta^{k+1} - u^k, lambda1 / mu),
	   where S() is a soft threshold operator of the L1-L2 norm penalty,
	3. update u by the multiplyer method: u^{k+1} = u^k + mu * (s^{k+1} - zeta^{l+1}),

	where a^k is the vector obtained by the k-th iteration.
	the steps 1, 2, and 3 of above are implemented by
	_update_b_ () and _update_zeta (),
	_update_s_ (), and _update_u_ ().

	This class allows to apply lower bound constraint
	
	zeta_min <= zeta_i
	
	This is implemented by add the following cycle to the ADMM iteration:
	4. update t: t_i^{k+1} = max(max (zeta_min, zeta_i^k - v_i^k)
	5. update v: v^{k+1} = v^k + nu * (t^{k+1} - zeta^{k+1})

	where t is a slack vector, v is a Lagrange vector, and nu is a penalty parameter
	for the lower bound constraint.
	step 4 and 5 are implemented by _update_t_ () and _update_v_ (),
	and one cycle of the update is implemented by _one_cycle_ ().
	
	The steps 1-6 are repeated until solution converged, and this process is
	implemented by start ().

***/
class ADMM {
protected:
	size_t	_size1_;	// num of data (=dim(f)), row of the kernel matrix X
	size_t	_size2_;	// num of subsurface grid cells, nim of columns of X

	double	_lambda1_;	// regularization parameter for group lasso penalty
	double	_lambda2_;	// regularization parameter for L2 norm penalty

	double	_mu_;		// penalty parameter for regularization

	double	_nu_;		// penalty parameter for lower bounds
	double	*_lower_;	// lower bound
	bool	_apply_lower_bound_;	// if lower bound constraint is applied, takes true

	// input data
	double	*_f_;		// observed data

	// transform matrix
	double	*_X_;		// kernel matrix

	// weighting
	double	*_w_;		// reciprocal of sensitivity weighting

	double	*_zeta_;	// model vector

	double	*_zeta_prev_;	// backup for model vector

	/*** regularization penalty ***/
	double	*_s_;		// slack variable
	double	*_u_;		// Lagrange dual

	/**** lower bound constraint ***/
	double	*_t_;		// slack variable
	double	*_v_;		// Lagrange dual

	double	*_c_;		// = X.T * f
	double	*_b_;		// = c + coef * (t + v)
	double	*_Ci_;		// = (X * X.T + coef * I)^-1

	double	_residual_;

public:
	ADMM () { __init__ (); }
	ADMM (double lambda1, double lambda2, double mu);

	// get model vector zeta
	double	*get_zeta ();

	// get depth weightings
	double	*get_w () { return _w_; } // get depth weighting of kernel matrix

	// set regulatization parameters
	void	set_params (double lambda1, double lambda2);
	void	simeq (size_t size1, size_t size2, double *f, double *X, bool normalize, double nu, double *lower);

	// start ADMM iteration until model converged
	size_t	start (const double tol, const size_t maxiter) { return start (tol, maxiter, false); };
	size_t	start (const double tol, const size_t maxiter, bool verbos);

	// get residual
	double	residual () { return _residual_; }

	// recover the input data
	double	*recover ();

protected:
	void	_initialize_ (); 	// initialize zeta, s, and u

	void	_calc_Ci_ ();  // compute CXi = inv(X * X.T + coef * I)

	// ADMM iterations
	void	_update_b_ (); // update b = X.T * f + coef * (s + u)
	void	_update_zeta_ ();
	void	_update_s_ ();
	void	_update_t_ ();

	// update duals
	void	_update_u_ ();
	void	_update_v_ ();

	void	_one_cycle_ ();	// perform ADMM iteration at once

	double	_eval_residuals_ (); // for convergence check

	// soft threshold for L2 + L1
	double	_soft_threshold_ (double gamma, double lambda);

	// normalize matrix K and return weighting
	double	*_normalize_ (size_t m, size_t n, double *K);
	// compute inv(C) using cholesky decomposition
	void	_cholinv_ (char uplo, size_t n, double *C);
	// compute inv(C) using LU decomposition
	void	_LUinv_ (size_t m, size_t n, double *C);
	// compute inv(X * X.T + coef * I)
	double	*_Cinv_SMW_ (double coef, size_t m, size_t n, double *K);
	// compute (I - X.T * Ci * X) * b / coef
	double	*_eval_zeta_SMW_ (double coef, size_t m, size_t n, double *K, double *Ci, double *b);

private:
	void	__init__ ();
};

#endif
