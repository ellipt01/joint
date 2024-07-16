#ifndef _ADMM_H_
#define _ADMM_H_

/***
	class implements the Alternative Direction Method of Multiplier (ADMM) 

	This class is designed for solving a linear problem
	with L1-L2 norm regularization:
	
	min (1/2) ||f - X * zeta||^2 + lambda1 * |zeta| + lambda2 ||zeta||^2/2

	where ||a|| indicates the Euclidean norm, and |a| is the absolute norm of a vector a.

	ADMM search an optimal solution of the above by replacing the problem as following:
	
	min (1/2) ||f - X * zeta||^2 + lambda1 * |s| + lambda2 ||s||^2/2 s.t. s = zeta
	
	where s is a slack variable(vector) and u is a Lagrange coefficient(vector).
	
	ADMM solves this problem by the following coordinate descent,
	e.g. repeat the following cycle until solution converged:
	
	1. update zeta: zeta^{k+1} = (X.T * X + mu * I)^-1 * (X.T * f + mu * (s^k - u^k)
	2. update s: s^{k+1} = (mu / (mu + lambda2)) * S(zeta^{k+1} - u^k, lambda1 / mu),
	   where S() is a soft threshold operator of the L1-L2 norm lenalty,
	3. update u by the multiplyer method: u^{k+1} = u^k + mu * (s^{k+1} - zeta^{l+1}),

	where a^k is the vector obtained by the k-th iteration.
	the steps 1, 2, and 3 of above are implemented by
	_update_b_ () and _update_zeta (),
	_update_s_ (), 
	and
	_update_u_ ().
***/
class ADMM {
protected:
	size_t	_size1_;	// num of row of the kernel matrix
	size_t	_size2_;	// num of columns

	double	_lambda1_;	// regularization parameter for group lasso penalty
	double	_lambda2_;	// regularization parameter for L2 norm penalty

	double	_mu_;		// penalty parameter for regularization

	double	_nu_;		// penalty parameter for lower bounds
	mm_real	*_lower_;	// lower bound
	bool	_apply_lower_bound_;	// if lower bound constraint is applied, takes true

	// input data
	mm_real	*_f_;		// observed data

	// transform matrix
	mm_real	*_X_;		// kernel matrix

	// weighting
	mm_real	*_w_;		// reciprocal of sensitivity weighting

	mm_real	*_zeta_;	// model vector

	mm_real	*_zeta_prev_;	// backup for model vector

	/*** regularization penalty ***/
	mm_real	*_s_;		// slack variable
	mm_real	*_u_;		// Lagrange dual

	/**** lower bound constraint ***/
	mm_real	*_t_;		// slack variable
	mm_real	*_v_;		// Lagrange dual

	mm_real	*_c_;		// = X.T * f
	mm_real	*_b_;		// = c + coef * (t + v)
	mm_real	*_Ci_;		// = (X * X.T + coef * I)^-1

	double	_residual_;

public:
	ADMM () { __init__ (); }
	ADMM (double lambda1, double lambda2, double mu);
	ADMM (double lambda1, double lambda2, double mu, double nu, mm_real *lower);

	// get model vector zeta
	mm_real	*get_zeta ();

	// get depth weightings
	mm_real	*get_w () { return _w_; } // weight for kernel matrix

	// set regulatization parameters
	void	set_params (double lambda1, double lambda2);
	void	simeq (mm_real *f, mm_real *X, bool normalize);

	// start ADMM iteration until model converged
	size_t	start (const double tol, const size_t maxiter) { return start (tol, maxiter, false); };
	size_t	start (const double tol, const size_t maxiter, bool verbos);
	size_t	restart (const double tol, const size_t maxiter);

	// get residual
	double	residual () { return _residual_; }

	// recover the input data
	mm_real	*recover ();

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
	mm_real	*_normalize_ (mm_real *K);
	// compute inv(C) using cholesky decomposition
	void	_cholinv_ (mm_real *C);
	// compute inv(C) using LU decomposition
	void	_LUinv_ (mm_real *C);
	// compute inv(X * X.T + coef * I)
	mm_real	*_Cinv_SMW_ (double coef, mm_real *K);
	// compute (I - X.T * Ci * X) * b / coef
	mm_real	*_inv_SMW_ (double coef, mm_real *X, mm_real *Ci, mm_real *b);

private:
	void	__init__ ();
};

#endif
