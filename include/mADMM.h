#ifndef _MADMM_H_
#define _MADMM_H_

/***
	mADMM: class implements the Alternative Direction Method of Multiplier (ADMM)
	for magnetic and gravity joint inversion

	This class is designed for solving a linear problem
	with L2 norm and group lasso combined regularization:
	
	min (1/2) ||b - Z * zeta||^2 + lambda1 * sum_j=0^M ||zeta_gj|| + lambda2 ||zeta||^2/2,
	
	where b = [f; g], and Z = [X;Y], and f and g are the observed magnetic and gravity data,
	and K and Y are the magnetic and gravity kernel matrices.
	zeta = [beta; rho] and beta and rho are magnetization and density model vectors, respectively.

	||a|| indicates the Euclidean norm of a vector a, and the second term of this
	objective function is the group lasso penalty.
	M is the number of the grid cells dividing the subsurface model space.

	mADMM search an optimal solution of the above by replacing the problem as following:
	
	min (1/2) ||b - Z * zeta||^2 + lambda1 * sum_j=0^M ||s_gj|| + lambda2 ||s||^2/2 s.t. s = zeta
	
	where s is a slack vector.
	The augmented Lagrange function of this problem is
	
	L = (1/2) ||b - Z * zeta||^2 + lambda1 * sum_j=0^M ||s_gj|| + lambda2 ||s||^2/2
		+ (mu / 2) * ||s - zeta + u||^2,

	where u is the Lagrange dual vector.
	
	mADMM solves this problem by the following coordinate descent,
	e.g. repeat the following cycle until solution converged:
	
	1. update zeta: zeta^{k+1} = (Z.T * Z + mu * I)^-1 * {Z.T * b + mu * (s^k - u^k)}
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
class mADMM : public ADMM {

	// input data
	mm_real	*_f_;		// magnetic anomaly
	mm_real	*_g_;		// gravity anomaly

	// transform matrix
	mm_real	*_X_;		// magnetic kernel
	mm_real	*_Y_;		// gravity kernel

	// weighting
	mm_real	*_wx_;
	mm_real	*_wy_;

	mm_real	*_beta_;	// magnetization
	mm_real	*_rho_;		// density

	mm_real	*_beta_prev_;	// backup magnetization
	mm_real	*_rho_prev_;	// backup density

	mm_real	*_cx_;		// = X.T * f
	mm_real	*_cy_;		// = Y.T * g
	mm_real	*_bx_;		// = cx + mu * (s + u)[1:size2] + nu * (t + v)[1:size2]
	mm_real	*_by_;		// = cy + mu * (s + u)[size2:2*size2] + nu * (t + v)[size2:2*size2]
	mm_real	*_CXi_;		// = (X * X.T + (mu + nu) * I)^-1
	mm_real	*_CYi_;		// = (Y * Y.T + (mu + nu) * I)^-1

	double	_residual_;

public:
	mADMM () { __init__ (); }
	mADMM (double lambda1, double lambda2, double mu);
	mADMM (double lambda1, double lambda2, double mu, double nu, mm_real *lower);

	mm_real	*get_beta ();
	mm_real	*get_rho ();

	mm_real	*get_wx () { return _wx_; } // weight for magkernel matrix
	mm_real	*get_wy () { return _wy_; } // weight for grvkernel matrix
	mm_real	*get_CXi () { return _CXi_; } // inverse for magkernel matrix
	mm_real	*get_CYi () { return _CYi_; } //             grvkernel matrix

	void	simeq (mm_real *f, mm_real *g, mm_real *X, mm_real *Y, bool normalize);

	size_t	start (const double tol, const size_t maxiter) { return start (tol, maxiter, false); };
	size_t	start (const double tol, const size_t maxiter, bool verbos);
	size_t	restart (const double tol, const size_t maxiter);

	double	residual () { return _residual_; }

	void	recover (mm_real *f, mm_real *g);

protected:
	void	_initialize_ (); // initialize zeta, t, and v

	void	_calc_Ci_ ();   /* compute CXi = (X * X.T + mu * I)^-1
						     and     CYi = (Y * Y.T + mu * I)^-1 */
	void	_update_bx_ (); // update bx = X.T * f + mu (t + v)[1:M]
	void	_update_by_ (); // update by = Y.T * g + mu (t + v)[M:2M]
	void	_update_b_ ();  // update bx and by

	// ADMM iteration
	void	_update_zeta_ ();
	void	_update_s_ ();
	void	_update_t_ ();

	void	_update_u_ ();
	void	_update_v_ ();

	void	_one_cycle_ ();	// perform ADMM iteration at once

	double	_eval_residuals_ (); // for convergence check

	// soft threshold for L2 + group Lasso
	double	_soft_threshold_ (double gamma, double lambda);

private:
	void	__init__ ();
};

#endif
