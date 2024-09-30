#ifndef _MADMM_H_
#define _MADMM_H_

/***
	mADMM: subclass of ADMM, and this class implements the Alternative Direction Method
	       of Multiplier (ADMM) algorithm

	This class is designed for solving a linear problem with L1-L2 norm regularization:
	
	min (1/2) ||b - Z * zeta||^2 + lambda1 * sum_j=0^M ||zeta_gj|| + lambda2 ||zeta||^2/2,
	
	where b = [f; g], and Z = [X;Y], and f and g are the observed magnetic and gravity data,
	and K and Y are the magnetic and gravity kernel matrices.
	zeta = [beta; rho] and beta and rho are magnetization and density model vectors, respectively.
	The second term of the objective function is group lasso penalty,
	and M is the number of the grid cells dividing the subsurface model space.

	Instead to minimize this objective function, mADMM searchs an optimal solution
	of the following constrained optimization problem:

	min (1/2) ||b - Z * zeta||^2 + lambda1 * sum_j=0^M ||s_gj|| + lambda2 ||s||^2/2 s.t. s = zeta
	
	where s is a slack vector.
	The augmented Lagrange function of this problem is
	
	L = (1/2) ||f - X * zeta||^2 + lambda1 * |s| + lambda2 ||s||^2/2
		+ (mu / 2) * ||s - zeta + u||^2,

	where u is the Lagrange dual vector, and mu is a penalty parameter.
	
	mADMM searhcs a model zeta which minimizes this function by the following
	coordinate descent method, e.g. repeat the following cycle until converged:
	
	1. update zeta: zeta^{k+1} = (X.T * X + mu * I)^-1 * {X.T * f + mu * (s^k - u^k)}
	2. update s: s^{k+1}_gj = C * max (0, 1 + lambda1 / (mu * || q^k_gj ||)) * q^k_gj,
	   where q^k = zeta^{k+1} - v^k, and C = mu / (mu + lambda2).
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
class mADMM : public ADMM {

	size_t	_size1_mag_;	// = dim(f)
	size_t	_size1_grv_;	// = dim(g)

	size_t	_m_;		// m = dim(f) + dim(g)
	size_t	_n_;		// n = size(X, 2) + size(Y, 2)

	// input data
	double	*_f_;		// magnetic anomaly
	double	*_g_;		// gravity anomaly

	// transform matrix
	double	*_X_;		// magnetic kernel matrix
	double	*_Y_;		// gravity kernel matrix

	// weighting
	double	*_wx_;
	double	*_wy_;

	double	*_beta_;	// magnetization
	double	*_rho_;		// density

	double	*_beta_prev_;	// backup magnetization
	double	*_rho_prev_;	// backup density

	double	*_cx_;		// = X.T * f
	double	*_cy_;		// = Y.T * g
	double	*_bx_;		// = cx + mu * (s + u)[1:size2] + nu * (t + v)[1:size2]
	double	*_by_;		// = cy + mu * (s + u)[size2:2*size2] + nu * (t + v)[size2:2*size2]
	double	*_CXi_;		// = (X * X.T + (mu + nu) * I)^-1
	double	*_CYi_;		// = (Y * Y.T + (mu + nu) * I)^-1

	double	_residual_;

public:
	mADMM () { __init__ (); }
	mADMM (double lambda1, double lambda2, double mu);

	double	*get_beta ();
	double	*get_rho ();

	double	*get_wx () { return _wx_; } // get depth weighting of magnetic kernel matrix
	double	*get_wy () { return _wy_; } //                        gravity kernel matrix

	void	simeq (size_t size1_mag, size_t size1_grv, size_t size2,
				   double *f, double *g, double *X, double *Y, bool normalize, double nu, double *lower);

	// start ADMM iteration until model converged
	size_t	start (const double tol, const size_t maxiter) { return start (tol, maxiter, false); };
	size_t	start (const double tol, const size_t maxiter, bool verbos);

	double	residual () { return _residual_; }

	void	recover (double *f, double *g);

protected:
	void	_initialize_ (); // initialize zeta, t, and v

	// ADMM iteration
	void	_update_bx_ (); // update bx = X.T * f + mu * (s + u)[1:M] + nu * (t + v)[1:M]
	void	_update_by_ (); // update by = Y.T * g + mu * (s + u)[M:2M] + nu * (t + v)[M:2M]
	void	_update_b_ ();  // update bx and by
	void	_update_zeta_ (); // update zeta = [beta; rho]
	void	_update_s_ (); // update slack vector
	void	_update_t_ (); // update slack vector for bound constraint
	void	_update_u_ (); // update Lagrange dual
	void	_update_v_ (); // update Lagrange dual for bound constraint

	void	_one_cycle_ ();	// perform ADMM iteration at once: updates b, zeta, s, t, t and v

	double	_eval_residuals_ (); // evaluate maximum of primal and dual residuals

	void	_calc_Ci_ ();   // compute CXi = (X * X.T + coef * I)^-1 and CYi = (Y * Y.T + coef * I)^-1
	// compute beta = (I - X.T * CXi * X) * bx / coef
	double	*_eval_beta_SMW_ (double coef, size_t m, size_t n, double *X, double *CXi, double *bx);
	// compute rho = (I - Y.T * CYi * Y) * by / coef
	double	*_eval_rho_SMW_ (double coef, size_t m, size_t n, double *Y, double *CYi, double *by);

	// soft threshold for L2 + group Lasso
	double	_soft_threshold_ (double gamma, double lambda);

private:
	void	__init__ ();
};

#endif
