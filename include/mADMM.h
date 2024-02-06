#ifndef _MADMM_H_
#define _MADMM_H_

class mADMM : public ADMM {

	size_t	_size1_;	// num of row of the kernel matrix
	size_t	_size2_;	// num of columns

	double	_alpha_;	// mixing ratio of L2 / group lasso
	double	_lambda_;	// regularization parameter

	double	_mu_;		// penalty parameter for regularization

	double	_nu_;		// penalty parameter for lower bound
	mm_real	*_lower_;	// lower bound
	bool	_apply_lower_bound_;

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

	mm_real	*_s_;		// slack variable for regularization
	mm_real	*_t_;		// slack variable for lower bound

	mm_real	*_u_;		// Lagrange dual
	mm_real	*_v_;

	mm_real	*_cx_;		// = X.T * f
	mm_real	*_cy_;		// = Y.T * g
	mm_real	*_bx_;		// = cx + mu * (s + u)[1:size2] + nu * (t + v)[1:size2]
	mm_real	*_by_;		// = cy + mu * (s + u)[size2:2*size2] + nu * (t + v)[size2:2*size2]
	mm_real	*_CXi_;		// = (X * X.T + (mu + nu) * I)^-1
	mm_real	*_CYi_;		// = (Y * Y.T + (mu + nu) * I)^-1

	double	_residual_;

public:
	mADMM (double alpha, double lambda, double mu, double nu, mm_real *lower);

	void	set_params (double alpha, double lambda);

	void	simeq (mm_real *f, mm_real *g, mm_real *X, mm_real *Y, bool normalize);

	mm_real	*get_beta ();
	mm_real	*get_rho ();

	mm_real	*get_t () { return _s_; }
	mm_real	*get_v () { return _u_; }
	mm_real	*get_CXi () { return _CXi_; } // magnetic field
	mm_real	*get_CYi () { return _CYi_; } // gravity

	size_t	start (const double tol, const size_t maxiter);
	size_t	restart (const double tol, const size_t maxiter);

	double	residual () { return _residual_; }

	void	recover (mm_real *f, mm_real *g);

	void	fread_Cinv (FILE *fp);
	void	fwrite_Cinv (FILE *fp);

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
