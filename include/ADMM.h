#include <mkl_blas.h>
#include <mkl_lapack.h>
#include <mmreal.h>

/*** A class implements the Alternative Direction Method of Multiplier (ADMM) ***/
class ADMM {

	size_t	_size1_;	// num of row of the kernel matrix
	size_t	_size2_;	// num of columns

	double	_alpha_;	// mixing ratio of L2 / group lasso
	double	_lambda_;	// regularization parameter

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
	mm_real	*_b_;		// = c + mu * (t + v)[1:size2]
	mm_real	*_Ci_;		// = (X * X.T + mu * I)^-1

	double	_residual_;

public:
	ADMM () {}
	ADMM (double alpha, double lambda, double mu);
	ADMM (double alpha, double lambda, double mu, double nu, mm_real *lower);

	void	set_params (double alpha, double lambda);

	void	simeq (mm_real *f, mm_real *X, bool normalize);

	mm_real	*get_zeta ();

	mm_real	*get_s () { return _s_; }
	mm_real	*get_t () { return _t_; }
	mm_real	*get_u () { return _u_; }
	mm_real	*get_v () { return _v_; }
	mm_real	*get_Ci () { return _Ci_; } // magnetic field

	size_t	start (const double tol, const size_t maxiter);
	size_t	restart (const double tol, const size_t maxiter);

	double	residual () { return _residual_; }

	mm_real	*recover ();

	void	fread_Cinv (FILE *fp);
	void	fwrite_Cinv (FILE *fp);

protected:
	void	_initialize_ (); 	// initialize zeta, s, and u

	void	_calc_Ci_ ();  // compute CXi = inv(X * X.T + mu * I)
	void	_update_b_ (); // update b = X.T * f + mu (s + u)

	// ADMM iterations
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
	// compute inv(X * X.T + mu * I)
	mm_real	*_Cinv_SMW_ (double mu, mm_real *K);
	// compute (I - X.T * Ci * X) * b / mu
	mm_real	*_inv_SMW_ (double mu, mm_real *X, mm_real *Ci, mm_real *b);

private:
	void	__init__ ();
};

