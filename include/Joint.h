#include "Kernel.h"
#include "mADMM.h"

class Joint
{
	double		_alpha_;
	double		_lambda_;
	double		_mu_;
	double		_nu_;
	mm_real		*_lower_;

	double		_exf_inc_;
	double		_exf_dec_;

	double		_mag_inc_;
	double		_mag_dec_;

	size_t		_nx_;
	size_t		_ny_;
	size_t		_nz_;

	double		*_xrange_;
	double		*_yrange_;
	double		*_zrange_;
	double		*_zsurf_;

	data_array	*_magdata_;
	MagKernel	*_magker_; 

	data_array	*_grvdata_;
	GravKernel	*_grvker_; 

	mm_real		*_f_;
	mm_real		*_K_;

	double		_scale_; // scale for gravity data
	mm_real		*_g_;
	mm_real		*_G_;

	size_t		_niter_;
	mADMM		*_admm_;

public:
	Joint (double alpha, double lambda, double mu, double nu, mm_real *lower);

	void	set_params (double alpha, double lambda, double mu);
	void	set_range (size_t nx, size_t ny, size_t nz, double x[], double y[], double z[]);
	void	set_surface (size_t c, double *zsurf);

	void	set_mag (double inc, double dec, data_array *data);
	void	set_mag (double exf_inc, double exf_dec, double mag_inc, double mag_dec, data_array *data);
	void	set_grv (data_array *data);

	size_t	start (const double tol, size_t maxiter, bool normalize);
	size_t	restart (const double tol, size_t maxiter);
	double	residual ();

	void	recover (mm_real *f, mm_real *g);

	mm_real	*get_beta ();
	mm_real	*get_rho ();

	mm_real	*get_f () { return _f_; }
	mm_real	*get_K () { return _K_; }
	mm_real	*get_g () { return _g_; }
	mm_real	*get_G () { return _G_; }

	void	fwrite (FILE *fp);

private:
	void	__init__ ();
};

