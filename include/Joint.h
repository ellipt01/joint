#ifndef _JOINT_H_
#define _JOINT_H_

typedef enum {
	DATA_TYPE_MAGNETIC = 0,
	DATA_TYPE_GRAVITY = 1,
	DATA_TYPE_NONE = -1
} data_type;

class Joint
{
	char		*_toolname_;

	// settings file name
	char		_fn_settings_[256];
	// terrain file name
	char		*_fn_ter_;

	// external field inclination and declination
	double		_exf_inc_;
	double		_exf_dec_;

	// magnetization vector inclination and declination
	double		_mgz_inc_;
	double		_mgz_dec_;

	// input file names
	char		_fn_mag_[256];
	char		_fn_grv_[256];

	// mixing ratio of L1-L2 norm penalty
	double		_alpha_;
	// regularization parameter
	double		_lambda_;

	// regularization parameters for L1 and L2 norm
	double		_lambda1_;
	double		_lambda2_;

	// number of subsurface grid
	size_t		_nx_;
	size_t		_ny_;
	size_t		_nz_;

	// range of the model space
	double		_xrange_[2]; // Eastward positive
	double		_yrange_[2]; // Northward positive
	double		_zrange_[2]; // Upward positive
	double		*_zsurf_; // store terrain data, read from _fn_ter_

	// penalty parameter for s = zeta
	double		_mu_;
	// penalty parameter for lower bound constraint
	double		_nu_;
	double		_beta_lower_; // magnetization lower bound
	double		_rho_lower_; // density lower bound
	mm_real		*_lower_; // store the above bounds

	// magnetic and gravity observed data
	data_array	*_magdata_;
	data_array	*_grvdata_;

	// magnetic data vector and kernel matrix
	mm_real		*_f_;
	mm_real		*_K_;

	// gravity data vector and kernel matrix
	mm_real		*_g_;
	mm_real		*_G_;

	// scale for magnetic and gravity data
	double		_scale_; // = |g|_{infinity} / |f|_{infinity}

	// kernel matrix calculator
	MagKernel	*_magker_;
	GravKernel	*_grvker_;

	// ADMM object
	mADMM		*_admm_;

	// tolerance and number of max iteration
	double		_tolerance_;
	size_t		_maxiter_;

	// counter for number of iterations
	size_t		_niter_;

	bool		_export_matrix_;
	bool		_verbos_;

public:
	Joint () { __init__ (); }

	double		get_alpha () { return _alpha_; }
	double		get_lambda () { return _lambda_; }

	double		get_lambda1 () { return _lambda1_; }
	double		get_lambda2 () { return _lambda2_; }

	double		get_tolerance () { return _tolerance_; }
	size_t		get_maxiter () { return _maxiter_; }

	void		usage ();
	
	void		prepare (int argc, char **argv);

	size_t		start (bool normalize);
	size_t		restart ();

	double		residual ();
	void		recover (mm_real *f, mm_real *g);

	mm_real		*get_beta (); // = zeta[:m]
	mm_real		*get_rho ();  // = zeta[m:]

	void		fwrite_inline (FILE *stream);
	void		fwrite_settings (FILE *stream);

	void		export_results ();

protected:
	void		read_inline (int argc, char **argv);
	void		fread_settings (FILE *stream);

	void		read_data ();
	void		set_lower_bounds ();
	void		set_surface (size_t c, double *zsurf);

	void		simeq ();

	void		set_mag (double exf_inc, double exf_dec, double mgz_inc, double mgz_dec, data_array *data);
	void		set_grv (data_array *data);

	void		fwrite_model (FILE *fp);
	void		export_weight ();
	void		export_matrix ();
	
private:
	void		__init__ ();

};

#endif
