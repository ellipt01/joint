#ifndef _JOINT_H_
#define _JOINT_H_

typedef enum {
	DATA_TYPE_MAGNETIC = 0,
	DATA_TYPE_GRAVITY = 1,
	DATA_TYPE_NONE = -1
} data_type;

/***
	Joint: class of magnetic and gravity joint inversion

	This class provides APIs for performing magnetic and gravity inversion
	based on L2 norm and group lasso regularization.

	*** flow of the processing ***
	call the following methods:

	prepare (): processing inline options and reading settings file.
	start (): run joint inversion.

	The following methods exports calculation results:
	export_results (): output derived magnetization and density models,
					   and recoverd magnetic and gravity anomalies.
***/

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

	// get instances
	double		get_alpha () { return _alpha_; }
	double		get_lambda () { return _lambda_; }

	double		get_lambda1 () { return _lambda1_; }
	double		get_lambda2 () { return _lambda2_; }

	double		get_tolerance () { return _tolerance_; }
	size_t		get_maxiter () { return _maxiter_; }

	// display usage
	void		usage ();

	// processing inline options and reading settings file
	void		prepare (int argc, char **argv);

	// performs inversion
	size_t		start (bool normalize);
	size_t		restart ();

	// return major and dual residuals
	double		residual ();

	// recover input anomalies using derived model
	void		recover (mm_real *f, mm_real *g);

	// gtet model instances derived by the inversion
	mm_real		*get_beta (); // = zeta[:m]
	mm_real		*get_rho ();  // = zeta[m:]

	// displays the settings specified
	// in the inline options and the settings file.
	void		fwrite_inline (FILE *stream);
	void		fwrite_settings (FILE *stream);

	// export calculation results
	void		export_results ();

protected:
	// read inline options
	void		_read_inline_ (int argc, char **argv);
	// read settings file
	void		_fread_settings_ (FILE *stream);

	// read data files
	void		_read_data_ ();
	// register lower bound constraint
	void		_set_lower_bounds_ ();
	// register terrai data
	void		_set_surface_ (size_t c, double *zsurf);

	void		_simeq_ ();

	void		_set_mag_ (double exf_inc, double exf_dec, double mgz_inc, double mgz_dec, data_array *data);
	void		_set_grv_ (data_array *data);

	void		_fwrite_model_ (FILE *fp);
	void		_export_weights_ ();
	void		_export_matrices_ ();
	
private:
	// initialize
	void		__init__ ();
	// count number of data
	size_t		__count__ (FILE *fp);
	// read terrain file
	double		*__read_terrain__ (FILE *fp, const size_t c);

};

#endif
