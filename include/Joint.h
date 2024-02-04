#include "Kernel.h"
#include "mADMM.h"

typedef enum {
	DATA_TYPE_MAGNETIC = 0,
	DATA_TYPE_GRAVITY = 1,
	DATA_TYPE_NONE = -1
} data_type;

class Joint
{

	char		*_toolname_;

	char		_fn_mag_[256];
	char		_fn_grv_[256];
	char		_fn_settings_[256];
	char		*_fn_ter_;

	double		_alpha_;
	double		_log10_lambda_;
	double		_lambda_;
	double		_mu_;

	double		_nu_;
	double		_beta_lower_;
	double		_rho_lower_;
	mm_real		*_lower_;

	double		_exf_inc_;
	double		_exf_dec_;

	double		_mgz_inc_;
	double		_mgz_dec_;

	size_t		_nx_;
	size_t		_ny_;
	size_t		_nz_;

	double		_xrange_[2];
	double		_yrange_[2];
	double		_zrange_[2];
	double		*_zsurf_;

	double		_tolerance_;
	size_t		_maxiter_;

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

	bool		_export_matrix_;

public:
	Joint () { __init__ (); }
	
	void	usage ();
	
	void	prepare (int argc, char **argv);

	size_t	start (bool normalize);
	size_t	restart ();

	double	residual ();
	void	recover (mm_real *f, mm_real *g);

	double	get_tolerance () { return _tolerance_; }
	size_t	get_maxiter () { return _maxiter_; }

	mm_real	*get_beta ();
	mm_real	*get_rho ();

	void	fwrite_inline (FILE *stream);
	void	fwrite_settings (FILE *stream);

	void	export_results ();

protected:
	void	read_inline (int argc, char **argv);
	void	fread_settings (FILE *stream);

	void	read_data ();
	void	set_lower_bounds ();
	void	set_surface (size_t c, double *zsurf);

	void	simeq ();

	void	set_mag (double inc, double dec, data_array *data);
	void	set_mag (double exf_inc, double exf_dec, double mag_inc, double mag_dec, data_array *data);
	void	set_grv (data_array *data);

	void	fwrite_model (FILE *fp);
	void	export_matrix ();
	
private:
	void	__init__ ();

};

