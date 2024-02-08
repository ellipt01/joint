#ifndef _INVERSION_H_
#define _INVERSION_H_

typedef enum {
	DATA_TYPE_MAGNETIC = 0,
	DATA_TYPE_GRAVITY = 1,
	DATA_TYPE_NONE = -1
} data_type;

class Inversion
{
protected:
	char		*_toolname_;

	char		_fn_settings_[256];
	char		*_fn_ter_;

	double		_alpha_;
	double		_log10_lambda_;
	double		_lambda_;
	double		_mu_;

	double		_nu_;
//	double		_zeta_lower_;
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

	MagKernel	*_magker_; 
	GravKernel	*_grvker_; 

	size_t		_niter_;

	bool		_export_matrix_;

public:
	Inversion () { __init__ (); }
	
	double	get_tolerance () { return _tolerance_; }
	size_t	get_maxiter () { return _maxiter_; }

	virtual void	usage () = 0;
	virtual void	prepare (int argc, char **argv) = 0;

	virtual void	fwrite_inline (FILE *stream) = 0;
	virtual void	fwrite_settings (FILE *stream) = 0;

protected:
	void	set_surface (size_t c, double *zsurf);

	virtual void	read_inline (int argc, char **argv) = 0;
	virtual void	fread_settings (FILE *stream) = 0;

	virtual void	simeq () = 0;

private:
	void	__init__ ();

};

#endif
