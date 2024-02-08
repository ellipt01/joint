#ifndef _JOINT_H_
#define _JOINT_H_

class Joint : public Inversion
{

	char		_fn_mag_[256];
	char		_fn_grv_[256];

	double		_beta_lower_;
	double		_rho_lower_;

	data_array	*_magdata_;
	data_array	*_grvdata_;

	mm_real		*_f_;
	mm_real		*_K_;

	double		_scale_; // scale for gravity data
	mm_real		*_g_;
	mm_real		*_G_;

	mADMM		*_admm_;

public:
	Joint () { __init__ (); }
	
	void	usage ();
	
	void	prepare (int argc, char **argv);

	size_t	start (bool normalize);
	size_t	restart ();

	double	residual ();
	void	recover (mm_real *f, mm_real *g);

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

	void	simeq ();

	void	set_mag (double exf_inc, double exf_dec, double mgz_inc, double mgz_dec, data_array *data);
	void	set_grv (data_array *data);

	void	fwrite_model (FILE *fp);
	void	export_matrix ();
	
private:
	void	__init__ ();

};

#endif
