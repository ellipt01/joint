#include <mgcal.h>
#include <mmreal.h>

/*** super class of kernel ***/
class Kernel {

protected:
	size_t		_nx_;
	size_t		_ny_;
	size_t		_nz_;

	double		*_xx_;
	double		*_yy_;
	double		*_zz_;

	grid		*_grd_;
	data_array	*_data_;

	mm_real		*_K_;

	mgcal_func	*_func_;

	vector3d	*_exf_;
	vector3d	*_mgz_;

public:
	Kernel () { _init_ (); }

	void	set_range (size_t nx, size_t ny, size_t nz, double *xx, double *yy, double *zz);
	void	set_range (size_t nx, size_t ny, size_t nz, double *xx, double *yy, double *zz, const double ll);
	void	set_surface (double *zsurf) { grid_set_surface (_grd_, zsurf); }
	void	set_data (data_array *array);

	grid	*get_grid () { return _grd_; }

	void	eval ();
	mm_real	*get ();

	void	fwrite (FILE *stream, mm_real *model);
	void	fwrite (FILE *stream, mm_real *model, const char *format);

protected:
	void	_init_ ();
};

/*** magnetic kernel ***/
class MagKernel : public Kernel
{

public:
	MagKernel (double inc, double dec);
	MagKernel (double exf_inc, double exf_dec, double mgz_inc, double mgz_dec);

	void	set_exf (double inc, double dec);
	void	set_mgz (double inc, double dec);
};

/*** gravity kernel ***/
class GravKernel : public Kernel
{

public:
	GravKernel ();
};

