#ifndef _KERNEL_H_
#define _KERNEL_H_

/*** super class of kernel ***/
class Kernel {

protected:
	// number of grid
	size_t		_nx_;
	size_t		_ny_;
	size_t		_nz_;

	// range of model space
	double		*_xx_;
	double		*_yy_;
	double		*_zz_;

	// grid object for model space
	grid		*_grd_;
	// observed data
	data_array	*_data_;

	// kernel matrix
	mm_real		*_K_;

	// mgcal func for computing kernel matrix
	mgcal_func	*_func_;

	// unit vector parallel to external field and magnetization
	vector3d	*_exf_;
	vector3d	*_mgz_;

public:
	Kernel () { __init__ (); }

	// set range and creates grid for model space
	void	set_range (size_t nx, size_t ny, size_t nz, double *xx, double *yy, double *zz);
	void	set_range (size_t nx, size_t ny, size_t nz, double *xx, double *yy, double *zz, const double ll);
	// set surface topography
	void	set_surface (double *zsurf) { grid_set_surface (_grd_, zsurf); }
	// set observed data
	void	set_data (data_array *array);

	// return grid object for model space
	grid	*get_grid () { return _grd_; }

	// compute and return kernel matrix
	mm_real	*get ();

	// fwrite model
	void	fwrite (FILE *stream, mm_real *model);
	void	fwrite (FILE *stream, mm_real *model, const char *format);

protected:
	// evaluate kernel matrix
	void	_eval_ ();

private:
	// initializer
	void	__init__ ();
};

/*** magnetic kernel ***/
class MagKernel : public Kernel
{

public:
	MagKernel (double inc, double dec);
	MagKernel (double exf_inc, double exf_dec, double mgz_inc, double mgz_dec);

	// set unit vector parallel to external field and magnetization
	void	set_exf (double inc, double dec);
	void	set_mgz (double inc, double dec);
};

/*** gravity kernel ***/
class GravKernel : public Kernel
{

public:
	GravKernel ();
};

#endif
