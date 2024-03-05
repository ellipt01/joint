#include <iostream>

#include "mgcal.h"
#include "mmreal.h"
#include "Kernel.h"
#include "gravity.h"

/*** public methods ***/
// set range and creates grid for model space
void
Kernel::set_range (size_t nx, size_t ny, size_t nz, double *xx, double *yy, double *zz)
{
	_nx_ = nx;
	_ny_ = ny;
	_nz_ = nz;
	_grd_ = grid_new (_nx_, _ny_, _nz_, xx, yy, zz);
	grid_stretch_at_edge (_grd_, 1000.);
}

void
Kernel::set_range (size_t nx, size_t ny, size_t nz, double *xx, double *yy, double *zz, const double ll)
{
	_nx_ = nx;
	_ny_ = ny;
	_nz_ = nz;
	_grd_ = grid_new (_nx_, _ny_, _nz_, xx, yy, zz);
	if (ll > 0.) grid_stretch_at_edge (_grd_, ll);
}

// set observed data
void
Kernel::set_data (data_array *data)
{
	_data_ = data;
}

// compute and return kernel matrix
mm_real *
Kernel::get ()
{
	if (_K_ == NULL) _eval_ ();
	return _K_;
}

// fwrite model
void
Kernel::fwrite (FILE *stream, mm_real *model)
{
	fwrite_grid_with_data (stream, _grd_, model->data, "%.4e\t%.4e\t%.4e\t%.8e");
}

void
Kernel::fwrite (FILE *stream, mm_real *model, const char *format)
{
	fwrite_grid_with_data (stream, _grd_, model->data, format);
}

/*** protected methods ***/
// evaluate kernel matrix
void
Kernel::_eval_ ()
{
	if (_data_ == NULL) throw std::runtime_error ("please set data_array.");
	if (_grd_ == NULL)   throw std::runtime_error ("please set range.");

	size_t	m = _data_->n;
	size_t	n = _grd_->n;
	_K_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, n, m * n);
	kernel_matrix_set (_K_->data, _data_, _grd_, _mgz_, _exf_, _func_);
}

/*** private methods ***/
// initializer
void
Kernel::__init__ ()
{
	_nx_ = -1;
	_ny_ = -1;
	_nz_ = -1;

	_grd_ = NULL;
	_data_ = NULL;

	_K_ = NULL;

	_func_ = NULL;

	_exf_ = NULL;
	_mgz_ = NULL;
}

/*** magnetic kernel ***/
MagKernel::MagKernel (double inc, double dec)
{
	_exf_ = vector3d_new_with_geodesic_poler (1., inc, dec);
	_mgz_ = vector3d_new_with_geodesic_poler (1., inc, dec);
	// set mgcal_func to total_force_prism
	_func_ = mgcal_func_new (total_force_prism, NULL);
}

MagKernel::MagKernel (double exf_inc, double exf_dec, double mgz_inc, double mgz_dec)
{
	_exf_ = vector3d_new_with_geodesic_poler (1., exf_inc, exf_dec);
	_mgz_ = vector3d_new_with_geodesic_poler (1., mgz_inc, mgz_dec);
	_func_ = mgcal_func_new (total_force_prism, NULL);
}

void
MagKernel::set_exf (double inc, double dec)
{
	_exf_ = vector3d_new_with_geodesic_poler (1., inc, dec);
}

void
MagKernel::set_mgz (double inc, double dec)
{
	_mgz_ = vector3d_new_with_geodesic_poler (1., inc, dec);
}

/*** gravity kernel ***/
GravKernel::GravKernel ()
{
	// set mgcal func to grav_prism
	_func_ = mgcal_func_new (grav_prism, NULL);
}

