#include <iostream>

#include "mgcal.h"
#include "gravity.h"
#include "Kernel.h"

/*** public methods ***/
// Destructor
Kernel::~Kernel ()
{
	if (grd_) grid_free (grd_);
	delete [] K_;
}

// Sets the model space dimensions and creates the grid.
void
Kernel::setRange (size_t nx, size_t ny, size_t nz, double *xx, double *yy, double *zz, const double ll)
{
	grd_ = grid_new (nx, ny, nz, xx, yy, zz);
	n_ = grd_->n;
	if (ll > 0.) grid_stretch_at_edge (grd_, ll);
}

// Sets the observed measurement data.
void
Kernel::setData (data_array *data)
{
	data_ = data;
	m_ = data_->n;
}

/*** magnetic kernel ***/
MagKernel::MagKernel (double inc, double dec)
{
	exf_ = vector3d_new_with_geodesic_poler (1., inc, dec);
	mgz_ = vector3d_new_with_geodesic_poler (1., inc, dec);
	// Set the kernel function to total_force_prism.
	func_ = mgcal_func_new (total_force_prism, NULL);
}

MagKernel::MagKernel (double exf_inc, double exf_dec, double mgz_inc, double mgz_dec)
{
	exf_ = vector3d_new_with_geodesic_poler (1., exf_inc, exf_dec);
	mgz_ = vector3d_new_with_geodesic_poler (1., mgz_inc, mgz_dec);
	func_ = mgcal_func_new (total_force_prism, NULL);
}

// Computes and returns the magnetic kernel matrix.
double *
MagKernel::get ()
{
	if (data_ == NULL)
		throw std::runtime_error ("Observed data has not been set. Call setData() before get().");
	if (grd_ == NULL)
		throw std::runtime_error ("Model range has not been set. Call setRange() before get().");
	if (exf_ == NULL)
		throw std::runtime_error ("External field direction has not been set.");
	if (mgz_ == NULL)
		throw std::runtime_error ("Magnetization direction has not been set.");

	if (K_ != NULL) delete [] K_;
	K_ = new double [m_ * n_];
	kernel_matrix_set (K_, data_, grd_, mgz_, exf_, func_);

	return K_;
}

/*** gravity kernel ***/
GravKernel::GravKernel ()
{
	// Set the kernel function to grav_prism.
	func_ = mgcal_func_new (grav_prism, NULL);
}

// Computes and returns the gravity kernel matrix.
double *
GravKernel::get ()
{
	if (data_ == NULL)
		throw std::runtime_error ("Observed data has not been set. Call setData() before get().");
	if (grd_ == NULL)
		throw std::runtime_error ("Model range has not been set. Call setRange() before get().");

	if (K_ != NULL) delete [] K_;
	K_ = new double [m_ * n_];
	kernel_matrix_set (K_, data_, grd_, NULL, NULL, func_);

	return K_;
}

