#ifndef _KERNEL_H_
#define _KERNEL_H_

/***
	@class Kernel
	@brief Base class for calculating kernel matrices.

	This class provides an interface for computing the kernel matrices
	used in magnetic and gravity inversion.

	@section workflow Typical Workflow
	1. set_range(): Specify the model space dimensions and discretization.
	2. set_data():  Register the observed data (locations, altitudes, anomalies).
	3. get():       Compute the kernel matrix (if not already computed) and return it.
***/
class Kernel {

protected:
	// Grid dimensions
	size_t	nx_ = 0;
	size_t	ny_ = 0;
	size_t	nz_ = 0;

	// Model space coordinates
	double	*xx_ = NULL;
	double	*yy_ = NULL;
	double	*zz_ = NULL;

	// Grid object representing the model space
	grid			*grd_ = NULL;
	// Observed data points
	data_array	*data_ = NULL;

	// The computed kernel matrix (m x n)
	size_t	m_ = 0;
	size_t	n_ = 0;
	double	*K_ = NULL;

	// Pointer to the specific kernel computation function
	mgcal_func	*func_ = NULL;

public:
	Kernel () { }
	~Kernel ();

	// Sets the model space dimensions and creates the grid.
	void		setRange (size_t nx, size_t ny, size_t nz, double *xx, double *yy, double *zz, const double ll = 0.);
	// Sets the surface topography for the model grid.
	void		setSurface (double *zsurf) { grid_set_surface (grd_, zsurf); }
	// Sets the observed measurement data.
	void		setData (data_array *array);

	// Computes and returns the kernel matrix.
	virtual double	*get () = 0;

	// Returns a pointer to the grid object.
	grid		*getGrid () const { return grd_; }

	// Writes the model data to a file stream.
	void		fwrite (FILE *stream, double *model);
	void		fwrite (FILE *stream, double *model, const char *format);

private:

};

/***
	@class MagKernel
	@brief Computes the kernel matrix for magnetic inversion.
***/
class MagKernel : public Kernel
{
	// Unit vectors for external magnetic field and magnetization direction
	vector3d	*exf_ = NULL;
	vector3d	*mgz_ = NULL;

public:
	MagKernel (double inc, double dec);
	MagKernel (double exf_inc, double exf_dec, double mgz_inc, double mgz_dec);
	double	*get ();
};

/***
	@class GravKernel
	@brief Computes the kernel matrix for gravity inversion.
***/
class GravKernel : public Kernel
{

public:
	GravKernel ();
	double	*get ();

};

#endif
