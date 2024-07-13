#include <iostream>
#include <cstring>
#include <unistd.h>
#include <cfloat>

#include "mgcal.h"
#include "mmreal.h"

#include "Kernel.h"
#include "Inversion.h"

/*** public methods ***/

/*** protected methods ***/
void
Inversion::set_surface (size_t c, double *zsurf)
{
	if (_nx_ <= 0 || _ny_ <= 0)
		throw std::runtime_error ("range dose not specified. Call set_range() before.");
	if (_nx_ * _ny_ != c)
		throw std::runtime_error ("dim(terrain) is not match with dim(f) and dim(g)");
	_zsurf_ = new double [c];
	for (size_t i = 0; i < c; i++) _zsurf_[i] = zsurf[i];
}

/*** private methods ***/
void
Inversion::__init__ ()
{
	strcpy (_fn_settings_, "settings.par");

	_fn_ter_ = NULL;

	_nx_ = -1;
	_ny_ = -1;
	_nz_ = -1;

	_mu_ = 1.;

	_zsurf_ = NULL;

	_tolerance_ = 1.e-3;
	_maxiter_ = 1000;

	_export_matrix_ = false;
}


