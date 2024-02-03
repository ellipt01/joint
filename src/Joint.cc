#include <iostream>
#include <cstring>
#include <unistd.h>

#include "Joint.h"

/*** public methods ***/
Joint::Joint (double alpha, double lambda, double mu, double nu, mm_real *lower)
{
	__init__ ();
	set_params (alpha, lambda, mu);

	_nu_ = nu;
	_lower_ = lower;
}

void
Joint::set_params (double alpha, double lambda, double mu)
{
	_alpha_ = alpha;
	_lambda_ = lambda;
	_mu_ = mu;
}

void
Joint::set_range (size_t nx, size_t ny, size_t nz, double x[], double y[], double z[])
{
	_nx_ = nx;
	_ny_ = ny;
	_nz_ = nz;

	_xrange_ = new double [2];
	_xrange_[0] = x[0];
	_xrange_[1] = x[1];

	_yrange_ = new double [2];
	_yrange_[0] = y[0];
	_yrange_[1] = y[1];

	_zrange_ = new double [2];
	_zrange_[0] = z[0];
	_zrange_[1] = z[1];
}

void
Joint::set_surface (size_t c, double *zsurf)
{
	if (_nx_ <= 0 || _ny_ <= 0) throw std::runtime_error ("range dose not specified. Call set_range() before.");
	if (_nx_ * _ny_ != c) throw std::runtime_error ("dim(terrain) is not match with dim(f) and dim(g)");
	_zsurf_ = new double [c];
	for (size_t i = 0; i < c; i++) _zsurf_[i] = zsurf[i];
}

void
Joint::set_mag (double inc, double dec, data_array *data)
{
	if (_xrange_ == NULL || _yrange_ == NULL || _zrange_ == NULL)
		throw std::runtime_error ("range is not specified. Call set_range()");

	_exf_inc_ = inc;
	_exf_dec_ = dec;
	_mag_inc_ = inc;
	_mag_dec_ = dec;

	_magdata_ = data;
	size_t	n = _magdata_->n;
	_f_ = mm_real_view_array (MM_REAL_DENSE, MM_REAL_GENERAL, n, 1, n, _magdata_->data);

	_magker_ = new MagKernel (_exf_inc_, _exf_dec_, _mag_inc_, _mag_dec_);
	_magker_->set_range (_nx_, _ny_, _nz_, _xrange_, _yrange_, _zrange_, 1000.);
	if (_zsurf_) _magker_->set_surface (_zsurf_);
	_magker_->set_data (_magdata_);	
	_K_ = _magker_->get ();
}

void
Joint::set_mag (double exf_inc, double exf_dec, double mag_inc, double mag_dec, data_array *data)
{
	if (_xrange_ == NULL || _yrange_ == NULL || _zrange_ == NULL)
		throw std::runtime_error ("range is not specified. Call set_range()");

	_exf_inc_ = exf_inc;
	_exf_dec_ = exf_dec;
	_mag_inc_ = mag_inc;
	_mag_dec_ = mag_dec;

	_magdata_ = data;
	size_t	n = _magdata_->n;
	_f_ = mm_real_view_array (MM_REAL_DENSE, MM_REAL_GENERAL, n, 1, n, _magdata_->data);

	_magker_ = new MagKernel (_exf_inc_, _exf_dec_, _mag_inc_, _mag_dec_);
	_magker_->set_range (_nx_, _ny_, _nz_, _xrange_, _yrange_, _zrange_, 1000.);
	if (_zsurf_) _magker_->set_surface (_zsurf_);
	_magker_->set_data (_magdata_);	
	_K_ = _magker_->get ();
}

void
Joint::set_grv (data_array *data)
{
	if (_xrange_ == NULL || _yrange_ == NULL || _zrange_ == NULL)
		throw std::runtime_error ("range is not specified. please start set_range()");

	_grvdata_ = data;
	size_t	n = _grvdata_->n;
	_g_ = mm_real_view_array (MM_REAL_DENSE, MM_REAL_GENERAL, n, 1, n, _grvdata_->data);

	_grvker_ = new GravKernel ();
	_grvker_->set_range (_nx_, _ny_, _nz_, _xrange_, _yrange_, _zrange_, 1000.);
	if (_zsurf_) _grvker_->set_surface (_zsurf_);
	_grvker_->set_data (_grvdata_);	
	_G_ = _grvker_->get ();
}

size_t
Joint::start (const double tol, size_t maxiter, bool normalize)
{
	if (_f_ == NULL || _K_ == NULL) throw std::runtime_error ("magnetic equation is not specified. Call set_grv()");
	if (_g_ == NULL || _G_ == NULL) throw std::runtime_error ("gravity equation is not specified. Call set_grv()");

	// scale gravity anomaly
	size_t	ifamax = mm_real_iamax (_f_);
	double	famax = fabs (_f_->data[ifamax]);
	size_t	igamax = mm_real_iamax (_g_);
	double	gamax = fabs (_g_->data[igamax]);
	_scale_ = famax / gamax;
std::cout << "scale = " << _scale_ << std::endl;
	mm_real_xj_scale (_g_, 0, _scale_);

	if (!_admm_) {
		_admm_ = new mADMM (_alpha_, _lambda_, _mu_, _nu_, _lower_);
		_admm_->simeq (_f_, _g_, _K_, _G_, normalize);
	}
	return _admm_->start (tol, maxiter);
}

size_t
Joint::restart (const double tol, size_t maxiter)
{
	if (!_admm_) throw std::runtime_error ("admm object has not yet been instantiated. Call start() first.");

	_admm_->set_params (_alpha_, _lambda_);
	return _admm_->restart (tol, maxiter);
}

double
Joint::residual ()
{
	if (!_admm_) throw std::runtime_error ("admm object has not yet been instantiated. Call start() first.");
	return _admm_->residual ();
}

void
Joint::recover (mm_real *f, mm_real *g)
{
	_admm_->recover (f, g);
	mm_real_xj_scale (g, 0, 1. / _scale_);
}

mm_real	*
Joint::get_beta ()
{
	mm_real	*_beta_ = _admm_->get_beta ();
	mm_real	*beta = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _beta_->m, 1, _beta_->m);
	mm_real_memcpy (beta, _beta_);
	return beta;
}

mm_real *
Joint::get_rho ()
{
	mm_real	*_rho_ = _admm_->get_rho ();
	mm_real	*rho = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _rho_->m, 1, _rho_->m);
	mm_real_memcpy (rho, _rho_);
	mm_real_xj_scale (rho, 0, 1. / _scale_);
	return rho;
}

void
Joint::fwrite (FILE *fp)
{
	if (_magker_ == NULL) throw std::runtime_error ("magnetic kernel is not specified. Call set_mag()");
	if (_grvker_ == NULL) throw std::runtime_error ("gravity kernel is not specified. Call set_grv()");

	mm_real	*beta = get_beta ();
	mm_real	*rho = get_rho ();

	grid		*grd = _magker_->get_grid ();
	vector3d	*pos = vector3d_new (0., 0., 0.);
	for (size_t k = 0; k < grd->n; k++) {
		grid_get_nth (grd, k, pos, NULL);
		fprintf (fp, "%.4e\t%.4e\t%.4e\t%.8e\t%.8e\n", pos->x, pos->y, pos->z, beta->data[k], rho->data[k]);
	}

	mm_real_free (beta);
	mm_real_free (rho);

#ifdef DEBUG
	FILE	*fp_grd = fopen ("grid.data", "w");
	if (fp) {
		fwrite_grid (fp_grd, grd);
		fclose (fp_grd);
	}
#endif

	delete [] pos;
}

/*** private methods ***/
void
Joint::__init__ ()
{
	_nx_ = -1;
	_ny_ = -1;
	_nz_ = -1;

	_xrange_ = NULL;
	_yrange_ = NULL;
	_zrange_ = NULL;
	_zsurf_ = NULL;

	_magdata_ = NULL;
	_grvdata_ = NULL;

	_magker_ = NULL; 
	_grvker_ = NULL; 

	_f_ = NULL;
	_K_ = NULL;

	_scale_ = 1.;
	_g_ = NULL;
	_G_ = NULL;

	_admm_ = NULL;

}


