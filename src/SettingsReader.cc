#include <iostream>
#include <cstring>

#include "SettingsReader.h"

void
SettingsReader::_init_ (void)
{
	_ngrid_specified_ = false;
	_range_specified_ = false;
	_incdec_specified_ = false;
	_invparams_specified_ = false;
	_pparam_specified_ = false;
	_apply_lower_bound_ = false;
	// default values
	_tolerance_ = 1.e-5;
	_maxiter_ = 10000;
}

SettingsReader::SettingsReader ()
{
	_init_ ();
}

void
SettingsReader::fread (FILE *stream)
{
	char	buf[BUFSIZ];

	while (fgets (buf, BUFSIZ, stream) != NULL) {
		if (buf[0] == '#') continue;
		if (buf[0] == ' ' || buf[0] == '\t' || buf[0] == '\n') continue;

		char	*p = strrchr (buf, ':');
		if (p == NULL) continue;
		do {
			p++;
		} while (*p == ' ' || *p == '\t');

		switch (buf[0]) {
			case '1': // number of grids
				sscanf (p, "%lu, %lu, %lu", &_nx_, &_ny_, &_nz_);
				_ngrid_specified_ = true;
				break;
			case '2': // x,y,z range
				sscanf (p, "%lf,%lf,%lf,%lf,%lf,%lf",
					&_xrange_[0], &_xrange_[1], &_yrange_[0], &_yrange_[1], &_zrange_[0], &_zrange_[1]);
				_range_specified_ = true;
				break;
			case '3': // inc, dec
				sscanf (p, "%lf,%lf,%lf,%lf", &_exf_inc_, &_exf_dec_, &_mgz_inc_, &_mgz_dec_);
				_incdec_specified_ = true;
				break;
			case '4': // tolerance, maxiter
				sscanf (p, "%lf,%lu", &_tolerance_, &_maxiter_);
				_invparams_specified_ = true;
				break;
			case '5': // mu
				_mu_ = (double) atof (p);
				_pparam_specified_ = true;
				break;
			case '6': // mu
				sscanf (p, "%lf,%lf,%lf", &_nu_, &_beta0_, &_rho0_);
				if (_nu_ > 0.) _apply_lower_bound_ = true;
				break;
			default:
				break;
		}
	}
}

void
SettingsReader::fwrite (FILE *stream)
{
	fprintf (stream, "\n");
	if (_ngrid_specified_) fprintf (stream, "number of grid:\t%ld/%ld/%ld\n", _nx_, _ny_, _nz_);
	if (_range_specified_)
		fprintf (stream, "x, y, z range:\t%.4f/%.4f/%.4f/%.4f/%.4f/%.4f\n", 
				_xrange_[0], _xrange_[1], _yrange_[0], _yrange_[1], _zrange_[0], _zrange_[1]);
	if (_incdec_specified_) {
		fprintf (stream, "exf:inc, dec:\t%.4f/%.4f\n", _exf_inc_, _exf_dec_);
		fprintf (stream, "mag:inc, dec:\t%.4f/%.4f\n", _mgz_inc_, _mgz_dec_);
	}
	if (_invparams_specified_) fprintf (stream, "tol, maxiter:\t%.2e/%ld\n", _tolerance_, _maxiter_);
	if (_pparam_specified_) fprintf (stream, "mu:\t\t%.4f\n", _mu_);
	if (_apply_lower_bound_) fprintf (stream, "nu:\t\t%.4f: lower bounds = %.4f, %.4f\n",
		_nu_, _beta0_, _rho0_);
}

void
SettingsReader::ngrid (size_t *nx, size_t *ny, size_t *nz)
{
	if (!_ngrid_specified_) throw std::runtime_error ("ngrid is not specified");
	if (nx) *nx = _nx_;
	if (ny) *ny = _ny_;
	if (nz) *nz = _nz_;
}

void
SettingsReader::range (double x[], double y[], double z[])
{
	if (!_range_specified_) throw std::runtime_error ("number of grid is not specified");
	if (x) {
		x[0] = _xrange_[0];
		x[1] = _xrange_[1];
	}
	if (y) {
		y[0] = _yrange_[0];
		y[1] = _yrange_[1];
	}
	if (z) {
		z[0] = _zrange_[0];
		z[1] = _zrange_[1];
	}
}

void
SettingsReader::incdec (double *exf_inc, double *exf_dec, double *mgz_inc, double *mgz_dec)
{
	if (!_incdec_specified_)
		throw std::runtime_error ("inclination and declination are not specified");
	if (exf_inc) *exf_inc = _exf_inc_;
	if (exf_dec) *exf_dec = _exf_dec_;
	if (mgz_inc) *mgz_inc = _mgz_inc_;
	if (mgz_dec) *mgz_dec = _mgz_dec_;
}

void
SettingsReader::invparams (double *tol, size_t *maxiter)
{
	if (tol) *tol = _tolerance_;
	if (maxiter) *maxiter = _maxiter_;
}

void
SettingsReader::pparam (double *mu)
{
	if (!_pparam_specified_)
		throw std::runtime_error ("penalty parameter is not specified");
	if (mu) *mu = _mu_;
}

void
SettingsReader::lower_bound (double *nu, double *beta0, double *rho0)
{
	if (!_apply_lower_bound_) return;
	if (nu) *nu = _nu_;
	if (beta0) *beta0 = _beta0_;
	if (rho0) *rho0 = _rho0_;
}


