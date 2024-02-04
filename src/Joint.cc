#include <iostream>
#include <cstring>
#include <unistd.h>

#include "Joint.h"
#include "util.h"

/*** public methods ***/
void
Joint::usage ()
{
	fprintf (stderr, "USAGE: %s\n", _toolname_);
	fprintf (stderr, "       -f <magnetic anomaly filename>\n");
	fprintf (stderr, "       -g <gravitic anomaly filename>\n");
	fprintf (stderr, "       -l <log10(lambda)>\n");
	fprintf (stderr, "       -a <alpha>>\n");
	fprintf (stderr, "[optional]\n");
	fprintf (stderr, "       -t <terrain filename>\n");
	fprintf (stderr, "       -s <setting filename:default is settings.par>\n");
	fprintf (stderr, "       -x (export kernel matrices)\n");
	fprintf (stderr, "       -h (show this message)\n");
}

void
Joint::prepare (int argc, char **argv)
{
	_toolname_ = get_toolname (argv[0]);

	read_inline (argc, argv);

	FILE	*fp = fopen (_fn_settings_, "r");
	if (!fp) {
		char	msg[128];
		sprintf (msg, "cannot open setting file %s", _fn_settings_);
		throw std::runtime_error (msg);
	}
	fread_settings (fp);
	fclose (fp);
}

size_t
Joint::start (bool normalize)
{
	if (_f_ == NULL || _K_ == NULL || _g_ == NULL || _G_ == NULL) simeq ();

	// scale gravity anomaly
	size_t	ifamax = mm_real_iamax (_f_);
	double	famax = fabs (_f_->data[ifamax]);
	size_t	igamax = mm_real_iamax (_g_);
	double	gamax = fabs (_g_->data[igamax]);

	_scale_ = famax / gamax;
	std::cerr << "scale = " << _scale_ << std::endl;

	mm_real_xj_scale (_g_, 0, _scale_);

	if (!_admm_) {
		_admm_ = new mADMM (_alpha_, _lambda_, _mu_, _nu_, _lower_);
		_admm_->simeq (_f_, _g_, _K_, _G_, normalize);
	}
	return _admm_->start (_tolerance_, _maxiter_);
}

size_t
Joint::restart ()
{
	if (!_admm_) throw std::runtime_error ("admm object has not yet been instantiated. Call start() first.");

	_admm_->set_params (_alpha_, _lambda_);
	return _admm_->restart (_tolerance_, _maxiter_);
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
Joint::fwrite_inline (FILE *stream)
{
	fprintf (stream, "\n");
	fprintf (stream, "inut file:\t%s\t%s\n", _fn_mag_, _fn_grv_);
	if (_fn_ter_) fprintf (stream, "terrain file:\t%.s\n", _fn_ter_); 
	fprintf (stream, "lambda:\t%.4f\n", _lambda_);
	fprintf (stream, "alpha :\t%.4f\n", _alpha_);
	fprintf (stream, "export matrix: ");
	(_export_matrix_) ? fprintf (stream, "true\n") : fprintf (stream, "false\n");
}

void
Joint::fwrite_settings (FILE *stream)
{
	fprintf (stream, "\n");
	fprintf (stream, "number of grid:\t%ld/%ld/%ld\n", _nx_, _ny_, _nz_);
	fprintf (stream, "x, y, z range:\t%.4f/%.4f/%.4f/%.4f/%.4f/%.4f\n", 
			 _xrange_[0], _xrange_[1], _yrange_[0], _yrange_[1], _zrange_[0], _zrange_[1]);
	fprintf (stream, "exf:inc, dec:\t%.4f/%.4f\n", _exf_inc_, _exf_dec_);
	fprintf (stream, "mag:inc, dec:\t%.4f/%.4f\n", _mgz_inc_, _mgz_dec_);
	fprintf (stream, "tol, maxiter:\t%.2e/%ld\n", _tolerance_, _maxiter_);
	fprintf (stream, "mu:\t\t%.4f\n", _mu_);
	if (_nu_ > 0.) fprintf (stream, "nu:\t\t%.4f: lower bounds = %.4f, %.4f\n",
		_nu_, _beta_lower_, _rho_lower_);
}

void
Joint::export_results ()
{
	FILE *fp;
	
	fp = fopen ("model.data", "w");
	if (!fp) throw std::runtime_error ("cannot open file model.data");

	fwrite_model (fp);
	fclose (fp);

	mm_real	*fr = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _magdata_->n, 1, _magdata_->n);
	mm_real	*gr = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, _grvdata_->n, 1, _grvdata_->n);

	recover (fr, gr);

	fp = fopen ("recover_mag.data", "w");
	if (!fp) throw std::runtime_error ("cannot open file recover_mag.data");

	fwrite_data_array_with_data (fp, _magdata_, fr->data, "%.4f\t%.4f\t%.4f\t%.4f");
	mm_real_free (fr);
	fclose (fp);

	fp = fopen ("recover_grv.data", "w");
	if (!fp) throw std::runtime_error ("cannot open file recover_grv.data");

	fwrite_data_array_with_data (fp, _grvdata_, gr->data, "%.4f\t%.4f\t%.4f\t%.4f");
	mm_real_free (gr);
	fclose (fp);

	return;
}

/*** protected methods ***/
void
Joint::read_inline (int argc, char **argv)
{
	bool	fn_mag_specified = false;
	bool	fn_grv_specified = false;
	bool	lambda_specified = false;
	bool	alpha_specified = false;

	char	opt;
	while ((opt = getopt (argc, argv, ":f:g:l:a:s:t:xh")) != -1) {
		switch (opt) {
			case 'f':
				strcpy (_fn_mag_, optarg);
				fn_mag_specified = true;
				break;

			case 'g':
				strcpy (_fn_grv_, optarg);
				fn_grv_specified = true;
				break;

			case 'l':
				_log10_lambda_ = (double) atof (optarg);
				_lambda_ = pow (10., _log10_lambda_);
				lambda_specified = true;
				break;

			case 'a':
				_alpha_ = (double) atof (optarg);
				alpha_specified = true;
				break;

			case 's':
				strcpy (_fn_settings_, optarg);
				break;

			case 't':
				_fn_ter_ = new char [256];
				strcpy (_fn_ter_, optarg);
				
				break;

			case 'x':
				_export_matrix_ = true;
				break;

			case 'h':
				usage ();

			default:
				break;
		}
	}
	if (!fn_mag_specified) throw std::runtime_error ("magnetic anomaly data filename is not specified");
	if (!fn_grv_specified) throw std::runtime_error ("gravitic anomaly data filename is not specified");
	if (!lambda_specified) throw std::runtime_error ("lambda is not specified");
	if (!alpha_specified) throw std::runtime_error ("alpha is not specified");
}

void
Joint::fread_settings (FILE *stream)
{
	bool	ngrid_specified = false;
	bool	range_specified = false;
	bool	incdec_specified = false;

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
				ngrid_specified = true;
				break;
			case '2': // x,y,z range
				sscanf (p, "%lf,%lf,%lf,%lf,%lf,%lf",
					&_xrange_[0], &_xrange_[1], &_yrange_[0], &_yrange_[1], &_zrange_[0], &_zrange_[1]);
				range_specified = true;
				break;
			case '3': // inc, dec
				sscanf (p, "%lf,%lf,%lf,%lf", &_exf_inc_, &_exf_dec_, &_mgz_inc_, &_mgz_dec_);
				incdec_specified = true;
				break;
			case '4': // tolerance, maxiter: default is 1.e-3, 1000
				sscanf (p, "%lf,%lu", &_tolerance_, &_maxiter_);
				break;
			case '5': // mu: default is 1.0
				_mu_ = (double) atof (p);
				break;
			case '6': // lower bounds: default is nu = -1 (lower bounds not applied) 
				sscanf (p, "%lf,%lf,%lf", &_nu_, &_beta_lower_, &_rho_lower_);
				break;
			default:
				break;
		}
	}
	if (!ngrid_specified) throw std::runtime_error ("number of grid not specified");
	if (!range_specified) throw std::runtime_error ("range of model space not specified");
	if (!incdec_specified) throw std::runtime_error ("inclinations and declinations are not specified");
}

void
Joint::read_data ()
{
	FILE	*fp = fopen (_fn_mag_, "r");
	if (!fp) {
		char	msg[128];
		sprintf (msg, "ERROR: cannot open file %s", _fn_mag_);
		throw std::runtime_error (msg);
	}
	_magdata_ = fread_data_array (fp);
	fclose (fp);

	fp = fopen (_fn_grv_, "r");
	if (!fp) {
		char	msg[128];
		sprintf (msg, "ERROR: cannot open file %s", _fn_grv_);
		throw std::runtime_error (msg);
	}
	_grvdata_ = fread_data_array (fp);
	fclose (fp);
}

void
Joint::set_lower_bounds ()
{
	size_t	m = _nx_ * _ny_ * _nz_;
	_lower_ = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, 2 * m, 1, 2 * m);
	for (size_t i = 0; i < m; i++) {
		_lower_->data[i] = _beta_lower_;
		_lower_->data[m + i] = _rho_lower_;
	}
}

void
Joint::set_surface (size_t c, double *zsurf)
{
	if (_nx_ <= 0 || _ny_ <= 0)
		throw std::runtime_error ("range dose not specified. Call set_range() before.");
	if (_nx_ * _ny_ != c)
		throw std::runtime_error ("dim(terrain) is not match with dim(f) and dim(g)");
	_zsurf_ = new double [c];
	for (size_t i = 0; i < c; i++) _zsurf_[i] = zsurf[i];
}

void
Joint::simeq ()
{
	read_data ();
	if (_nu_ > 0.) set_lower_bounds ();

	if (_fn_ter_ != NULL) {
		size_t	c = count (_fn_ter_);
		double	*zsurf = read_terrain (c, _fn_ter_);
		set_surface (c, zsurf);
		delete zsurf;
	}
	set_mag (_exf_inc_, _exf_dec_, _mgz_inc_, _mgz_dec_, _magdata_);
	set_grv (_grvdata_);

	if (_export_matrix_) export_matrix ();

}

void
Joint::set_mag (double inc, double dec, data_array *data)
{
	_exf_inc_ = inc;
	_exf_dec_ = dec;
	_mgz_inc_ = inc;
	_mgz_dec_ = dec;

	_magdata_ = data;
	size_t	n = _magdata_->n;
	_f_ = mm_real_view_array (MM_REAL_DENSE, MM_REAL_GENERAL, n, 1, n, _magdata_->data);

	_magker_ = new MagKernel (_exf_inc_, _exf_dec_, _mgz_inc_, _mgz_dec_);
	_magker_->set_range (_nx_, _ny_, _nz_, _xrange_, _yrange_, _zrange_, 1000.);
	if (_zsurf_) _magker_->set_surface (_zsurf_);
	_magker_->set_data (_magdata_);	
	_K_ = _magker_->get ();
}

void
Joint::set_mag (double exf_inc, double exf_dec, double mag_inc, double mag_dec, data_array *data)
{
	_exf_inc_ = exf_inc;
	_exf_dec_ = exf_dec;
	_mgz_inc_ = mag_inc;
	_mgz_dec_ = mag_dec;

	_magdata_ = data;
	size_t	n = _magdata_->n;
	_f_ = mm_real_view_array (MM_REAL_DENSE, MM_REAL_GENERAL, n, 1, n, _magdata_->data);

	_magker_ = new MagKernel (_exf_inc_, _exf_dec_, _mgz_inc_, _mgz_dec_);
	_magker_->set_range (_nx_, _ny_, _nz_, _xrange_, _yrange_, _zrange_, 1000.);
	if (_zsurf_) _magker_->set_surface (_zsurf_);
	_magker_->set_data (_magdata_);	
	_K_ = _magker_->get ();
}

void
Joint::set_grv (data_array *data)
{
	_grvdata_ = data;
	size_t	n = _grvdata_->n;
	_g_ = mm_real_view_array (MM_REAL_DENSE, MM_REAL_GENERAL, n, 1, n, _grvdata_->data);

	_grvker_ = new GravKernel ();
	_grvker_->set_range (_nx_, _ny_, _nz_, _xrange_, _yrange_, _zrange_, 1000.);
	if (_zsurf_) _grvker_->set_surface (_zsurf_);
	_grvker_->set_data (_grvdata_);	
	_G_ = _grvker_->get ();
}

void
Joint::fwrite_model (FILE *fp)
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

void
Joint::export_matrix ()
{
	FILE	*fp = fopen ("K.mat", "w");
	if (fp) {
		mm_real_fwrite (fp, _K_, "%.8f");
		fclose (fp);
	}
	fp = fopen ("G.mat", "w");
	if (fp) {
		mm_real_fwrite (fp, _G_, "%.8f");
		fclose (fp);
	}
}

/*** private methods ***/
void
Joint::__init__ ()
{
	strcpy (_fn_settings_, "settings.par");

	_nx_ = -1;
	_ny_ = -1;
	_nz_ = -1;

	_mu_ = 1.;
	_nu_ = -1.;
	_lower_ = NULL;

	_zsurf_ = NULL;

	_tolerance_ = 1.e-3;
	_maxiter_ = 1000;

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

	_export_matrix_ = false;
}


