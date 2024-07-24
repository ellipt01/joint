#include <iostream>
#include <cstring>
#include <unistd.h>
#include <cfloat>

#include "mgcal.h"

#include "Kernel.h"
#include "ADMM.h"
#include "mADMM.h"
#include "Joint.h"

#include "_version_.h"

#ifdef __cplusplus
	extern "C" {
#endif // __cplusplus

static size_t	ione = 1;
static double	dzero = 0.;
static double	done = 1.;
static double	dmone = -1.;

size_t	idamax_ (size_t *n, double *x, size_t *inc);
void	dscal_ (size_t *n, double *scale, double *x, size_t *inc);

static char *
get_toolname (char *str)
{
	char	*p = strrchr (str, '/');
	if (p == NULL) p = str;
	else p++;
	return p;
}

/*** public methods ***/
void
Joint::usage ()
{
	fprintf (stderr, "=====   PROGRAM   =====\n");
	fprintf (stderr, "%s version %s\n\n", _toolname_, _version_info_);

	fprintf (stderr, "===== DESCRIPTION =====\n");
	fprintf (stderr, "This program performs joint inversion of magnetic and gravity data\n");
	fprintf (stderr, "based on L2 norm and group lasso combined regularization.\n\n");	

	fprintf (stderr, "===== USAGE: %s =====\n", _toolname_);
	fprintf (stderr, "       -f <magnetic anomaly filename>\n");
	fprintf (stderr, "       -g <gravitic anomaly filename>\n");
	fprintf (stderr, "       -l <log10(lambda1):log10(lambda2)>\n");
	fprintf (stderr, "       -a <alpha:log10(lambda)>\n");
	fprintf (stderr, "[optional]\n");
	fprintf (stderr, "       -t <terrain filename>\n");
	fprintf (stderr, "       -s <setting filename:default is settings.par>\n");
	fprintf (stderr, "       -v (verbos mode)\n");
	fprintf (stderr, "       -h (show this message)\n");
	fprintf (stderr, "[notice]\n");
	fprintf (stderr, "       If log10(lambda) <= -16, lambda is set to 0.\n");
	fprintf (stderr, "       The -l and -a options are exclusive and cannot be used simultaneously.\n\n");
}

// read inline options and settings file
void
Joint::prepare (int argc, char **argv)
{
	_toolname_ = get_toolname (argv[0]);

	_read_inline_ (argc, argv);

	FILE	*fp = fopen (_fn_settings_, "r");
	if (!fp) {
		char	msg[128];
		sprintf (msg, "cannot open setting file %s", _fn_settings_);
		throw std::runtime_error (msg);
	}
	_fread_settings_ (fp);
	fclose (fp);
}

// start inversion
size_t
Joint::start (bool normalize)
{
	if (_f_ == NULL || _K_ == NULL || _g_ == NULL || _G_ == NULL) _simeq_ ();

	// scale gravity anomaly
	size_t	ifamax = idamax_ (&_size1_, _f_, &ione);
	double	famax = fabs (_f_[ifamax - 1]);
	size_t	igamax = idamax_ (&_size1_, _g_, &ione);
	double	gamax = fabs (_g_[igamax - 1]);

	_scale_ = famax / gamax;
	std::cerr << "scale = " << _scale_ << std::endl;
	FILE	*fp = fopen ("scale.data", "w");
	if (fp) {
		fprintf (fp, "%.8e\n", _scale_);
		fclose (fp);
	}
	dscal_ (&_size1_, &_scale_, _g_, &ione);

	if (!_admm_) {
		_admm_ = new mADMM (_lambda1_, _lambda2_, _mu_);
		_admm_->simeq (_size1_, _size2_, _f_, _g_, _K_, _G_, normalize, _nu_, _lower_);
	}
	// export weight for kernel matrix
	_export_weights_ ();
	return _admm_->start (_tolerance_, _maxiter_, _verbos_);
}

// return residual (max of primal and dual residuals)
double
Joint::residual ()
{
	if (!_admm_) throw std::runtime_error ("admm object has not yet been instantiated. Call start() first.");
	return _admm_->residual ();
}

// recover magnetic and gravity anomalies
// and store them into f and g
void
Joint::recover (double *f, double *g)
{
	_admm_->recover (f, g);
	double	scale = 1. / _scale_;
	dscal_ (&_size1_, &scale, g, &ione);
}

// get instance of magnetization model
// depth weighting is removed
double *
Joint::get_beta ()
{
	double	*_beta_ = _admm_->get_beta ();
	double	*beta = new double [_size2_];
	for (size_t j = 0; j < _size2_; j++) beta[j] = _beta_[j];
	return beta;
}

// get instance of density model
// depth weighting is removed
double *
Joint::get_rho ()
{
	double	*_rho_ = _admm_->get_rho ();
	double	*rho = new double [_size2_];
	for (size_t j = 0; j < _size2_; j++) rho[j] = _rho_[j];
	double	scale = 1. / _scale_;
	dscal_ (&_size2_, &scale, rho, &ione);
	return rho;
}

// fwrite settings specified by inline options
void
Joint::fwrite_inline (FILE *stream)
{
	fprintf (stream, "\n");
	fprintf (stream, "input file:\t%s\t%s\n", _fn_mag_, _fn_grv_);
	if (_fn_ter_) fprintf (stream, "terrain file:\t%s\n", _fn_ter_); 
	if (_alpha_ > 0.) fprintf (stream, "alpha,lambda:\t%.4e,%.4e\n", _alpha_, _lambda_);
	else fprintf (stream, "lambda1,lambda2:\t%.4e,%.4e\n", _lambda1_, _lambda2_);
}

// fwrite setting specified by settings file
void
Joint::fwrite_settings (FILE *stream)
{
	fprintf (stream, "\n");
	fprintf (stream, "number of grid:\t%ld/%ld/%ld\n", _nx_, _ny_, _nz_);
	fprintf (stream, "x,y,z range:\t%.4f/%.4f,%.4f/%.4f,%.4f/%.4f\n", 
			 _xrange_[0], _xrange_[1], _yrange_[0], _yrange_[1], _zrange_[0], _zrange_[1]);
	fprintf (stream, "exf:inc,dec:\t%.4f,%.4f\n", _exf_inc_, _exf_dec_);
	fprintf (stream, "mag:inc,dec:\t%.4f,%.4f\n", _mgz_inc_, _mgz_dec_);
	fprintf (stream, "tol,maxiter:\t%.2e,%ld\n", _tolerance_, _maxiter_);
	fprintf (stream, "mu:\t\t%.4f\n", _mu_);
	if (_nu_ > DBL_EPSILON)
		fprintf (stream, "nu:\t\t%.4f:lower bounds = %.4f, %.4f\n", _nu_, _beta_lower_, _rho_lower_);
}

// export calculation results
void
Joint::export_results ()
{
	FILE *fp;

	// export derived models
	fp = fopen ("model.data", "w");
	if (!fp) throw std::runtime_error ("cannot open file model.data");

	_fwrite_model_ (fp);
	fclose (fp);

	double	*fr = new double [_size1_];
	double	*gr = new double [_size1_];

	// recover input magnetic and gravity anomalies
	recover (fr, gr);

	fp = fopen ("recover_mag.data", "w");
	if (!fp) throw std::runtime_error ("cannot open file recover_mag.data");

	fwrite_data_array_with_data (fp, _magdata_, fr, "%.4f\t%.4f\t%.4f\t%.4f");
	delete [] fr;
	fclose (fp);

	fp = fopen ("recover_grv.data", "w");
	if (!fp) throw std::runtime_error ("cannot open file recover_grv.data");

	fwrite_data_array_with_data (fp, _grvdata_, gr, "%.4f\t%.4f\t%.4f\t%.4f");
	delete [] gr;
	fclose (fp);

	return;
}

/*** protected methods ***/
// read inline options
void
Joint::_read_inline_ (int argc, char **argv)
{
	bool	fn_mag_specified = false;
	bool	fn_grv_specified = false;
	bool	lambda_specified = false;
	bool	alpha_specified = false;

	double	log10_lambda, log10_lambda1, log10_lambda2;

	char	opt;
	while ((opt = getopt (argc, argv, ":f:g:l:a:s:t:vh")) != -1) {
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
				sscanf (optarg, "%lf:%lf", &log10_lambda1, &log10_lambda2);
				if (log10_lambda1 > -16) _lambda1_ = pow (10., log10_lambda1);
				else _lambda1_ = 0.;
				if (log10_lambda2 > -16) _lambda2_ = pow (10., log10_lambda2);
				else _lambda2_ = 0.;
				lambda_specified = true;
				break;

			case 'a':
				sscanf (optarg, "%lf:%lf", &_alpha_, &log10_lambda);
				if (log10_lambda1 > -16) _lambda_ = pow (10., log10_lambda);
				else _lambda_ = 0.;
				_lambda1_ = _alpha_ * _lambda_;
				_lambda2_ = (1. - _alpha_) * _lambda_;
				alpha_specified = true;
				break;

			case 's':
				strcpy (_fn_settings_, optarg);
				break;

			case 't':
				_fn_ter_ = new char [256];
				strcpy (_fn_ter_, optarg);
				
				break;

			case 'v':
				_verbos_ = true;
				break;

			case 'h':
				usage ();
				exit (1);

			case '?':
			default:
				throw std::runtime_error ("unknown option is specified");
				break;
		}
	}
	if (!fn_mag_specified) throw std::runtime_error ("magnetic anomaly data filename is not specified");
	if (!fn_grv_specified) throw std::runtime_error ("gravitic anomaly data filename is not specified");
	if (!lambda_specified && !alpha_specified)
		throw std::runtime_error ("lambda1:lambda2 or alpha:lambda must be specified");
	if (lambda_specified && alpha_specified)
		throw std::runtime_error ("-l and -a options cannot be used simultaneously.");
}

// read settings file
void
Joint::_fread_settings_ (FILE *stream)
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
	_size2_ = _nx_ * _ny_ * _nz_;
}

// read data files
void
Joint::_read_data_ ()
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

	if (_magdata_->n != _grvdata_->n) throw std::runtime_error ("numbers of magnetic and gravity data incompatible.");
	_size1_ = _magdata_->n;
}

// register lower bound constraint
void
Joint::_set_lower_bounds_ ()
{
	_lower_ = new double [2 * _size2_];
	for (size_t j = 0; j < _size2_; j++) {
		_lower_[j] = _beta_lower_;
		_lower_[j + _size2_] = _rho_lower_;
	}
}

// register terrai data
void
Joint::_set_surface_ (size_t c, double *zsurf)
{
	if (_nx_ <= 0 || _ny_ <= 0)
		throw std::runtime_error ("range dose not specified. Call set_range() before.");
	if (_nx_ * _ny_ != c)
		throw std::runtime_error ("dim(terrain) is not match with dim(f) and dim(g)");
	_zsurf_ = new double [c];
	for (size_t i = 0; i < c; i++) _zsurf_[i] = zsurf[i];
}

// set simultaneous equation for joint inversion
void
Joint::_simeq_ ()
{
	_read_data_ ();
	if (_nu_ > DBL_EPSILON) _set_lower_bounds_ ();

	if (_fn_ter_ != NULL) {
		FILE	*fp = fopen (_fn_ter_, "r");
		if (!fp) {
			char	msg[80];
			sprintf (msg, "cannot open terrain file: %s", _fn_ter_);
			throw std::runtime_error (msg);
		}
		size_t	c = __count__ (fp);
		double	*zsurf = __read_terrain__ (fp, c);
		fclose (fp);
		_set_surface_ (c, zsurf);
		delete zsurf;
	}
	_set_mag_ (_exf_inc_, _exf_dec_, _mgz_inc_, _mgz_dec_, _magdata_);
	_set_grv_ (_grvdata_);
}

// set magnetic data, inclination, and declination
void
Joint::_set_mag_ (double exf_inc, double exf_dec, double mgz_inc, double mgz_dec, data_array *data)
{
	_f_ = data->data;

	_magker_ = new MagKernel (exf_inc, exf_dec, mgz_inc, mgz_dec);
	_magker_->set_range (_nx_, _ny_, _nz_, _xrange_, _yrange_, _zrange_, 1000.);
	if (_zsurf_) _magker_->set_surface (_zsurf_);
	_magker_->set_data (data);	
	_K_ = _magker_->get ();
}

// set gravity data
void
Joint::_set_grv_ (data_array *data)
{
	_g_ = data->data;

	_grvker_ = new GravKernel ();
	_grvker_->set_range (_nx_, _ny_, _nz_, _xrange_, _yrange_, _zrange_, 1000.);
	if (_zsurf_) _grvker_->set_surface (_zsurf_);
	_grvker_->set_data (data);	
	_G_ = _grvker_->get ();
}

// write derived model into a file
void
Joint::_fwrite_model_ (FILE *fp)
{
	if (_magker_ == NULL) throw std::runtime_error ("magnetic kernel is not specified. Call set_mag()");
	if (_grvker_ == NULL) throw std::runtime_error ("gravity kernel is not specified. Call set_grv()");

	double	*beta = get_beta ();
	double	*rho = get_rho ();

	grid		*grd = _magker_->get_grid ();
	vector3d	*pos = vector3d_new (0., 0., 0.);
	for (size_t k = 0; k < grd->n; k++) {
		grid_get_nth (grd, k, pos, NULL);
		fprintf (fp, "%.4e\t%.4e\t%.4e\t%.8e\t%.8e\n", pos->x, pos->y, pos->z, beta[k], rho[k]);
	}

	delete [] beta;
	delete [] rho;

#ifdef DEBUG
	FILE	*fp_grd = fopen ("grid.data", "w");
	if (fp) {
		fwrite_grid (fp_grd, grd);
		fclose (fp_grd);
	}
#endif

	delete [] pos;
}

// export depth weightings
void
Joint::_export_weights_ ()
{
	FILE	*fp = fopen ("wx.vec", "w");
	if (fp) {
		double	*wx = _admm_->get_wx ();
		for (size_t j = 0; j < _size2_; j++) fprintf (fp, "%.8e\n", wx[j]);
		fclose (fp);
	}
	fp = fopen ("wy.vec", "w");
	if (fp) {
		double	*wy = _admm_->get_wy ();
		for (size_t j = 0; j < _size2_; j++) fprintf (fp, "%.8e\n", wy[j]);
		fclose (fp);
	}
}

/*** private methods ***/
// initialize
void
Joint::__init__ ()
{
	strcpy (_fn_settings_, "settings.par");
	_fn_ter_ = NULL;

	_alpha_  = -1.;
	_lambda_ = -1.;

	_fn_ter_ = NULL;

	_nx_ = 0;
	_ny_ = 0;
	_nz_ = 0;

	_zsurf_ = NULL;

	_mu_ = 1.;
	_nu_ = 0.;
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

	_tolerance_ = 1.e-3;
	_maxiter_ = 1000;

	_admm_ = NULL;

	_verbos_ = false;
}

// count number of data
size_t
Joint::__count__ (FILE *fp)
{
	int		c = 0;
	char	buf[BUFSIZ];
	while (fgets (buf, BUFSIZ, fp) != NULL) c++;
	fseek (fp, SEEK_SET, 0L);
	return c;
}

// read terrain file
double *
Joint::__read_terrain__ (FILE *fp, const size_t c)
{
	double	*zsurf = new double [c];
	char	buf[BUFSIZ];
	int		k = 0;
	while (fgets (buf, BUFSIZ, fp) != NULL) {
		double	x, y, z;
		sscanf (buf, "%lf\t%lf\t%lf", &x, &y, &z);
		zsurf[k] = z;
		if (++k >= c) break;
	}
	return zsurf;
}

#ifdef __cplusplus
	}
#endif // __cplusplus

