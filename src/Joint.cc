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

#ifdef __cplusplus
	}
#endif // __cplusplus

#define bufsiz 512

static char *
get_toolname (char *str)
{
	char	*p = strrchr (str, '/');
	if (p == NULL) p = str;
	else p++;
	return p;
}

/*** public methods ***/
// Destructor
Joint::~Joint ()
{

	delete [] fn_ter_;
	delete [] zsurf_;
	if (magdata_) data_array_free (magdata_);
	if (grvdata_) data_array_free (grvdata_);
}

void
Joint::print_usage ()
{
	fprintf (stderr, "=====   PROGRAM   =====\n");
	fprintf (stderr, "%s version %s\n\n", toolname_, _version_info_);

	fprintf (stderr, "===== DESCRIPTION =====\n");
	fprintf (stderr, "This program performs joint inversion of magnetic and gravity data\n");
	fprintf (stderr, "based on L2 norm and group lasso combined regularization.\n\n");	

	fprintf (stderr, "===== USAGE: %s =====\n", toolname_);
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
	toolname_ = get_toolname (argv[0]);

	parse_command_line (argc, argv);

	FILE	*fp = fopen (fn_settings_, "r");
	if (!fp) {
		char	msg[bufsiz];
		sprintf (msg, "cannot open setting file %s", fn_settings_);
		throw std::runtime_error (msg);
	}
	read_settings_file (fp);
	fclose (fp);
}

// start inversion
size_t
Joint::start (bool normalize)
{
	if (f_ == NULL || K_ == NULL || g_ == NULL || G_ == NULL) setup_problem ();

	// scale gravity anomaly
	size_t	ifamax = idamax_ (&size1_mag_, f_, &ione);
	double	famax = fabs (f_[ifamax - 1]);
	size_t	igamax = idamax_ (&size1_grv_, g_, &ione);
	double	gamax = fabs (g_[igamax - 1]);

	scale_ = famax / gamax;
	std::cerr << "scale = " << scale_ << std::endl;
	FILE	*fp = fopen ("scale.data", "w");
	if (fp) {
		fprintf (fp, "%.8e\n", scale_);
		fclose (fp);
	}
	dscal_ (&size1_grv_, &scale_, g_, &ione);

	if (!admm_) {
		admm_ = new mADMM (lambda1_, lambda2_, mu_);
		admm_->setup_problem (size1_mag_, size1_grv_, size2_, f_, g_, K_, G_, normalize, nu_, lower_);
	}
	// export weight for kernel matrix
	export_depth_weights ();
	return admm_->solve (tolerance_, maxiter_, verbos_);
}

// return residual (max of primal and dual residuals)
double
Joint::get_residual ()
{
	if (!admm_) throw std::runtime_error ("admm object has not yet been instantiated. Call start() first.");
	return admm_->get_residual ();
}

// recover magnetic and gravity anomalies
// and store them into f and g
void
Joint::recover_data (double *f, double *g)
{
	admm_->recover_data (f, g);
	double	scale = 1. / scale_;
	dscal_ (&size1_grv_, &scale, g, &ione);
}

// get instance of magnetization model
// depth weighting is removed
double *
Joint::get_magnetization_model ()
{
	double	*beta_ = admm_->get_magnetization ();
	double	*beta = new double [size2_];
	for (size_t j = 0; j < size2_; j++) beta[j] = beta_[j];
	return beta;
}

// get instance of density model
// depth weighting is removed
double *
Joint::get_gravity_model ()
{
	double	*rho_ = admm_->get_density ();
	double	*rho = new double [size2_];
	for (size_t j = 0; j < size2_; j++) rho[j] = rho_[j];
	double	scale = 1. / scale_;
	dscal_ (&size2_, &scale, rho, &ione);
	return rho;
}

// fwrite settings specified by inline options
void
Joint::print_command_line_options (FILE *stream) const
{
	fprintf (stream, "\n");
	fprintf (stream, "input file:\t%s\t%s\n", fn_mag_, fn_grv_);
	if (fn_ter_) fprintf (stream, "terrain file:\t%s\n", fn_ter_); 
	if (alpha_ > 0.) fprintf (stream, "alpha,lambda:\t%.4e,%.4e\n", alpha_, lambda_);
	else fprintf (stream, "lambda1,lambda2:\t%.4e,%.4e\n", lambda1_, lambda2_);
}

// fwrite setting specified by settings file
void
Joint::print_settings (FILE *stream) const
{
	fprintf (stream, "\n");
	fprintf (stream, "number of grid:\t%ld/%ld/%ld\n", nx_, ny_, nz_);
	fprintf (stream, "x,y,z range:\t%.4f/%.4f,%.4f/%.4f,%.4f/%.4f\n", 
			 xrange_[0], xrange_[1], yrange_[0], yrange_[1], zrange_[0], zrange_[1]);
	fprintf (stream, "exf:inc,dec:\t%.4f,%.4f\n", exf_inc_, exf_dec_);
	fprintf (stream, "mag:inc,dec:\t%.4f,%.4f\n", mgz_inc_, mgz_dec_);
	fprintf (stream, "tol,maxiter:\t%.2e,%ld\n", tolerance_, maxiter_);
	fprintf (stream, "mu:\t\t%.4f\n", mu_);
	if (nu_ > DBL_EPSILON)
		fprintf (stream, "nu:\t\t%.4f:lower bounds = %.4f, %.4f\n", nu_, beta_lower_, rho_lower_);
}

// export calculation results
void
Joint::export_results ()
{
	FILE *fp;

	// export derived models
	fp = fopen ("model.data", "w");
	if (!fp) throw std::runtime_error ("cannot open file model.data");

	write_model_to_file (fp);
	fclose (fp);

	double	*fr = new double [size1_mag_];
	double	*gr = new double [size1_grv_];

	// recover input magnetic and gravity anomalies
	recover_data (fr, gr);

	fp = fopen ("recover_mag.data", "w");
	if (!fp) throw std::runtime_error ("cannot open file recover_mag.data");

	fwrite_data_array_with_data (fp, magdata_, fr, "%.4f\t%.4f\t%.4f\t%.4f");
	delete [] fr;
	fclose (fp);

	fp = fopen ("recover_grv.data", "w");
	if (!fp) throw std::runtime_error ("cannot open file recover_grv.data");

	fwrite_data_array_with_data (fp, grvdata_, gr, "%.4f\t%.4f\t%.4f\t%.4f");
	delete [] gr;
	fclose (fp);

	return;
}

/*** protected methods ***/
// read inline options
void
Joint::parse_command_line (int argc, char **argv)
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
				strcpy (fn_mag_, optarg);
				fn_mag_specified = true;
				break;

			case 'g':
				strcpy (fn_grv_, optarg);
				fn_grv_specified = true;
				break;

			case 'l':
				sscanf (optarg, "%lf:%lf", &log10_lambda1, &log10_lambda2);
				if (log10_lambda1 > -16) lambda1_ = pow (10., log10_lambda1);
				else lambda1_ = 0.;
				if (log10_lambda2 > -16) lambda2_ = pow (10., log10_lambda2);
				else lambda2_ = 0.;
				lambda_specified = true;
				break;

			case 'a':
				sscanf (optarg, "%lf:%lf", &alpha_, &log10_lambda);
				if (log10_lambda1 > -16) lambda_ = pow (10., log10_lambda);
				else lambda_ = 0.;
				lambda1_ = alpha_ * lambda_;
				lambda2_ = (1. - alpha_) * lambda_;
				alpha_specified = true;
				break;

			case 's':
				strcpy (fn_settings_, optarg);
				break;

			case 't':
				fn_ter_ = new char [256];
				strcpy (fn_ter_, optarg);
				
				break;

			case 'v':
				verbos_ = true;
				break;

			case 'h':
				print_usage ();
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
Joint::read_settings_file (FILE *stream)
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
				sscanf (p, "%lu, %lu, %lu", &nx_, &ny_, &nz_);
				ngrid_specified = true;
				break;
			case '2': // x,y,z range
				sscanf (p, "%lf,%lf,%lf,%lf,%lf,%lf",
					&xrange_[0], &xrange_[1], &yrange_[0], &yrange_[1], &zrange_[0], &zrange_[1]);
				range_specified = true;
				break;
			case '3': // inc, dec
				sscanf (p, "%lf,%lf,%lf,%lf", &exf_inc_, &exf_dec_, &mgz_inc_, &mgz_dec_);
				incdec_specified = true;
				break;
			case '4': // tolerance, maxiter: default is 1.e-3, 1000
				sscanf (p, "%lf,%lu", &tolerance_, &maxiter_);
				break;
			case '5': // mu: default is 1.0
				mu_ = (double) atof (p);
				break;
			case '6': // lower bounds: default is nu = -1 (lower bounds not applied) 
				sscanf (p, "%lf,%lf,%lf", &nu_, &beta_lower_, &rho_lower_);
				break;
			default:
				break;
		}
	}
	if (!ngrid_specified) throw std::runtime_error ("number of grid not specified");
	if (!range_specified) throw std::runtime_error ("range of model space not specified");
	if (!incdec_specified) throw std::runtime_error ("inclinations and declinations are not specified");
	size2_ = nx_ * ny_ * nz_;
}

// read data files
void
Joint::load_data_files ()
{
	FILE	*fp = fopen (fn_mag_, "r");
	if (!fp) {
		char	msg[bufsiz];
		sprintf (msg, "ERROR: cannot open file %s", fn_mag_);
		throw std::runtime_error (msg);
	}
	magdata_ = fread_data_array (fp);
	fclose (fp);

	fp = fopen (fn_grv_, "r");
	if (!fp) {
		char	msg[bufsiz];
		sprintf (msg, "ERROR: cannot open file %s", fn_grv_);
		throw std::runtime_error (msg);
	}
	grvdata_ = fread_data_array (fp);
	fclose (fp);

	//if (magdata_->n != grvdata_->n) throw std::runtime_error ("numbers of magnetic and gravity data incompatible.");
	size1_mag_ = magdata_->n;
	size1_grv_ = grvdata_->n;

	if (verbos_) {
		std::cerr << "number of mag data: " << size1_mag_ << std::endl;
		std::cerr << "number of grv data: " << size1_grv_ << std::endl;
	}
}

// register lower bound constraint
void
Joint::setup_lower_bounds ()
{
	lower_ = new double [2 * size2_];
	for (size_t j = 0; j < size2_; j++) {
		lower_[j] = beta_lower_;
		lower_[j + size2_] = rho_lower_;
	}
}

// register terrai data
void
Joint::load_surface_topography (size_t c, double *zsurf)
{
	if (nx_ <= 0 || ny_ <= 0)
		throw std::runtime_error ("range dose not specified. Call set_range() before.");
	if (nx_ * ny_ != c)
		throw std::runtime_error ("dim(terrain) is not match with dim(f) and dim(g)");
	zsurf_ = new double [c];
	for (size_t i = 0; i < c; i++) zsurf_[i] = zsurf[i];
}

// set simultaneous equation for joint inversion
void
Joint::setup_problem ()
{
	load_data_files ();
	if (nu_ > DBL_EPSILON) setup_lower_bounds ();

	if (fn_ter_ != NULL) {
		FILE	*fp = fopen (fn_ter_, "r");
		if (!fp) {
			char	msg[bufsiz];
			sprintf (msg, "cannot open terrain file: %s", fn_ter_);
			throw std::runtime_error (msg);
		}
		size_t	c = count_lines_in_file_ (fp);
		double	*zsurf = read_terrain_file (fp, c);
		fclose (fp);
		load_surface_topography (c, zsurf);
		delete zsurf;
	}
	setup_magnetic_problem (exf_inc_, exf_dec_, mgz_inc_, mgz_dec_, magdata_);
	setup_gravity_problem (grvdata_);
}

// set magnetic data, inclination, and declination
void
Joint::setup_magnetic_problem (double exf_inc, double exf_dec, double mgz_inc, double mgz_dec, data_array *data)
{
	f_ = data->data;

	magker_ = new MagKernel (exf_inc, exf_dec, mgz_inc, mgz_dec);
	magker_->set_range (nx_, ny_, nz_, xrange_, yrange_, zrange_, 1000.);
	if (zsurf_) magker_->set_surface (zsurf_);
	magker_->set_data (data);	
	K_ = magker_->get ();
}

// set gravity data
void
Joint::setup_gravity_problem (data_array *data)
{
	g_ = data->data;

	grvker_ = new GravKernel ();
	grvker_->set_range (nx_, ny_, nz_, xrange_, yrange_, zrange_, 1000.);
	if (zsurf_) grvker_->set_surface (zsurf_);
	grvker_->set_data (data);	
	G_ = grvker_->get ();
}

// write derived model into a file
void
Joint::write_model_to_file (FILE *fp)
{
	if (magker_ == NULL) throw std::runtime_error ("magnetic kernel is not specified. Call set_mag()");
	if (grvker_ == NULL) throw std::runtime_error ("gravity kernel is not specified. Call set_grv()");

	double	*beta = get_magnetization_model ();
	double	*rho = get_gravity_model ();

	grid		*grd = magker_->get_grid ();
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
Joint::export_depth_weights ()
{
	FILE	*fp = fopen ("wx.vec", "w");
	if (fp) {
		double	*wx = admm_->get_magnetic_depth_weights ();
		for (size_t j = 0; j < size2_; j++) fprintf (fp, "%.8e\n", wx[j]);
		fclose (fp);
	}
	fp = fopen ("wy.vec", "w");
	if (fp) {
		double	*wy = admm_->get_gravity_depth_weights ();
		for (size_t j = 0; j < size2_; j++) fprintf (fp, "%.8e\n", wy[j]);
		fclose (fp);
	}
}

// read terrain file
double *
Joint::read_terrain_file (FILE *fp, const size_t c)
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

/*** private methods ***/
// initialize
void
Joint::initialize_ ()
{
	strcpy (fn_settings_, "settings.par");
	fn_ter_ = NULL;

	alpha_  = -1.;
	lambda_ = -1.;

	fn_ter_ = NULL;

	nx_ = 0;
	ny_ = 0;
	nz_ = 0;

	zsurf_ = NULL;

	mu_ = 1.;
	nu_ = 0.;
	lower_ = NULL;

	zsurf_ = NULL;

	tolerance_ = 1.e-3;
	maxiter_ = 1000;

	magdata_ = NULL;
	grvdata_ = NULL;

	magker_ = NULL; 
	grvker_ = NULL; 

	f_ = NULL;
	K_ = NULL;

	scale_ = 1.;
	g_ = NULL;
	G_ = NULL;

	tolerance_ = 1.e-3;
	maxiter_ = 1000;

	admm_ = NULL;

	verbos_ = false;
}

// count number of data
size_t
Joint::count_lines_in_file_ (FILE *fp)
{
	int		c = 0;
	char	buf[BUFSIZ];
	while (fgets (buf, BUFSIZ, fp) != NULL) c++;
	fseek (fp, SEEK_SET, 0L);
	return c;
}

