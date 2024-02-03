/*** jinv main ***/

#include <iostream>
#include <cstring>
#include <unistd.h>

#include "Joint.h"
#include "SettingsReader.h"
#include "util.h"
#include "inv.h"

#include <gsl/gsl_statistics.h>

char	fn_par[256];
char	fn_mag[256];
char	fn_grv[256];
char	*fn_ter = NULL;

bool	output_matrix = false;

static void
usage (void)
{
	fprintf (stderr, "USAGE: %s\n", toolname);
	fprintf (stderr, "       -f <magnetic anomaly filename>\n");
	fprintf (stderr, "       -g <gravitic anomaly filename>\n");
	fprintf (stderr, "       -l <log10(lambda)>\n");
	fprintf (stderr, "       -a <alpha>>\n");
	fprintf (stderr, "[optional]\n");
	fprintf (stderr, "       -t <terrain file>\n");
	fprintf (stderr, "       -p <parameter filename:default is settings.par>\n");
	fprintf (stderr, "       -x (output kernel matrices)\n");
	fprintf (stderr, "       -h (show this message)\n");
	exit (1);
}

static void
read_params (int argc, char **argv)
{

	bool	fn_mag_specified = false;
	bool	fn_grv_specified = false;
	bool	lambda_specified = false;
	bool	alpha_specified = false;
	bool	parfile_specified = false;
	double	log10_lambda = 0.;
	char	opt;
	while ((opt = getopt (argc, argv, ":f:g:l:a:p:t:xh")) != -1) {
		switch (opt) {
			case 'f':
				strcpy (fn_mag, optarg);
				fn_mag_specified = true;
				break;

			case 'g':
				strcpy (fn_grv, optarg);
				fn_grv_specified = true;
				break;

			case 'l':
				log10_lambda = (double) atof (optarg);
				lambda = pow (10., log10_lambda);
				lambda_specified = true;
				break;

			case 'a':
				alpha = (double) atof (optarg);
				alpha_specified = true;
				break;

			case 'p':
				strcpy (fn_par, optarg);
				parfile_specified = true;
				break;

			case 't':
				fn_ter = new char [128];
				strcpy (fn_ter, optarg);
				break;

			case 'x':
				output_matrix = true;
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
	if (!parfile_specified) strcpy (fn_par, fn_params);
}

static void
output_results (Joint *joint, data_array *mag_data, data_array *grv_data)
{
	FILE *fp = fopen ("model.data", "w");
	if (!fp) {
		fprintf (stderr, "cannot open file model.data\n");
		exit (1);
	}
	joint->fwrite (fp);
	fclose (fp);

	mm_real	*fr = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, mag_data->n, 1, mag_data->n);
	mm_real	*gr = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, grv_data->n, 1, grv_data->n);

	joint->recover (fr, gr);

	fp = fopen ("recover_mag.data", "w");
	if (!fp) {
		fprintf (stderr, "cannot open file recover_mag.data\n");
		exit (1);
	}
	fwrite_data_array_with_data (fp, mag_data, fr->data, "%.4f\t%.4f\t%.4f\t%.4f");
	mm_real_free (fr);
	fclose (fp);

	fp = fopen ("recover_grv.data", "w");
	if (!fp) {
		fprintf (stderr, "cannot open file recover_grv.data\n");
		exit (1);
	}
	fwrite_data_array_with_data (fp, grv_data, gr->data, "%.4f\t%.4f\t%.4f\t%.4f");
	mm_real_free (gr);
	fclose (fp);

	return;
}

int
main (int argc, char **argv)
{
	toolname = get_toolname (argv[0]);

	try {
		read_params (argc, argv);
	} catch (std::runtime_error e) {
		std::cerr << "ERROR: " << toolname << ": ";
		std::cerr << e.what () << std::endl;
		usage ();
		exit (1);
	}

	SettingsReader	reader;
	FILE	*fp = fopen (fn_par, "r");
	if (!fp) {
		fprintf (stderr, "ERROR: cannot open file %s\n", fn_par);
		exit (1);
	}
	reader.fread (fp);
	fclose (fp);

	try {
		reader.ngrid (&nx, &ny, &nz);
		reader.range (xrange, yrange, zrange);
		reader.incdec (&exf_inc, &exf_dec, &mgz_inc, &mgz_dec);
		reader.invparams (&tol, &maxiter);
		reader.pparam (&mu);
		if (reader.apply_lower_bound ()) reader.lower_bound (&nu, &beta0, &rho0);
	} catch (std::runtime_error e) {
		std::cerr << "ERROR: " << toolname << ": ";
		std::cerr << e.what () << std::endl;
		exit (1);
	}
	reader.fwrite (stderr);

	fp = fopen (fn_mag, "r");
	if (!fp) {
		fprintf (stderr, "ERROR: cannot open file %s\n", fn_mag);
		exit (1);
	}
	data_array	*mag_data = fread_data_array (fp);
	fclose (fp);

	fp = fopen (fn_grv, "r");
	if (!fp) {
		fprintf (stderr, "ERROR: cannot open file %s\n", fn_grv);
		exit (1);
	}
	data_array	*grv_data = fread_data_array (fp);
	fclose (fp);

	if (nu > 0) {
		size_t	m = nx * ny * nz;
		lower = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, 2 * m, 1, 2 * m);
		for (size_t i = 0; i < m; i++) {
			lower->data[i] = beta0;
			lower->data[m + i] = rho0;
		}
	}

	Joint	joint (alpha, lambda, mu, nu, lower);
	try {
		joint.set_range (nx, ny, nz, xrange, yrange, zrange);
		if (fn_ter) {
			size_t	c = count (fn_ter);
			double	*zsurf = read_terrain (c, fn_ter);
			joint.set_surface (c, zsurf);
			delete zsurf;
		}
		joint.set_mag (exf_inc, exf_dec, mgz_inc, mgz_dec, mag_data);
		joint.set_grv (grv_data);
	} catch (std::runtime_error e) {
		std::cerr << "ERROR: " << toolname << ": " << e.what () << std::endl;
		exit (1);
	}

	if (output_matrix) {
		fp = fopen ("K.mat", "w");
		if (fp) {
			mm_real_fwrite (fp, joint.get_K (), "%.8f");
			fclose (fp);
		}
		fp = fopen ("G.mat", "w");
		if (fp) {
			mm_real_fwrite (fp, joint.get_G (), "%.8f");
			fclose (fp);
		}
	}

	try {
		size_t	niter = joint.start (tol, maxiter, true);
		std::cerr << "number of iterations = " << niter << " / " << maxiter << std::endl;
	} catch (std::runtime_error e) {
		std::cerr << "ERROR: " << toolname << ": " << e.what () << std::endl;
		exit (1);
	}

	output_results (&joint, mag_data, grv_data);

	return EXIT_SUCCESS;
}
