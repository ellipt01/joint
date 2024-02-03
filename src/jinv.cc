/*** jinv main ***/

#include <iostream>

#include "jinv.h"

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
	OptReader	opt (get_toolname (argv[0]));

	try {
		opt.read (argc, argv);
	} catch (std::runtime_error e) {
		std::cerr << "ERROR: " << opt.toolname () << ": ";
		std::cerr << e.what () << std::endl;
		opt.usage ();
		exit (1);
	}
	opt.fwrite (stderr);

	SettingsReader	reader;
	FILE	*fp = fopen (opt.fn_settings (), "r");
	if (!fp) {
		fprintf (stderr, "ERROR: cannot open file %s\n", opt.fn_settings ());
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
		std::cerr << "ERROR: " << opt.toolname () << ": ";
		std::cerr << e.what () << std::endl;
		exit (1);
	}
	reader.fwrite (stderr);

	fp = fopen (opt.fn_mag (), "r");
	if (!fp) {
		fprintf (stderr, "ERROR: cannot open file %s\n", opt.fn_mag ());
		exit (1);
	}
	data_array	*mag_data = fread_data_array (fp);
	fclose (fp);

	fp = fopen (opt.fn_grv (), "r");
	if (!fp) {
		fprintf (stderr, "ERROR: cannot open file %s\n", opt.fn_grv ());
		exit (1);
	}
	data_array	*grv_data = fread_data_array (fp);
	fclose (fp);

	if (nu > 0.) {
		size_t	m = nx * ny * nz;
		lower = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, 2 * m, 1, 2 * m);
		for (size_t i = 0; i < m; i++) {
			lower->data[i] = beta0;
			lower->data[m + i] = rho0;
		}
	}

	Joint	joint (opt.alpha (), opt.lambda (), mu, nu, lower);
	try {
		joint.set_range (nx, ny, nz, xrange, yrange, zrange);
		if (opt.fn_ter ()) {
			size_t	c = count (opt.fn_ter ());
			double	*zsurf = read_terrain (c, opt.fn_ter ());
			joint.set_surface (c, zsurf);
			delete zsurf;
		}
		joint.set_mag (exf_inc, exf_dec, mgz_inc, mgz_dec, mag_data);
		joint.set_grv (grv_data);
	} catch (std::runtime_error e) {
		std::cerr << "ERROR: " << opt.toolname () << ": " << e.what () << std::endl;
		exit (1);
	}

	if (opt.output_matrix ()) {
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
		std::cerr << "ERROR: " << opt.toolname () << ": " << e.what () << std::endl;
		exit (1);
	}

	output_results (&joint, mag_data, grv_data);

	return EXIT_SUCCESS;
}
