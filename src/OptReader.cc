#include <iostream>
#include <cmath>
#include <cstring>
#include <unistd.h>

#include "OptReader.h"

void
OptReader::_init_ ()
{
	_fn_mag_specified_ = false;
	_fn_grv_specified_ = false;

	_fn_ter_ = NULL;

	_lambda_specified_ = false;
	_log10_lambda_ = 0.;

	_alpha_specified_ = false;
	_alpha_ = 0.;

	_output_matrix_ = false;

	_type_ = DATA_TYPE_NONE;

	strcpy (_fn_settings_, "settings.par");
}

OptReader::OptReader (char *toolname)
{
	_toolname_ = toolname;
	_init_ ();
}

void
OptReader::usage ()
{
	fprintf (stderr, "USAGE: %s\n", _toolname_);
	fprintf (stderr, "       -f <magnetic anomaly filename>\n");
	fprintf (stderr, "       -g <gravitic anomaly filename>\n");
	fprintf (stderr, "       -l <log10(lambda)>\n");
	fprintf (stderr, "       -a <alpha>>\n");
	fprintf (stderr, "[optional]\n");
	fprintf (stderr, "       -t <terrain filename>\n");
	fprintf (stderr, "       -p <setting filename:default is settings.par>\n");
	fprintf (stderr, "       -x (output kernel matrices)\n");
	fprintf (stderr, "       -h (show this message)\n");
}

void
OptReader::read (int argc, char **argv)
{
	char	opt;
	while ((opt = getopt (argc, argv, ":f:g:l:a:p:t:xh")) != -1) {
		switch (opt) {
			case 'f':
				strcpy (_fn_mag_, optarg);
				_fn_mag_specified_ = true;
				break;

			case 'g':
				strcpy (_fn_grv_, optarg);
				_fn_grv_specified_ = true;
				break;

			case 'l':
				_log10_lambda_ = (double) atof (optarg);
				_lambda_ = pow (10., _log10_lambda_);
				_lambda_specified_ = true;
				break;

			case 'a':
				_alpha_ = (double) atof (optarg);
				_alpha_specified_ = true;
				break;

			case 'p':
				strcpy (_fn_settings_, optarg);
				break;

			case 't':
				_fn_ter_ = new char [256];
				strcpy (_fn_ter_, optarg);
				break;

			case 'x':
				_output_matrix_ = true;
				break;

			case 'h':
				usage ();

			default:
				break;
		}
	}
	if (!_fn_mag_specified_) throw std::runtime_error ("magnetic anomaly data filename is not specified");
	if (!_fn_grv_specified_) throw std::runtime_error ("gravitic anomaly data filename is not specified");
	if (!_lambda_specified_) throw std::runtime_error ("lambda is not specified");
	if (!_alpha_specified_) throw std::runtime_error ("alpha is not specified");
}

void
OptReader::fwrite (FILE *stream)
{
	fprintf (stream, "\n");
	fprintf (stream, "inut file:\t%s\t%s\n", _fn_mag_, _fn_grv_);
	if (_fn_ter_) fprintf (stream, "terrain file:\t%.s\n", _fn_ter_); 
	fprintf (stream, "lambda:\t%.4f\n", _lambda_);
	fprintf (stream, "alpha :\t%.4f\n", _alpha_);
	fprintf (stream, "output matrix: ");
	(_output_matrix_) ? fprintf (stream, "true\n") : fprintf (stream, "false\n");
}

