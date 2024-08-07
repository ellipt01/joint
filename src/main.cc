#include <iostream>

#include "mgcal.h"

#include "Kernel.h"
#include "ADMM.h"
#include "mADMM.h"
#include "Joint.h"

int
main (int argc, char **argv)
{
	Joint	joint;

	try {
		joint.prepare (argc, argv);
	} catch (std::runtime_error e) {
		std::cerr << e.what () << std::endl;
		joint.usage ();
		exit (1);
	}

	joint.fwrite_inline (stderr);
	joint.fwrite_settings (stderr);

	try {
		size_t	niter = joint.start (true);
		std::cerr << "number of iterations = " << niter << " / " << joint.get_maxiter () << std::endl;
		joint.export_results ();
	} catch (std::runtime_error e) {
		std::cerr << "ERROR: " << ": " << e.what () << std::endl;
		exit (1);
	}

	return EXIT_SUCCESS;
}

