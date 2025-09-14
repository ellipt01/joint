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
	} catch (const std::exception &e) {
		std::cerr << e.what () << std::endl;
		joint.printUsage ();
		exit (1);
	}

	joint.printCommandLineOptions (stderr);
	joint.printSettings (stderr);

	try {
		size_t	niter = joint.start (true);
		std::cerr << "number of iterations = " << niter << " / " << joint.getMaxIterations () << std::endl;
		joint.exportResults ();
	} catch (const std::exception &e) {
		std::cerr << "ERROR: " << ": " << e.what () << std::endl;
		exit (1);
	}

	return EXIT_SUCCESS;
}

