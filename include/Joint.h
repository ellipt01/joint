#ifndef JOINT_H_
#define JOINT_H_

/***
	@class Joint
	@brief Performs a joint inversion of magnetic and gravity data.

	This class provides a high-level interface for performing joint inversion
	using a combination of L2-norm and group-lasso regularization.

	@section workflow Typical Workflow
	1. prepare():   Parses command-line arguments and reads a settings file.
	2. start():     Executes the joint inversion algorithm.
	3. export_results(): Exports the resulting models and recovered data.
***/

class Joint
{
	char		*toolname_;

	// Path to the settings file
	char		fn_settings_[256];
	// Path to the terrain data file
	char		*fn_ter_;

	// External magnetic field direction (inclination, declination)
	double	exf_inc_;
	double	exf_dec_;

	// Magnetization vector direction (inclination, declination)
	double	mgz_inc_;
	double	mgz_dec_;

	// Input data file paths
	char		fn_mag_[256];
	char		fn_grv_[256];

	// Regularization parameter 
	double	lambda_;
	// Mixing ratio of L1-L2 norm penalty
	double	alpha_;

	// Regularization parameters for L1 and L2 norm penalties
	double	lambda1_;
	double	lambda2_;

	// Number of grid cells in each dimension
	size_t	nx_;
	size_t	ny_;
	size_t	nz_;

	// Number of magnetic/gravity data points
	size_t	size1_mag_;	// = dim(f)
	size_t	size1_grv_;	// = dim(g)
	// Total number of model parameters (grid cells)
	size_t	size2_;	// = size(X, 2) = size(Y, 2)

	// Total number of data points (rows in the system matrix)
	size_t	m_;
	// Total number of variables (2 * model parameters)
	size_t	n_;

	// Physical boundaries of the model space
	double	xrange_[2]; // Eastward positive
	double	yrange_[2]; // Northward positive
	double	zrange_[2]; // Upward positive
	double	*zsurf_; // store terrain data, read from fn_ter_

	// ADMM penalty parameter for the main consensus constraint (s = zeta)
	double	mu_;
	// ADMM penalty parameter for the lower bound constraint
	double	nu_;
	double	beta_lower_; // Lower bound for magnetization
	double	rho_lower_;  // Lower bound for density
	double	*lower_; // Combined lower bounds for the solver

	// Magnetic and gravity observed data
	data_array	*magdata_;
	data_array	*grvdata_;

	// Magnetic data vector and kernel matrix
	double	*f_;
	double	*K_;

	// Gravity data vector and kernel matrix
	double	*g_;
	double	*G_;

	// Scaling factor to balance magnetic and gravity data contributions
	double	scale_; // = |g|_infinity / |f|_infinity

	// Kernel matrix generator objects
	MagKernel	*magker_;
	GravKernel	*grvker_;

	// The ADMM solver instance
	mADMM		*admm_;

	// Convergence tolerance and maximum iterations for the solver
	double	tolerance_;
	size_t	maxiter_;

	// Actual number of iterations performed
	size_t	niter_;

	bool		export_matrix_;
	bool		verbos_;

public:
	Joint () { initialize_ (); }
	~Joint ();

	// Accessor methods for parameters
	double	getAlpha () const { return alpha_; }
	double	getLambda () const { return lambda_; }

	double	getLambda1 () const { return lambda1_; }
	double	getLambda2 () const { return lambda2_; }

	double	getTolerance () const { return tolerance_; }
	size_t	getMaxIterations () const { return maxiter_; }

	// Displays program usage instructions.
	void		printUsage ();

	// Parses command-line arguments and reads the settings file to configure the inversion.
	void		prepare (int argc, char **argv);

	// Executes the joint inversion algorithm.
	size_t	start (bool normalize);

	// Returns the final residual of the ADMM solver.
	double	getResidual ();

	// Recovers the data anomalies using the computed model.
	double	*recoverData (FieldType type);

	// Returns the resulting model vectors from the inversion.
	double	*getModel (ModelType type);

	// Prints the settings parsed from command-line options and the settings file.
	void		printCommandLineOptions (FILE *stream) const;
	void		printSettings (FILE *stream) const;

	// Exports all primary results (models, recovered data) to output files.
	void		exportResults ();

protected:
	// Parses command-line arguments.
	void		parse_command_line (int argc, char **argv);
	// Reads parameters from the settings file.
	void		read_settings_file (FILE *stream);

	// Loads magnetic and gravity data from files.
	void		load_data_files ();
	// Configures the lower bound constraints for the solver.
	void		setup_lower_bounds ();
	// Loads and applies the surface topography data.
	void		load_surface_topography (size_t c, double *zsurf);

	// Sets up the matrices and vectors for the joint inversion problem.
	void		setup_problem ();

	// Configures the magnetic kernel and data.
	void		setup_magnetic_problem (double exf_inc, double exf_dec, double mgz_inc, double mgz_dec, data_array *data);
	// Configures the gravity kernel and data.
	void		setup_gravity_problem (data_array *data);

	// Writes a model vectors to a file stream.
	void		write_model_to_file (FILE *fp);
	// Exports depth weighting function to a file.
	void		export_depth_weights ();
	// Reads terrain data from a file stream.
	double	*read_terrain_file (FILE *fp, const size_t c);
	
private:
	// Initializes all member variables to default values.
	void		initialize_ ();
	// Counts the number of lines in a file to determine data size.
	size_t	count_lines_in_file_ (FILE *fp);
};

#endif
