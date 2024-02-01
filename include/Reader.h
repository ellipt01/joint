class Reader
{

	size_t	_nx_;
	size_t	_ny_;
	size_t	_nz_;

	double	_xrange_[2];
	double	_yrange_[2];
	double	_zrange_[2];

	double	_exf_inc_;
	double	_exf_dec_;

	double	_mgz_inc_;
	double	_mgz_dec_;

	double	_mu_;
	
	bool	_apply_lower_bound_;
	double	_nu_;
	double	_beta0_;
	double	_rho0_;

	double	_tolerance_;
	size_t	_maxiter_;

	bool	_ngrid_specified_;
	bool	_range_specified_;
	bool	_incdec_specified_;
	bool	_invparams_specified_;
	bool	_pparam_specified_;

public:
	Reader ();
	void	fread (FILE *stream);
	void	fwrite (FILE *stream);

	void	ngrid (size_t *nx, size_t *ny, size_t *nz);
	void	range (double x[], double y[], double z[]);
	void	incdec (double *exf_inc, double *exf_dec, double *mgz_inc, double *mgz_dec);
	void	invparams (double *tol, size_t *maxiter);
	void	pparam (double *mu);
	void	lower_bound (double *nu, double *beta0, double *rho0);
	bool	apply_lower_bound () { return _apply_lower_bound_; }

private:
	void	_init_ (void);
};

