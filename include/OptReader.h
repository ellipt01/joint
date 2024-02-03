
typedef enum {
	DATA_TYPE_MAGNETIC = 0,
	DATA_TYPE_GRAVITY = 1,
	DATA_TYPE_NONE = -1
} data_type;


class OptReader
{
	char	_fn_mag_[256];
	char	_fn_grv_[256];
	char	*_fn_ter_;

	double	_log10_lambda_;
	double	_lambda_;
	double	_alpha_;

	char	_fn_settings_[256];

	bool	_output_matrix_;

	data_type	_type_;

	char	*_toolname_;

	bool	_fn_mag_specified_;
	bool	_fn_grv_specified_;
	bool	_lambda_specified_;
	bool	_alpha_specified_;

public:
	OptReader (char *toolname);

	void	usage ();
	void	read (int argc, char **argv);

	char	*toolname () { return _toolname_; }

	char	*fn_mag () { return _fn_mag_; }
	char	*fn_grv () { return _fn_grv_; }
	char	*fn_ter () { return _fn_ter_; }
	char	*fn_settings () { return _fn_settings_; }

	double	lambda () { return _lambda_; }
	double	alpha () { return _alpha_; }
	
	bool	output_matrix () { return _output_matrix_; }

	void	fwrite (FILE *stream);

private:
	void	_init_ ();
};
