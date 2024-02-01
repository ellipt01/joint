#ifndef _INV_H_
#define _INV_H_

#ifndef _DATA_TYPE_ENUM_
#define _DATA_TYPE_ENUM_
typedef enum {
	DATA_TYPE_MAGNETIC = 0,
	DATA_TYPE_GRAVITY = 1,
	DATA_TYPE_NONE = -1
} data_type;
#endif

static data_type	type = DATA_TYPE_NONE;

static char			*toolname;

static size_t		nx;
static size_t		ny;
static size_t		nz;

static double		xrange[2];
static double		yrange[2];
static double		zrange[2];

static double		exf_inc;
static double		exf_dec;
static double		mgz_inc;
static double		mgz_dec;

static double		alpha;
static double		log10_lambda;
static double		lambda;
static double		mu = 1.;

static double		nu = -1.;
static double		beta0;
static double		rho0;
static mm_real		*lower = NULL;


static double		tol;
static size_t		maxiter;

const char			fn_params[] = "settings.par";

#endif
