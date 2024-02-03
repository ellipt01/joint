#include "Joint.h"
#include "OptReader.h"
#include "SettingsReader.h"
#include "util.h"

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

static double		mu = 1.;

static double		nu = -1.;
static double		beta0;
static double		rho0;
static mm_real		*lower = NULL;

static double		tol;
static size_t		maxiter;

