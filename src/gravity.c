#include <stdio.h>
#include <math.h>
#include <float.h>

#include <mgcal.h>

#include "gravity.h"

#define SIGN(a) ((a) < 0. ? -1. : +1.)

const double	G = 6.6743015;

void
source_set_density (source *s, double rho)
{
	source_set_magnetization (s, rho, 0., 0.);
}

static double
masspoint_kernel (const double x, const double y, const double z)
{
	double	r = sqrt (pow (x, 2.) + pow (y, 2.) + pow (z, 2.));
	return - G * z / pow (r, 3.);
}

double
grav_masspoint (const vector3d *obs, const source *s, void *data)
{
	double		x, y, z;
	double		x0, y0, z0;
	double		tmp;
	source_item	*cur;

	double		g = 0.;
#ifdef IN_MGCAL
	if (!obs) error_and_exit ("grav_masspoint", "vector3d *obs is empty.", __FILE__, __LINE__);
	if (!s) error_and_exit ("grav_masspoint", "source *s is empty.", __FILE__, __LINE__);
#endif
	x0 = obs->x;
	y0 = obs->y;
	z0 = obs->z;

	cur = s->begin;
	while (cur) {
		double	dv = 1.0;
#ifdef IN_MGCAL
		if (!cur->pos) error_and_exit ("grav_masspoint", "position of source item is empty.", __FILE__, __LINE__);
#endif
		if (cur->dim) {
			double	dx = cur->dim->x;
			double	dy = cur->dim->y;
			double	dz = cur->dim->z;
			double	dd = dx * dy;
			if (fabs (dz) > __DBL_EPSILON__) dd *= dz;
			dv = fabs (dd);
		}
		x = cur->pos->x;
		y = cur->pos->y;
		z = cur->pos->z;
		tmp = masspoint_kernel (x - x0, y - y0, z - z0);
		g += dv * tmp;
		cur = cur->next;
	}
	return g;
}

static double
prism_kernel (const double x, const double y, const double z)
{
	double	r = sqrt (pow (x, 2.) + pow (y, 2.) + pow (z, 2.));
	double	g = x * log (r + y) + y * log (r + x) - z * atan (x * y / (z * r));
	return G * g;
}

double
grav_prism (const vector3d *obs, const source *s, void *data)
{
	double		a[2], b[2], c[2];
	double		x, y, z;
	double		x0, y0, z0;
	double		tmp[8];
	source_item	*cur;
	
	double		g = 0.;

#ifdef IN_MGCAL
	if (!obs) error_and_exit ("grav_masspoint", "vector3d *obs is empty.", __FILE__, __LINE__);
	if (!s) error_and_exit ("grav_masspoint", "source *s is empty.", __FILE__, __LINE__);
#endif
	x0 = obs->x;
	y0 = obs->y;
	z0 = obs->z;

	cur = s->begin;
	while (cur) {
		double	dx, dy, dz;
		double	flag;
#ifdef IN_MGCAL
		if (!cur->pos) error_and_exit ("grav_masspoint", "position of source item is empty.", __FILE__, __LINE__);
#endif
		dx = cur->dim->x;
		dy = cur->dim->y;
		dz = cur->dim->z;
		flag = SIGN (dx) * SIGN (dy) * SIGN (dz);

		x = cur->pos->x;
		y = cur->pos->y;
		z = cur->pos->z;

		a[0] = x - 0.5 * dx - x0;
		b[0] = y - 0.5 * dy - y0;
		c[0] = z - 0.5 * dz - z0;

		a[1] = a[0] + dx;
		b[1] = b[0] + dy;
		c[1] = c[0] + dz;

		tmp[0] = prism_kernel (a[1], b[1], c[1]);
		tmp[2] = prism_kernel (a[1], b[0], c[1]);
		tmp[4] = prism_kernel (a[0], b[1], c[1]);
		tmp[6] = prism_kernel (a[0], b[0], c[1]);

		if (fabs (dz) < __DBL_EPSILON__) {
			tmp[1] = 0.;
			tmp[3] = 0.;
			tmp[5] = 0.;
			tmp[7] = 0.;
		} else {
			tmp[1] = prism_kernel (a[1], b[1], c[0]);
			tmp[3] = prism_kernel (a[1], b[0], c[0]);
			tmp[5] = prism_kernel (a[0], b[1], c[0]);
			tmp[7] = prism_kernel (a[0], b[0], c[0]);
		}

		g += flag * (tmp[0] - tmp[1] - tmp[2] + tmp[3]
			- tmp[4] + tmp[5] + tmp[6] - tmp[7]);
		cur = cur->next;
	}
	return g;
}


static double
masspoint_yz_kernel (const double y, const double z)
{
	double	r2 = pow (y, 2.) + pow (z, 2.);
	return - G * 2. * z / r2;
}

double
grav_masspoint_yz (const vector3d *obs, const source *s, void *data)
{
	double		y, z;
	double		y0, z0;
	double		tmp;
	source_item	*cur;

	double		g = 0.;
#ifdef IN_MGCAL
	if (!obs) error_and_exit ("grav_masspoint", "vector3d *obs is empty.", __FILE__, __LINE__);
	if (!s) error_and_exit ("grav_masspoint", "source *s is empty.", __FILE__, __LINE__);
#endif
	y0 = obs->y;
	z0 = obs->z;

	cur = s->begin;
	while (cur) {
		double	dv = 1.0;
#ifdef IN_MGCAL
		if (!cur->pos) error_and_exit ("grav_masspoint", "position of source item is empty.", __FILE__, __LINE__);
#endif
		if (cur->dim) {
			double	dy = cur->dim->y;
			double	dz = cur->dim->z;
			double	dd = dy;
			if (fabs (dz) > __DBL_EPSILON__) dd *= dz;
			dv = fabs (dd);
		}
		y = cur->pos->y;
		z = cur->pos->z;
		tmp = masspoint_yz_kernel (y - y0, z - z0);
		g += dv * tmp;
		cur = cur->next;
	}
	return g;
}

static double
prism_yz_kernel (const double y, const double z)
{
	double	r2 = pow (y, 2.) + pow (z, 2.);
	double	g = - y * log (r2) - 2. * z * atan (y / z);
	return G * g;
}

double
grav_prism_yz (const vector3d *obs, const source *s, void *data)
{
	double		b[2], c[2];
	double		y, z;
	double		y0, z0;
	double		tmp[8];
	source_item	*cur;
	
	double		g = 0.;

#ifdef IN_MGCAL
	if (!obs) error_and_exit ("grav_masspoint", "vector3d *obs is empty.", __FILE__, __LINE__);
	if (!s) error_and_exit ("grav_masspoint", "source *s is empty.", __FILE__, __LINE__);
#endif
	y0 = obs->y;
	z0 = obs->z;

	cur = s->begin;
	while (cur) {
		double	dy, dz;
		double	flag;
#ifdef IN_MGCAL
		if (!cur->pos) error_and_exit ("grav_masspoint", "position of source item is empty.", __FILE__, __LINE__);
#endif
		dy = cur->dim->y;
		dz = cur->dim->z;
		flag = SIGN (dy) * SIGN (dz);

		y = cur->pos->y;
		z = cur->pos->z;

		b[0] = y - 0.5 * dy - y0;
		c[0] = z - 0.5 * dz - z0;

		b[1] = b[0] + dy;
		c[1] = c[0] + dz;

		tmp[0] = prism_yz_kernel (b[1], c[1]);
		tmp[1] = prism_yz_kernel (b[0], c[1]);

		if (fabs (dz) < __DBL_EPSILON__) {
			tmp[2] = 0.;
			tmp[3] = 0.;
		} else {
			tmp[2] = prism_yz_kernel (b[1], c[0]);
			tmp[3] = prism_yz_kernel (b[0], c[0]);
		}

		g += flag * (tmp[0] - tmp[1] - tmp[2] + tmp[3]);
		cur = cur->next;
	}
	return g;
}

