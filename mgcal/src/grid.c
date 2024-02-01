#include <stdlib.h>
#include <stdbool.h>

#include "../include/vector3d.h"
#include "private/util.h"
#include "grid.h"

static grid *
grid_alloc (void)
{
	grid	*g = (grid *) malloc (sizeof (grid));
	g->n = 0;
	g->nh = 0;

	g->nx = 0;
	g->ny = 0;
	g->nz = 0;

	set_range (g->xrange, 0., 0.);
	set_range (g->yrange, 0., 0.);
	set_range (g->zrange, 0., 0.);

	g->x = NULL;
	g->y = NULL;
	g->z = NULL;
	g->z1 = NULL;

	g->dx = NULL;
	g->dy = NULL;
	g->dz = NULL;

	g->data = NULL;

	return g;
}

static double
grid_array (const int n, const double t0, const double *dt, double *s)
{
	int		i;
	double	t = t0;
	for (i = 0; i < n; i++) {
		s[i] = t + 0.5 * dt[i];
		t += dt[i];
	}
	return t;
}

static bool
grid_set_surface_0 (grid *g, const double *z1)
{
	if (g->z1) free (g->z1);
	g->z1 = (double *) malloc (g->nh * sizeof (double));
	if (!array_copy (g->nh, g->z1, z1)) return false;
	return true;
}


static grid *
grid_new_0 (const int nx, const int ny, const int nz, const double x[], const double y[], const double z[], const double *dx, const double *dy, const double *dz, const double *z1)
{
	double	xx0, xx1, yy0, yy1, zz0, zz1;
	grid	*g = grid_alloc ();

	xx0 = x[0];
	yy0 = y[0];
	zz0 = z[0];

	xx1 = (nx == 1) ? xx0 : x[1];
	yy1 = (ny == 1) ? yy0 : y[1];
	zz1 = (nz == 1) ? zz0 : z[1];

	g->nx = nx;
	g->ny = ny;
	g->nz = nz;
	g->nh = g->nx * g->ny;
	g->n = g->nh * g->nz;

	g->dx = (double *) malloc (nx * sizeof (double));
	if (dx) array_copy (nx, g->dx, dx);
	else {
		double	incx = (xx1 - xx0) / (double) nx;
		array_set_all (nx, g->dx, incx);
	}
	g->x = (double *) malloc (g->nx * sizeof (double));
	xx1 = grid_array (g->nx, xx0, g->dx, g->x);

	g->dy = (double *) malloc (ny * sizeof (double));
	if (dy) array_copy (ny, g->dy, dy);
	else {
		double	incy = (yy1 - yy0) / (double) ny;
		array_set_all (ny, g->dy, incy);
	}
	g->y = (double *) malloc (g->ny * sizeof (double));
	yy1 = grid_array (g->ny, yy0, g->dy, g->y);

	g->dz = (double *) malloc (nz * sizeof (double));
	if (dz) array_copy (nz, g->dz, dz);
	else {
		double	incz = (zz1 - zz0) / (double) nz;
		array_set_all (nz, g->dz, incz);
	}
	g->z = (double *) malloc (g->nz * sizeof (double));
	zz1 = grid_array (g->nz, zz0, g->dz, g->z);

	set_range (g->xrange, xx0, xx1);
	set_range (g->yrange, yy0, yy1);
	set_range (g->zrange, zz0, zz1);

	if (z1) {
		if (!grid_set_surface_0 (g, z1)) return NULL;
	}

	return g;
}

grid *
grid_new (const int nx, const int ny, const int nz, const double x[], const double y[], const double z[])
{
	grid	*g;
	if (nx <= 0 || ny <= 0 || nz <= 0) error_and_exit ("grid_new", "nx, ny, nz must be >= 1.", __FILE__, __LINE__);
	g = grid_new_0 (nx, ny, nz, x, y, z, NULL, NULL, NULL, NULL);
	if (!g) error_and_exit ("grid_new", "failed to create grid object.", __FILE__, __LINE__);
	return g;
}

grid *
grid_new_full (const int nx, const int ny, const int nz, const double x[], const double y[], const double z[], const double *dx, const double *dy, const double *dz, const double *z1)
{
	grid	*g;
	if (nx <= 0 || ny <= 0 || nz <= 0) error_and_exit ("grid_new_full", "nx, ny, nz must be >= 1.", __FILE__, __LINE__);
	g = grid_new_0 (nx, ny, nz, x, y, z, dx, dy, dz, z1);
	if (!g) error_and_exit ("grid_new_full", "failed to create grid object.", __FILE__, __LINE__);
	return g;
}

bool
grid_set_surface (grid *g, const double *z1)
{
	if (!g) error_and_exit ("grid_set_surface", "grid *g is empty.", __FILE__, __LINE__);
	return grid_set_surface_0 (g, z1);
}

void
grid_free (grid *g)
{
	if (g) {
		if (g->x) free (g->x);
		if (g->y) free (g->y);
		if (g->z) free (g->z);
		if (g->z1) free (g->z1);
		if (g->dx) free (g->dx);
		if (g->dy) free (g->dy);
		if (g->dz) free (g->dz);
		free (g);
	}
	return;
}

static bool
grid_check_index (const grid *g, const int i, const int j, const int k, const int h)
{
	if (i < 0 || g->nx <= i) return false;
	if (j < 0 || g->ny <= j) return false;
	if (k < 0 || g->nz <= k) return false;
	if (h < 0 || g->nh <= h) return false;
	return true;
}

static bool
grid_get_index_0 (const grid *g, const int n, int *i, int *j, int *k, int *h)
{
	int		_i, _j, _k, _h;
	if (g->n <= n) return false;

	_k = n / g->nh;
	_h = n % g->nh;
	_j = _h / g->nx;
	_i = _h % g->nx;

	if (i) *i = _i;
	if (j) *j = _j;
	if (k) *k = _k;
	if (h) *h = _h;

	return grid_check_index (g, _i, _j, _k, _h);
}

void
grid_get_index (const grid *g, const int n, int *i, int *j, int *k, int *h)
{
	if (!g) error_and_exit ("grid_get_index", "grid *g is empty.", __FILE__, __LINE__);
	if (!grid_get_index_0 (g, n, i, j, k, h)) error_and_exit ("grid_get_index", "index invalid.", __FILE__, __LINE__);
	return;
}

void
grid_get_nth (const grid *g, const int n, vector3d *center, vector3d *dim)
{
	int		i = 0, j = 0, k = 0, h = 0;
	if (!g) error_and_exit ("grid_get_nth", "grid *g is empty.", __FILE__, __LINE__);
	if (!grid_get_index_0 (g, n, &i, &j, &k, &h)) error_and_exit ("grid_get_nth", "index invalid.", __FILE__, __LINE__);
	if (center) {
		double	zk = g->z[k];
		if (g->z1) zk += g->z1[h];
		vector3d_set (center, g->x[i], g->y[j], zk);
	}
	if (dim) vector3d_set (dim, g->dx[i], g->dy[j], g->dz[k]);
	return;
}

void
grid_stretch_at_edge (grid *g, double l)
{
	if (!g) error_and_exit ("grid_stretch_at_edge", "grid *g is empty.", __FILE__, __LINE__);
	g->x[0] -= l / 2.;
	g->x[g->nx - 1] += l / 2.;
	g->y[0] -= l / 2.;
	g->y[g->ny - 1] += l / 2.;
	g->z[g->nz - 1] -= l / 2.;

	g->dx[0] += l;
	g->dx[g->nx - 1] += l;
	g->dy[0] += l;
	g->dy[g->ny - 1] += l;
	g->dz[g->nz - 1] -= l;
	return;
}

