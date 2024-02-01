/*
 * kernel.c
 *
 *  Created on: 2015/03/15
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../include/vector3d.h"
#include "source.h"
#include "data_array.h"
#include "grid.h"
#include "kernel.h"
#include "private/util.h"

static mgcal_func *
mgcal_func_alloc (void)
{
	mgcal_func	*f = (mgcal_func *) malloc (sizeof (mgcal_func));
	f->function = NULL;
	f->parameter = NULL;
	return f;
}

mgcal_func *
mgcal_func_new (const mgcal_theoretical func, void *data)
{
	mgcal_func	*f = mgcal_func_alloc ();
	f->function = func;
	f->parameter = data;
	return f;
}

void
mgcal_func_free (mgcal_func *f)
{
	if (f) free (f);
	return;
}

void
kernel_matrix_set (double *a, const data_array *array, const grid *g, const vector3d *mgz, const vector3d *exf, const mgcal_func *f)
{
	int		m;
	int		nx;
	int		ny;
	int		nz;
	int		nh;

	if (!a) error_and_exit ("kernel_matrix_set", "double *a is empty.", __FILE__, __LINE__);

	m = array->n;
	nx = g->nx;
	ny = g->ny;
	nz = g->nz;
	nh = g->nh;

#pragma omp parallel
	{
		size_t		i, j, k, l;
		double		*z1 = NULL;
		vector3d	*obs = vector3d_new (0., 0., 0.);
		source		*src = source_new ();
		if (exf) src->exf = vector3d_copy (exf);
		source_append_item (src);
		src->begin->pos = vector3d_new (0., 0., 0.);
		src->begin->dim = vector3d_new (0., 0., 0.);
		if (mgz) src->begin->mgz = vector3d_copy (mgz);

#pragma omp for
		for (k = 0; k < nz; k++) {
			double			*zk = g->z + k;	// for parallel calculation
			double			*dzk = g->dz + k;	// for parallel calculation
			double			*yj = g->y;
			double			*dyj = g->dy;
			unsigned long	offsetk = ((unsigned long) k) * ((unsigned long) nh) * ((unsigned long) m);
			double			*ak = a + offsetk;
			for (j = 0; j < ny; j++) {
				double			*xi = g->x;
				double			*dxi = g->dx;
				unsigned long	offsetj = ((unsigned long) j) * ((unsigned long) nx) * ((unsigned long) m);
				double			*aj = ak + offsetj;
				if (g->z1) z1 = g->z1 + j * nx;
				for (i = 0; i < nx; i++) {
					double	*al = aj + i * m;
					double	*xl = array->x;
					double	*yl = array->y;
					double	*zl = array->z;
					double	z1k = *zk;
					if (z1) z1k += z1[i];
					vector3d_set (src->begin->pos, *xi, *yj, z1k);
					vector3d_set (src->begin->dim, *dxi, *dyj, *dzk);
					for (l = 0; l < m; l++) {
						vector3d_set (obs, *xl, *yl, *zl);
						*al = f->function (obs, src, f->parameter);
						al++;
						xl++;
						yl++;
						zl++;
					}
					xi++;
					dxi++;
				}
				yj++;
				dyj++;
			}
		}
		vector3d_free (obs);
		source_free (src);
	}
	return;
}

double *
kernel_matrix (const data_array *array, const grid *g, const vector3d *mgz, const vector3d *exf, const mgcal_func *f)
{
	int				m, n;
	double			*a;
	unsigned long	size;

	m = array->n;
	n = g->n;
	size = ((long) m) * ((long) n);
	a = (double *) malloc (size * sizeof (double));
	if (!a) error_and_exit ("kernel_matrix", "failed to allocate memory of *a.", __FILE__, __LINE__);
	kernel_matrix_set (a, array, g, mgz, exf, f);
	return a;
}


void
kernel_matrix_jth_col_set (double *a, size_t stride, size_t j, const data_array *array, const grid *g, const vector3d *mgz, const vector3d *exf, const mgcal_func *f)
{
	int		m;

	if (g->n <= j) error_and_exit ("kernel_matrix_nth_grid", "n exceeds number of grid.", __FILE__, __LINE__);
	m = array->n;

	{
		size_t		i;
		vector3d	*obs = vector3d_new (0., 0., 0.);
		vector3d	*grd = vector3d_new (0., 0., 0.);
		vector3d	*dim = vector3d_new (0., 0., 0.);

		source		*src = source_new ();
		if (exf) src->exf = vector3d_copy (exf);
		source_append_item (src);
		src->begin->pos = vector3d_new (0., 0., 0.);
		src->begin->dim = vector3d_new (0., 0., 0.);
		if (mgz) src->begin->mgz = vector3d_copy (mgz);

		grid_get_nth (g, j, grd, dim);
		vector3d_set (src->begin->pos, grd->x, grd->y, grd->z);
		vector3d_set (src->begin->dim, dim->x, dim->y, dim->z);
		vector3d_free (grd);
		vector3d_free (dim);

		for (i = 0; i < m; i++) {
			vector3d_set (obs, array->x[i], array->y[i], array->z[i]);
			a[i * stride] = f->function (obs, src, f->parameter);
		}
		vector3d_free (obs);
		source_free (src);
	}
	return;
}

double *
kernel_matrix_jth_col_vector (size_t j, const data_array *array, const grid *g, const vector3d *mgz, const vector3d *exf, const mgcal_func *f)
{
	int		m = array->n;
	double	*a;

	if (g->n <= j) error_and_exit ("kernel_matrix_nth_grid", "n exceeds number of grid.", __FILE__, __LINE__);

	a = (double *) malloc (m * sizeof (double));
	kernel_matrix_jth_col_set (a, 1, j, array, g, mgz, exf, f);

	return a;
}

void
kernel_matrix_ith_row_set (double *a, size_t stride, size_t i, const data_array *array, const grid *g, const vector3d *mgz, const vector3d *exf, const mgcal_func *f)
{
	int		n;
	double	xobs, yobs, zobs;

	if (array->n <= i) error_and_exit ("kernel_matrix_mth_site", "m exceeds number of array.", __FILE__, __LINE__);
	n = g->n;

	xobs = array->x[i];
	yobs = array->y[i];
	zobs = array->z[i];
	{
		size_t		j;
		vector3d	*obs = vector3d_new (xobs, yobs, zobs);
		vector3d	*grd = vector3d_new (0., 0., 0.);
		vector3d	*dim = vector3d_new (0., 0., 0.);

		source		*src = source_new ();
		if (exf) src->exf = vector3d_copy (exf);
		source_append_item (src);
		src->begin->pos = vector3d_new (0., 0., 0.);
		src->begin->dim = vector3d_new (0., 0., 0.);
		if (mgz) src->begin->mgz = vector3d_copy (mgz);

		for (j = 0; j < n; j++) {
			grid_get_nth (g, j, grd, dim);
			vector3d_set (src->begin->pos, grd->x, grd->y, grd->z);
			vector3d_set (src->begin->dim, dim->x, dim->y, dim->z);
			a[j * stride] = f->function (obs, src, f->parameter);
		}
		vector3d_free (grd);
		vector3d_free (dim);
		vector3d_free (obs);
		source_free (src);
	}
	return;
}

double *
kernel_matrix_ith_row_vector (size_t i, const data_array *array, const grid *g, const vector3d *mgz, const vector3d *exf, const mgcal_func *f)
{
	int		n = g->n;
	double	*a;

	if (array->n <= i) error_and_exit ("kernel_matrix_mth_site", "m exceeds number of array.", __FILE__, __LINE__);

	a = (double *) malloc (n * sizeof (double));
	kernel_matrix_ith_row_set (a, 1, i, array, g, mgz, exf, f);

	return a;
}

