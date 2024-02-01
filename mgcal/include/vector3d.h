/*
 * vector3d.h
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#ifndef VECTOR3D_H_
#define VECTOR3D_H_

#ifdef __cplusplus
extern "C" {
#endif

#define	deg2rad(d)	((d) * M_PI / 180.)
#define	rad2deg(r)	((r) * 180. / M_PI)

/* cartesian vector */
typedef struct s_vector3d	vector3d;
/* poler vector */
typedef struct s_pvector	pvector;

struct s_vector3d {
	double	x;
	double	y;
	double	z;
};

typedef enum {
	XCOMP,
	YCOMP,
	ZCOMP
} CCOMPONENT;

vector3d	*vector3d_new (const double x, const double y, const double z);
vector3d	*vector3d_new_with_geodesic_poler (const double r, const double inc, const double dec);
void		vector3d_set (vector3d *cv, const double x, const double y, const double z);
void		vector3d_free (vector3d *cv);

void		vector3d_scale (vector3d *x, const double alpha);
vector3d	*vector3d_copy (const vector3d *src);
void		vector3d_axpy (const double alpha, const vector3d *x, vector3d *y);
double		vector3d_dot (const vector3d *x, const vector3d *y);
double		vector3d_nrm (const vector3d *c);

#ifdef __cplusplus
}
#endif

#endif /* VECTOR3D_H_ */
