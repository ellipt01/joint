#ifndef GRID_H
#define GRID_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct s_grid	grid;

struct s_grid {
	int			n;
	int			nh;
	int			nx;
	int			ny;
	int			nz;

	double		xrange[2];
	double		yrange[2];
	double		zrange[2];

	double		*x;
	double		*y;
	double		*z;
	double		*z1;	// irregular surface

	double		*dx;
	double		*dy;
	double		*dz;

	void		*data;
};

grid	*grid_new (const int nx, const int ny, const int nz, const double x[], const double y[], const double z[]);
grid	*grid_new_full (const int nx, const int ny, const int nz, const double x[], const double y[], const double z[], const double *dx, const double *dy, const double *dz, const double *z1);
bool	grid_set_surface (grid *g, const double *z1);
void	grid_free (grid *g);

void	grid_get_index (const grid *g, const int n, int *i, int *j, int *k, int *h);
void	grid_get_nth (const grid *g, const int n, vector3d *center, vector3d *dim);
void	grid_stretch_at_edge (grid *g, double l);

#ifdef __cplusplus
}
#endif

#endif /* GRID_H */
