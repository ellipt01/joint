/*
 * kernel.h
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#ifndef KERNEL_H_
#define KERNEL_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef double	(*mgcal_theoretical) (const vector3d *pos, const source *src, void *data);

typedef struct s_mgcal_func	mgcal_func;

struct s_mgcal_func
{
	mgcal_theoretical	function;
	void				*parameter;
};

mgcal_func	*mgcal_func_new (const mgcal_theoretical func, void *data);
void		mgcal_func_free (mgcal_func *f);
void		kernel_matrix_set (double *a, const data_array *array, const grid *g, const vector3d *mgz, const vector3d *exf, const mgcal_func *f);
double		*kernel_matrix (const data_array *array, const grid *g, const vector3d *mgz, const vector3d *exf, const mgcal_func *f);

void		kernel_matrix_jth_col_set (double *a, size_t stride, size_t j, const data_array *array, const grid *g, const vector3d *mgz, const vector3d *exf, const mgcal_func *f);
double		*kernel_matrix_jth_col_vector (size_t j, const data_array *array, const grid *g, const vector3d *mgz, const vector3d *exf, const mgcal_func *f);

void		kernel_matrix_ith_row_set (double *a, size_t stride, size_t i, const data_array *array, const grid *g, const vector3d *mgz, const vector3d *exf, const mgcal_func *f);
double		*kernel_matrix_ith_row_vector (size_t i, const data_array *array, const grid *g, const vector3d *mgz, const vector3d *exf, const mgcal_func *f);

#ifdef __cplusplus
}
#endif

#endif /* KERNEL_H_ */
