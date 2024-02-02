/*
 * calc.h
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#ifndef CALC_H_
#define CALC_H_

#ifdef __cplusplus
extern "C" {
#endif

vector3d	*dipole (const vector3d *obs, const source *s);
vector3d	*prism (const vector3d *obs, const source *s);

vector3d	*dipole_yz (const vector3d *obs, const source *s);
vector3d	*prism_yz (const vector3d *obs, const source *s);

double		total_force (const vector3d *exf, const vector3d *f);

double		x_component_dipole (const vector3d *obs, const source *src, void *data);
double		y_component_dipole (const vector3d *obs, const source *src, void *data);
double		z_component_dipole (const vector3d *obs, const source *src, void *data);
double		total_force_dipole (const vector3d *obs, const source *src, void *data);

double		x_component_prism (const vector3d *obs, const source *src, void *data);
double		y_component_prism (const vector3d *obs, const source *src, void *data);
double		z_component_prism (const vector3d *obs, const source *src, void *data);
double		total_force_prism (const vector3d *obs, const source *src, void *data);

double		y_component_dipole_yz (const vector3d *obs, const source *src, void *data);
double		z_component_dipole_yz (const vector3d *obs, const source *src, void *data);
double		total_force_dipole_yz (const vector3d *obs, const source *src, void *data);

double		y_component_prism_yz (const vector3d *obs, const source *src, void *data);
double		z_component_prism_yz (const vector3d *obs, const source *src, void *data);
double		total_force_prism_yz (const vector3d *obs, const source *src, void *data);

double		total_force_dipole_bh (const vector3d *obs, const source *src, void *data);
double		total_force_prism_bh (const vector3d *obs, const source *src, void *data);

#ifdef __cplusplus
}
#endif

#endif /* CALC_H_ */
