#ifndef _GRAVITY_H_
#define _GRAVITY_H_

#include <stdio.h>
#include <mgcal.h>

#ifdef __cplusplus
extern "C" {
#endif

void	source_set_density (source *s, double rho);

double	grav_masspoint (const vector3d *obs, const source *s, void *data);
double	grav_prism (const vector3d *obs, const source *s, void *data);

double	grav_masspoint_yz (const vector3d *obs, const source *s, void *data);
double	grav_prism_yz (const vector3d *obs, const source *s, void *data);

#ifdef __cplusplus
}
#endif

#endif
