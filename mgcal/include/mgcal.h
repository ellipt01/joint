#ifndef MGCAL_H
#define MGCAL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <stdbool.h>

extern double	scale_factor;

#include "vector3d.h"
#include "data_array.h"
#include "grid.h"
#include "source.h"
#include "calc.h"
#include "io.h"
#include "kernel.h"
#include "io.h"

void	mgcal_set_scale_factor (const double val);
double	mgcal_get_scale_factor (void);

#ifdef __cplusplus
}
#endif

#endif /* MGCAL_H */
