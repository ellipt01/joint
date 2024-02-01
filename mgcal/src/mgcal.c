/*
 * scale_factor.c
 *
 *  Created on: 2015/04/17
 *      Author: utsugi
 */

#include <stdio.h>
#include "mgcal.h"

void
mgcal_set_scale_factor (const double val)
{
	scale_factor = val;
	return;
}

double
mgcal_get_scale_factor (void)
{
	return scale_factor;
}

