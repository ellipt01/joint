/*
 * array_array.c
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#include <stdlib.h>

#include "data_array.h"
#include "private/util.h"

static data_array *
data_array_alloc (void)
{
	data_array	*array = (data_array *) malloc (sizeof (data_array));
	array->n = 0;
	array->x = NULL;
	array->y = NULL;
	array->z = NULL;
	array->data = NULL;
	return array;
}

data_array *
data_array_new (const int n)
{
	data_array	*array = data_array_alloc ();

	array->n = n;
	array->x = (double *) malloc (n * sizeof (double));
	array->y = (double *) malloc (n * sizeof (double));
	array->z = (double *) malloc (n * sizeof (double));
	array->data = (double *) malloc (n * sizeof (double));
	array_set_all (n, array->x, 0.);
	array_set_all (n, array->y, 0.);
	array_set_all (n, array->z, 0.);
	array_set_all (n, array->data, 0.);
	return array;
}

void
data_array_free (data_array *array)
{
	if (array) {
		if (array->x) free (array->x);
		if (array->y) free (array->y);
		if (array->z) free (array->z);
		if (array->data) free (array->data);
	}
	return;
}
