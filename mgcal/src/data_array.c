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
		free (array);
	}
	return;
}

void
data_array_ith_copy (data_array *dest, int j, data_array *src, int i)
{
	if (j < 0 || dest->n <= j || i < 0 || src->n <= i)
		error_and_exit ("data_array_ith_copy", "specified index out of range", __FILE__, __LINE__);
	dest->x[j] = src->x[i];
	dest->y[j] = src->y[i];
	dest->z[j] = src->z[i];
	dest->data[j] = src->data[i];


}

void
data_array_copy (data_array *dest, data_array *src)
{
	int		i;
	if (dest->n != src->n)
		error_and_exit ("data_array_copy", "size of the data array invalid", __FILE__, __LINE__);
	for (i = 0; i < src->n; i++) data_array_ith_copy (dest, i, src, i);
	return;
}
