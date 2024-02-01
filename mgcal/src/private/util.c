/*
 * util.c
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

void
error_and_exit (const char *funcname, const char *msg, const char *file, const int line)
{
	fprintf (stderr, "ERROR: %s : %s : %s %d\n", funcname, msg, file, line);
	exit (1);
}

bool
array_set_all (const int n, double *x, const double val)
{
	int		i;
	if (!x) return false;
	for (i = 0; i < n; i++) x[i] = val;
	return true;
}

bool
array_copy (const int n, double *dist, const double *src)
{
	int		i;
	if (!dist || !src) return false;
	for (i = 0; i < n; i++) dist[i] = src[i];
	return true;
}

void
set_range (double p[], const double x0, const double x1)
{
	p[0] = x0;
	p[1] = x1;
	return;
}
