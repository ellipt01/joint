/*
 * util.h
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#ifndef UTIL_H_
#define UTIL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

void	error_and_exit (const char *funcname, const char *msg, const char *file, const int line);
bool	array_set_all (const int n, double *x, const double val);
bool	array_copy (const int n, double *dist, const double *src);
void	set_range (double p[], const double x1, const double x2);

#ifdef __cplusplus
}
#endif

#endif /* UTIL_H_ */
