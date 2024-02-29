/*
 * data_array.h
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#ifndef DATA_ARRAY_H_
#define DATA_ARRAY_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct s_data_array	data_array;

struct s_data_array {
	int		n;
	double	*x;
	double	*y;
	double	*z;
	double	*data;
};

data_array	*data_array_new (const int n);
void		data_array_free (data_array *array);

void		data_array_ith_copy (data_array *dest, int j, data_array *src, int i);
void		data_array_copy (data_array *dest, data_array *src);

#ifdef __cplusplus
}
#endif

#endif /* DATA_ARRAY_H_ */
