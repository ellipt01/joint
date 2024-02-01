/*
 * source.h
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#ifndef SOURCE_H_
#define SOURCE_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct s_source_item	source_item;
typedef struct s_source			source;

struct s_source_item {
	vector3d	*mgz;	// magnetization
	vector3d	*pos;	// center of the magnetized body
	vector3d	*dim;	// dimension of source
	source_item	*next;
};

struct s_source {
	vector3d	*exf;	// external field
	source_item	*item;
	source_item	*begin;
	source_item	*end;
};

int		source_append_item (source *src);
source	*source_new (void);
void	source_free (source *src);
void	source_set_external (source *src, const double inc, const double dec);
void	source_set_position (source *src, const double x, const double y, const double z);
void	source_set_dimension (source *src, const double dx, const double dy, const double dz);
void	source_set_magnetization (source *src, const double mgz, const double inc, const double dec);

#ifdef __cplusplus
}
#endif

#endif /* SOURCE_H_ */
