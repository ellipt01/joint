/*
 * io.h
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#ifndef IO_H_
#define IO_H_

#ifdef __cplusplus
extern "C" {
#endif

data_array	*fread_data_array (FILE *stream);
void		fwrite_data_array_with_data (FILE *stream, const data_array *array, const double *data, const char *format);
void		fwrite_data_array (FILE *stream, const data_array *array, const char *format);
grid		*fread_grid (FILE *stream);
void		fwrite_grid (FILE *stream, const grid *g);
void		fwrite_grid_to_xyz (FILE *stream, const grid *g, const char *format);
void		fwrite_grid_with_data (FILE *stream, const grid *g, const double *data, const char *format);

#ifdef __cplusplus
}
#endif

#endif /* IO_H_ */
