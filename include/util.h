#ifndef _UTIL_H_
#define _UTIL_H_

void		mm_real_fprintf_matrix (FILE *stream, mm_real *a, const char *format);
char		*get_toolname (char *str);
void		check_mem (const char *header);
size_t		count (const char *fn);
double		*read_terrain (const size_t c, const char *fn);

#endif // _UTIL_H_
