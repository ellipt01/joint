#include <iostream>
#include <cstring>

#include "mmreal.h"

// fprintf mm_real in vector / matrix form
void
mm_real_fprintf_matrix (FILE *stream, mm_real *a, const char *format)
{
	for (size_t i = 0; i < a->m; i++) {
		for (size_t j = 0; j < a->n; j++) {
			fprintf (stream, format, a->data[i + j * a->m]);
			if (j < a->n - 1) fprintf (stream, " ");
		}
		fprintf (stream, "\n");
	}
}

char *
get_toolname (char *str)
{
	char	*p = strrchr (str, '/');
	if (p == NULL) p = str;
	else p++;
	return p;
}

size_t
count (const char *fn)
{
	FILE	*fp = fopen (fn, "r");
	if (!fp) {
		char	msg[80];
		sprintf (msg, "cannot open terrain file: %s", fn);
		throw std::runtime_error (msg);
	}

	int		c = 0;
	char	buf[BUFSIZ];
	while (fgets (buf, BUFSIZ, fp) != NULL) c++;
	fclose (fp);
	return c;
}

double *
read_terrain (const size_t c, const char *fn)
{
	FILE	*fp = fopen (fn, "r");
	if (!fp) {
		char	msg[80];
		sprintf (msg, "cannot open terrain file: %s", fn);
		throw std::runtime_error (msg);
	}

	double	*zsurf = new double [c];
	char	buf[BUFSIZ];
	int		k = 0;
	while (fgets (buf, BUFSIZ, fp) != NULL) {
		double	x, y, z;
		sscanf (buf, "%lf\t%lf\t%lf", &x, &y, &z);
		zsurf[k] = z;
		if (++k >= c) break;
	}
	fclose (fp);
	return zsurf;
}

