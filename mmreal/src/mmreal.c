/*
 * mmreal.c
 *
 *  Created on: 2014/06/25
 *      Author: utsugi
 */

#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "mmreal.h"
#include "_blas_.h"

static MM_INT	ione  =  1;

/* print an error message and exit */
static void
error_and_exit (const char *function_name, const char *error_msg, const char *file, int line)
{
	fprintf (stderr, "ERROR: %s: %s:%d: %s\n", function_name, file, line, error_msg);
	exit (1);
}

/* print warning message */
static void
printf_warning (const char *function_name, const char *error_msg, const char *file, int line)
{
	fprintf (stderr, "WARNING: %s: %s:%d: %s\n", function_name, file, line, error_msg);
	return;
}

/* mm_real supports real symmetric/general sparse/dense matrix */
static bool
is_type_supported (MM_typecode typecode)
{
	// invalid type
	if (!mm_is_valid (typecode)) return false;

	// pattern is not supported
	if (mm_is_pattern (typecode)) return false;

	// integer and complex matrix are not supported
	if (mm_is_integer (typecode) || mm_is_complex (typecode)) return false;

	// skew and hermitian are not supported
	if (mm_is_skew (typecode) || mm_is_hermitian (typecode)) return false;

	return true;
}

/* check format */
static bool
is_format_valid (MMRealFormat format) {
	return (format == MM_REAL_SPARSE || format == MM_REAL_DENSE);
}

/* check symmetric */
static bool
is_symm_valid (MMRealSymm symm)
{
	return (symm == MM_REAL_GENERAL || symm == MM_REAL_SYMMETRIC_UPPER
			|| symm == MM_REAL_SYMMETRIC_LOWER);
}

/* allocate mm_real */
static mm_real *
mm_real_alloc (void)
{
	mm_real	*x = (mm_real *) malloc (sizeof (mm_real));
	if (x == NULL) return NULL;

	x->m = 0;
	x->n = 0;
	x->nnz = 0;
	x->i = NULL;
	x->p = NULL;
	x->data = NULL;

	x->symm = MM_REAL_GENERAL;

	/* set typecode = "M_RG" : Matrix Real General */
	// typecode[3] = 'G' : General
	mm_initialize_typecode (&x->typecode);
	// typecode[0] = 'M' : Matrix
	mm_set_matrix (&x->typecode);
	// typecode[2] = 'R' : Real
	mm_set_real (&x->typecode);

	x->owner = true;

	return x;
}

/*** create new mm_real object
 * MMRealFormat	format: MM_REAL_DENSE or MM_REAL_SPARSE
 * MMRealSymm		symm  : MM_REAL_GENERAL, MM_REAL_SYMMETRIC_UPPER or MM_REAL_SYMMETRIC_LOWER
 * MM_INT			m, n  : rows and columns of the matrix
 * MM_INT			nnz   : number of nonzero elements of the matrix ***/
static void
mm_real_set (mm_real *x, MMRealFormat format, MMRealSymm symm, MM_INT m, MM_INT n, MM_INT nnz)
{
	bool	symmetric;

	if (x == NULL) error_and_exit ("mm_real_new", "failed to allocate object.", __FILE__, __LINE__);
	x->m = m;
	x->n = n;
	x->nnz = nnz;

	// typecode[1] = 'C' or 'A'
	if (format == MM_REAL_SPARSE) {
		mm_set_coordinate (&x->typecode);
		if (nnz > 0) {
			x->i = (MM_INT *) malloc (x->nnz * sizeof (MM_INT));
			if (x->i == NULL) error_and_exit ("mm_real_new", "cannot allocate memory.", __FILE__, __LINE__);
		}
		x->p = (MM_INT *) malloc ((x->n + 1) * sizeof (MM_INT));
		if (x->p == NULL) error_and_exit ("mm_real_new", "cannot allocate memory.", __FILE__, __LINE__);
		// initialize x->p[0]
		for (size_t j = 0; j <= n; j++) x->p[j] = 0;
	} else mm_set_array (&x->typecode);

	x->symm = symm;
	// typecode[3] = 'G' -> 'S'
	symmetric = symm & MM_SYMMETRIC;
	if (symmetric) mm_set_symmetric (&x->typecode);


	return;
}

/*** create mm_real object ***/
mm_real *
mm_real_new (MMRealFormat format, MMRealSymm symm, MM_INT m, MM_INT n, MM_INT nnz)
{
	mm_real	*x;

	if (!is_format_valid (format))
		error_and_exit ("mm_real_new", "invalid MMRealFormat format.", __FILE__, __LINE__);
	if (!is_symm_valid (symm))
		error_and_exit ("mm_real_new", "invalid MMRealSymm symm.", __FILE__, __LINE__);

	if ((symm & MM_SYMMETRIC) && m != n)
		error_and_exit ("mm_real_new", "symmetric matrix must be square.", __FILE__, __LINE__);

	x = mm_real_alloc ();
	if (x == NULL) error_and_exit ("mm_real_new", "failed to allocate object.", __FILE__, __LINE__);

	mm_real_set (x, format, symm, m, n, nnz);

	if (!is_type_supported (x->typecode)) {
		char	msg[128];
		sprintf (msg, "matrix type does not supported :[%s].", mm_typecode_to_str (x->typecode));
		error_and_exit ("mm_real_new", msg, __FILE__, __LINE__);
	}

	// allocate arrays
	if (nnz > 0) {
		x->data = (MM_DBL *) malloc (x->nnz * sizeof (MM_DBL));
		if (x->data == NULL) error_and_exit ("mm_real_new", "cannot allocate memory.", __FILE__, __LINE__);
	}
	return x;
}


/*** create mm_real view array object ***/
mm_real *
mm_real_view_array (MMRealFormat format, MMRealSymm symm, MM_INT m, MM_INT n, MM_INT nnz, MM_DBL *data)
{
	mm_real	*x;

	if (!is_format_valid (format))
		error_and_exit ("mm_real_view_array", "invalid MMRealFormat format.", __FILE__, __LINE__);
	if (!is_symm_valid (symm))
		error_and_exit ("mm_real_view_array", "invalid MMRealSymm symm.", __FILE__, __LINE__);

	if ((symm & MM_SYMMETRIC) && m != n)
		error_and_exit ("mm_real_view_array", "symmetric matrix must be square.", __FILE__, __LINE__);

	x = mm_real_alloc ();
	if (x == NULL) error_and_exit ("mm_real_view_array", "failed to allocate object.", __FILE__, __LINE__);

	mm_real_set (x, format, symm, m, n, nnz);

	if (!is_type_supported (x->typecode)) {
		char	msg[128];
		sprintf (msg, "matrix type does not supported :[%s].", mm_typecode_to_str (x->typecode));
		error_and_exit ("mm_real_view_array", msg, __FILE__, __LINE__);
	}

	x->data = data;
	x->owner = false;

	return x;
}

/*** free mm_real ***/
void
mm_real_free (mm_real *x)
{
	if (x) {
		if (x->p) free (x->p);
		if (x->i) free (x->i);
		if (x->owner && x->data) free (x->data);
		free (x);
	}
	return;
}

/*** reallocate mm_real ***/
bool
mm_real_realloc (mm_real *x, MM_INT nnz)
{
	if (x->nnz == nnz) return true;
	x->data = (MM_DBL *) realloc (x->data, nnz * sizeof (MM_DBL));
	if (x->data == NULL) return false;
	if (mm_real_is_sparse (x)) {
		x->i = (MM_INT *) realloc (x->i, nnz * sizeof (MM_INT));
		if (x->i == NULL) return false;
	}
	x->nnz = nnz;
	return true;
}

/* set to general */
static void
mm_real_set_general (mm_real *x)
{
	if (!mm_real_is_symmetric (x)) return;
	mm_set_general (&(x->typecode));
	x->symm = MM_REAL_GENERAL;
	return;
}

/* set to symmetric
 * by default, assume symmetric upper
 * i.e., x->symm is set to MM_SYMMETRIC | MM_UPPER */
static void
mm_real_set_symmetric (mm_real *x)
{
	if (mm_real_is_symmetric (x)) return;
	mm_set_symmetric (&(x->typecode));
	x->symm = MM_REAL_SYMMETRIC_UPPER;	// by default, assume symmetric upper
	return;
}

/* set to symmetric upper */
static void
mm_real_set_upper (mm_real *x)
{
	if (mm_real_is_upper (x)) return;
	x->symm = MM_REAL_SYMMETRIC_UPPER;
	return;
}

/* set to symmetric lower */
static void
mm_real_set_lower (mm_real *x)
{
	if (mm_real_is_lower (x)) return;
	x->symm = MM_REAL_SYMMETRIC_LOWER;
	return;
}

/* sort each columns of mm_sparse */
typedef struct {
	MM_INT	i;
	MM_DBL	data;
} matrix_element;

static int
compare_row_index (const void *a, const void *b)
{
	matrix_element	*_a = (matrix_element *) a;
	matrix_element	*_b = (matrix_element *) b;
	return (int) (_a->i - _b->i);
}

static void
mm_real_sort_sparse (mm_sparse *s)
{
	MM_INT	j;
	MM_INT	*si = s->i;
	MM_DBL	*sdata = s->data;
	MM_INT	*sp = s->p;
	matrix_element	t[s->m];

	for (j = 0; j < s->n; j++) {
		MM_INT		k;
		MM_INT	p = sp[j];
		MM_INT	n = sp[j + 1] - p;

		if (n < 2) continue;

		for (k = 0; k < n; k++) {
			t[k].i = si[k + p];
			t[k].data = sdata[k + p];
		}

		qsort (t, n, sizeof (matrix_element), compare_row_index);

		for (k = 0; k < n; k++) {
			si[k + p] = t[k].i;
			sdata[k + p] = t[k].data;
		}

	}
	return;
}

/*** sort mm_real ***/
void
mm_real_sort (mm_real *x)
{
	if (mm_real_is_sparse (x)) mm_real_sort_sparse (x);
	return;
}

/* memcpy sparse */
static void
mm_real_memcpy_sparse (mm_sparse *dest, mm_sparse *src)
{
	MM_INT		k;
	MM_INT		n = src->n;
	MM_INT		nnz = src->nnz;

	MM_INT		*si = src->i;
	MM_INT		*sp = src->p;
	MM_INT		*di = dest->i;
	MM_INT		*dp = dest->p;

	for (k = 0; k < nnz; k++) di[k] = si[k];
	for (k = 0; k <= n; k++) dp[k] = sp[k];
	dcopy_ (&nnz, src->data, &ione, dest->data, &ione);

	return;
}

/* memcpy dense */
static void
mm_real_memcpy_dense (mm_dense *dest, mm_dense *src)
{
	dcopy_ (&src->nnz, src->data, &ione, dest->data, &ione);
	return;
}

/*** memcpy mm_real ***/
void
mm_real_memcpy (mm_real *dest, mm_real *src)
{
	if (mm_real_is_sparse (src)) {
		if (!mm_real_is_sparse (dest))
			error_and_exit ("mm_real_memcpy", "destination matrix format does not match source matrix format.", __FILE__, __LINE__);
		mm_real_memcpy_sparse (dest, src);
	} else {
		if (!mm_real_is_dense (dest))
			error_and_exit ("mm_real_memcpy", "destination matrix format does not match source matrix format.", __FILE__, __LINE__);
		mm_real_memcpy_dense (dest, src);
	}
	return;
}


/* copy sparse */
static mm_sparse *
mm_real_copy_sparse (mm_sparse *src)
{
	mm_sparse	*dest = mm_real_new (MM_REAL_SPARSE, src->symm, src->m, src->n, src->nnz);
	mm_real_memcpy_sparse (dest, src);
	return dest;
}

/* copy dense */
static mm_dense *
mm_real_copy_dense (mm_dense *src)
{
	mm_dense	*dest = mm_real_new (MM_REAL_DENSE, src->symm, src->m, src->n, src->nnz);
	mm_real_memcpy_dense (dest, src);
	return dest;
}

/*** copy x ***/
mm_real *
mm_real_copy (mm_real *x)
{
	return (mm_real_is_sparse (x)) ? mm_real_copy_sparse (x) : mm_real_copy_dense (x);
}

/*** convert sparse -> dense ***/
mm_dense *
mm_real_copy_sparse_to_dense (mm_sparse *s)
{
	MM_INT		j, k;
	MM_INT		p, n;
	mm_dense	*d;

	MM_INT		*si;
	MM_DBL		*sd;
	MM_INT		*sp = s->p;
	MM_DBL		*dd;

	if (mm_real_is_dense (s)) return mm_real_copy (s);
	d = mm_real_new (MM_REAL_DENSE, s->symm, s->m, s->n, s->m * s->n);
	mm_real_set_all (d, 0.);

	dd = d->data;
	for (j = 0; j < s->n; j++) {
		p = sp[j];
		n = sp[j + 1] - p;
		si = s->i + p;
		sd = s->data + p;
		dd = d->data + j * s->m;
		for (k = 0; k < n; k++) dd[si[k]] = sd[k];
	}
	return d;
}

/*** convert dense -> sparse
 * if fabs (x->data[j]) < threshold, set to 0 ***/
mm_sparse *
mm_real_copy_dense_to_sparse (mm_dense *d, MM_DBL threshold)
{
	MM_INT		i, j, k;
	mm_sparse	*s;

	MM_INT		*si;
	MM_INT		*sp;
	MM_DBL		*sd;

	MM_INT		dm = d->m;
	MM_INT		dn = d->n;
	MM_DBL		*dd = d->data;
	MM_DBL		dij;

	if (mm_real_is_sparse (d)) return mm_real_copy (d);
	s = mm_real_new (MM_REAL_SPARSE, d->symm, d->m, d->n, d->nnz);

	si = s->i;
	sp = s->p;
	sd = s->data;

	k = 0;
	if (!mm_real_is_symmetric (d)) {
		for (j = 0; j < dn; j++) {
			for (i = 0; i < dm; i++) {
				dij = dd[i + j * dm];
				if (fabs (dij) >= threshold) {
					si[k] = i;
					sd[k] = dij;
					k++;
				}
			}
			sp[j + 1] = k;
		}
	} else {
		if (mm_real_is_upper (d)) {
			for (j = 0; j < dn; j++) {
				for (i = 0; i < j + 1; i++) {
					dij = dd[i + j * dm];
					if (fabs (dij) >= threshold) {
						si[k] = i;
						sd[k] = dij;
						k++;
					}
				}
				sp[j + 1] = k;
			}
		} else if (mm_real_is_lower (d)) {
			for (j = 0; j < dn; j++) {
				for (i = j; i < dm; i++) {
					dij = dd[i + j * dm];
					if (fabs (dij) >= threshold) {
						si[k] = i;
						sd[k] = dij;
						k++;
					}
				}
				sp[j + 1] = k;
			}
		}
	}
	if (s->nnz != k) {
		mm_real_realloc (s, k);
		s->nnz = k;
	}
	return s;
}

/* set all elements of array to val */
static void
mm_real_array_set_all (MM_INT nnz, MM_DBL *data, MM_DBL val)
{
	MM_INT		k;
	for (k = 0; k < nnz; k++) data[k] = val;
	return;
}

/*** set x->data to val ***/
void
mm_real_set_all (mm_real *x, MM_DBL val)
{
	mm_real_array_set_all (x->nnz, x->data, val);
	return;
}

/*** convert sparse to dense ***/
bool
mm_real_sparse_to_dense (mm_sparse *s)
{
	MM_INT	j, k;

	MM_INT	m;
	MM_INT	n;
	MM_INT	nnz;

	MM_DBL	*tmp_data;

	MM_DBL	*td;
	MM_INT	*si;
	MM_INT	*sp;
	MM_DBL	*sd;

	if (!mm_real_is_sparse (s)) return false;

	m = s->m;
	n = s->n;
	nnz = m * n;

	/* copy s->data */
	tmp_data = (MM_DBL *) malloc (s->nnz * sizeof (MM_DBL));
	for (k = 0; k < s->nnz; k++) tmp_data[k] = s->data[k];
	td = tmp_data;

	/* reallocate s->data */
	free (s->data);
	s->data = (MM_DBL *) malloc (nnz * sizeof (MM_DBL));
	mm_real_array_set_all (nnz, s->data, 0.);

	/* convert sparse to dense */
	si = s->i;
	sp = s->p;
	sd = s->data;
	for (j = 0; j < n; j++) {
		MM_INT	np = sp[1] - *sp;
		for (k = 0; k < np; k++) {
			sd[*si] = *td;
			si++;
			td++;
		}
		sd += s->m;
		sp++;
	}
	free (tmp_data);

	mm_set_array (&s->typecode);
	s->nnz = nnz;
	free (s->i);
	s->i = NULL;
	free (s->p);
	s->p = NULL;

	return true;
}

/*** convert dense to sparse ***/
bool
mm_real_dense_to_sparse (mm_dense *d, MM_DBL threshold)
{
	MM_INT	i, j, k;

	MM_INT	m;
	MM_INT	n;
	MM_INT	nnz;

	MM_INT	*di;
	MM_INT	*dp;
	MM_DBL	*dd;

	if (!mm_real_is_dense (d)) return false;

	m = d->m;
	n = d->n;
	nnz = d->nnz;

	/* convert dense to sparse */
	d->i = (MM_INT *) malloc (nnz * sizeof (MM_INT));
	d->p = (MM_INT *) malloc ((n + 1) * sizeof (MM_INT));
	di = d->i;
	dp = d->p;
	dd = d->data;

	k = 0;
	*dp = 0;

	if (!mm_real_is_symmetric (d)) {
		MM_DBL	*d0 = d->data;
		for (j = 0; j < n; j++) {
			for (i = 0; i < m; i++) {
				MM_DBL	val = *d0;
				if (fabs (val) > threshold) {
					*(di++) = i;
					*(dd++) = val;
					k++;
				}
				d0++;
			}
			*(++dp) = k;
		}
	} else {
		if (mm_real_is_upper (d)) {
			for (j = 0; j < n; j++) {
				for (i = 0; i <= j; i++) {
					MM_DBL	val = dd[i + j * m];
					if (fabs (val) > threshold) {
						di[k] = i;
						dd[k] = val;
						k++;
					}
				}
				dp[j + 1] = k;
			}
		} else if (mm_real_is_lower (d)) {
			for (j = 0; j < n; j++) {
				for (i = j; i < m; i++) {
					MM_DBL	val = dd[i + j * m];
					if (fabs (val) > threshold) {
						di[k] = i;
						dd[k] = val;
						k++;
					}
				}
				dp[j + 1] = k;
			}
		} else return false;
	}

	mm_set_coordinate (&d->typecode);
	if (k != d->nnz) {
		mm_real_realloc (d, k);
		d->nnz = k;
	}
	return true;
}

/* binary search: search the value key from array *s of length n.
 * if found, return its index otherwise -1 */
static MM_INT
bin_search (MM_INT key, MM_INT *s, MM_INT n)
{
	MM_INT	start = 0;
	MM_INT	end = n - 1;
	MM_INT	mid;
	while (start < end) {
		mid = (MM_INT) (start + end) / 2;
		if (key <= s[mid]) end = mid;
		else start = mid + 1;
	}
	return (key == s[end]) ? end : -1;
}

/* find element that s->i[l] = j in the k-th column of s and return its index l */
static MM_INT
find_row_element (MM_INT j, mm_sparse *s, MM_INT k)
{
	MM_INT	p = s->p[k];
	MM_INT	n = s->p[k + 1] - p;
	MM_INT	*si = s->i + p;
	MM_INT	res = bin_search (j, si, n);
	return (res < 0) ? -1 : res + p;
}

/* convert sparse symmetric -> sparse general */
static mm_sparse *
mm_real_symmetric_to_general_sparse (mm_sparse *x)
{
	MM_INT		i, j;
	mm_sparse	*s;

	MM_INT		*si;
	MM_INT		*sp;
	MM_DBL		*sd;
	MM_INT		xn = x->n;
	MM_INT		*xi = x->i;
	MM_INT		*xp = x->p;
	MM_DBL		*xd = x->data;

	if (!mm_real_is_symmetric (x)) return mm_real_copy (x);
	s = mm_real_new (MM_REAL_SPARSE, MM_REAL_GENERAL, x->m, x->n, 2 * x->nnz);

	si = s->i;
	sp = s->p;
	sd = s->data;

	i = 0;
	for (j = 0; j < xn; j++) {
		MM_INT		k;
		MM_INT		pend = xp[j + 1];
		if (mm_real_is_upper (x)) {
			for (k = xp[j]; k < pend; k++) {
				si[i] = xi[k];
				sd[i++] = xd[k];
			}
			for (k = j + 1; k < xn; k++) {
				MM_INT	l = find_row_element (j, x, k);
				// if found
				if (l >= 0) {
					si[i] = k;
					sd[i++] = xd[l];
				}
			}
		} else if (mm_real_is_lower (x)) {
			for (k = 0; k < j; k++) {
				MM_INT	l = find_row_element (j, x, k);
				// if found
				if (l >= 0) {
					si[i] = k;
					sd[i++] = xd[l];
				}
			}
			for (k = xp[j]; k < pend; k++) {
				si[i] = xi[k];
				sd[i++] = xd[k];
			}
		}
		sp[j + 1] = i;
	}
	if (s->nnz != i) mm_real_realloc (s, i);
	return s;
}

/* convert dense symmetric -> dense general */
static void
mm_real_symmetric_to_general_dense (mm_dense *d)
{
	MM_INT	j;
	if (mm_real_is_upper (d)) {
		for (j = 0; j < d->n; j++) {
			MM_INT	i0 = j + 1;
			MM_INT	n = d->m - i0;
			dcopy_ (&n, d->data + j + i0 * d->m, &d->m, d->data + i0 + j * d->m, &ione);
		}
	} else {
		for (j = 0; j < d->n; j++) dcopy_ (&j, d->data + j, &d->m, d->data + j * d->m, &ione);
	}
	mm_real_set_general (d);
	return;
}

/*** convert symmetric -> general ***/
bool
mm_real_symmetric_to_general (mm_real *x)
{
	if (!mm_real_is_symmetric (x)) return false;

	if (mm_real_is_sparse (x)) {
		mm_sparse	*s = mm_real_symmetric_to_general_sparse (x);
		if (s->nnz != x->nnz) mm_real_realloc (x, s->nnz);
		mm_real_memcpy_sparse (x, s);
		mm_real_free (s);
	} else {
		mm_real_symmetric_to_general_dense (x);
	}

	return true;
}

/* copy general to symmetric upper sparse */
static mm_sparse *
mm_real_general_to_symmetric_upper_sparse (mm_sparse *x)
{
	int		j, k;
	mm_real	*s;

	MM_INT	*xi = x->i;
	MM_INT	*xp = x->p;
	MM_DBL	*xd = x->data;

	MM_INT	*si;
	MM_INT	*sp;
	MM_DBL	*sd;

	MM_INT	nnz = 0;
	for (j = 0; j < x->n; j++) {
		int		p0 = xp[j];
		int		p1 = xp[j + 1];
		int		n;
		for (n = p0; n < p1; n++) {
			if (xi[n] > j) break;
			nnz++;
		}
	}
	s = mm_real_new (MM_REAL_SPARSE, MM_REAL_SYMMETRIC_UPPER, x->m, x->n, nnz);
	si = s->i;
	sd = s->data;
	sp = &s->p[0];

	k = 0;
	*sp = 0;
	for (j = 0; j < x->n; j++) {
		int		p0 = xp[j];
		int		p1 = xp[j + 1];
		int		n;
		for (n = p0; n < p1; n++) {
			if (xi[n] > j) break;
			si[k] = xi[n];
			sd[k] = xd[n];
			if (++k > nnz) break;
		}
		*(++sp) = k;
	}

	return s;
}

/* copy general to symmetric lower sparse */
static mm_sparse *
mm_real_general_to_symmetric_lower_sparse (mm_sparse *x)
{
	int		j, k;
	mm_real	*s;

	MM_INT	*xi = x->i;
	MM_INT	*xp = x->p;
	MM_DBL	*xd = x->data;

	MM_INT	*si;
	MM_INT	*sp;
	MM_DBL	*sd;

	MM_INT	nnz = 0;
	for (j = 0; j < x->n; j++) {
		int		p0 = xp[j];
		int		p1 = xp[j + 1];
		int		n;
		for (n = p0; n < p1; n++) {
			if (xi[n] >= j) nnz++;
		}
	}
	s = mm_real_new (MM_REAL_SPARSE, MM_REAL_SYMMETRIC_LOWER, x->m, x->n, nnz);
	si = s->i;
	sd = s->data;
	sp = &s->p[0];

	k = 0;
	*sp = 0;
	for (j = 0; j < x->n; j++) {
		int		p0 = xp[j];
		int		p1 = xp[j + 1];
		int		n;
		for (n = p0; n < p1; n++) {
			if (xi[n] >= j) {
				si[k] = xi[n];
				sd[k] = xd[n];
				if (++k > nnz) break;
			}
		}
		*(++sp) = k;
	}

	return s;
}

/* convert general to symmetric sparse */
static bool
mm_real_general_to_symmetric_sparse (const char uplo, mm_sparse *s)
{
	mm_sparse	*x;
	if (uplo == 'u' || uplo == 'U') x = mm_real_general_to_symmetric_upper_sparse (s);
	else if (uplo == 'l' || uplo == 'L') x = mm_real_general_to_symmetric_lower_sparse (s);
	else return false;

	mm_real_realloc (s, x->nnz);
	mm_real_memcpy (s, x);
	mm_real_free (x);

	mm_set_symmetric (&(s->typecode));
	if (uplo == 'u' || uplo == 'U') s->symm = MM_REAL_SYMMETRIC_UPPER;
	else if (uplo == 'l' || uplo == 'L') s->symm = MM_REAL_SYMMETRIC_LOWER;

	return true;
}

/* convert general to symmetric dense */
static bool
mm_real_general_to_symmetric_dense (const char uplo, mm_dense *d)
{
	mm_set_symmetric (&(d->typecode));
	if (uplo == 'u' || uplo == 'U') d->symm = MM_REAL_SYMMETRIC_UPPER;
	else if (uplo == 'l' || uplo == 'L') d->symm = MM_REAL_SYMMETRIC_LOWER;
	else return false;

	return true;
}

/*** convert general to symmetric ***/
bool
mm_real_general_to_symmetric (const char uplo, mm_real *x)
{
	bool	ret = true;
	if (mm_real_is_symmetric (x)) {
		printf_warning ("mm_real_general_to_symmetric", "matrix is symmetric", __FILE__, __LINE__);
		return false;
	}
	if (uplo != 'u' && uplo != 'U' && uplo != 'l' && uplo != 'L')
		error_and_exit ("mm_real_general_to_symmetric", "uplo must be l, L or u, U", __FILE__, __LINE__);
	if (mm_real_is_sparse (x)) ret = mm_real_general_to_symmetric_sparse (uplo, x);
	else if (mm_real_is_dense (x)) ret = mm_real_general_to_symmetric_dense (uplo, x);
	else return false;

	return ret;
}

/* iamax of sparse */
static MM_INT
mm_real_sj_iamax (mm_sparse *s, MM_INT j)
{
	MM_INT	n = s->p[j + 1] - s->p[j];
	if (n <= 0) return -1;
	return idamax_ (&n, s->data + s->p[j], &ione) - 1;
}

/* iamax of dense */
static MM_INT
mm_real_dj_iamax (mm_dense *d, MM_INT j)
{
	return idamax_ (&d->m, d->data + j * d->m, &ione) - 1;
}

/*** iamax ***/
MM_INT
mm_real_iamax (mm_real *x)
{
	return idamax_ (&x->nnz, x->data, &ione) - 1;
}

/*** iamax ***/
MM_INT
mm_real_xj_iamax (mm_real *x, MM_INT j)
{
	if (j < 0 || x->n <= j) error_and_exit ("mm_real_xj_iamax", "index out of range.", __FILE__, __LINE__);
	return (mm_real_is_sparse (x)) ? mm_real_sj_iamax (x, j) : mm_real_dj_iamax (x, j);
}

/* extract column of sparse matrix */
static mm_dense *
mm_real_sj_col (mm_sparse *s, MM_INT j)
{
	MM_INT	l, p, n;
	MM_INT	*si;
	MM_DBL	*sd;

	mm_dense	*dj = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, s->m, 1, s->m);
	mm_real_set_all (dj, 0.);
	p = *(s->p + j);
	n = *(s->p + j + 1) - p;
	si = s->i + p;
	sd = s->data + p;
	for (l = 0; l < n; l++, si++, sd++) dj->data[*si] = *sd;

	if (mm_real_is_symmetric (s)) {
		if (mm_real_is_upper (s)) {
#pragma omp parallel for
			for (l = j + 1; l < s->n; l++) {
				MM_INT	k = find_row_element (j, s, l);
				if (k >= 0) dj->data[l] = s->data[k];
			}
		} else if (mm_real_is_lower (s)) {
#pragma omp parallel for
			for (l = 0; l < j; l++) {
				MM_INT	k = find_row_element (j, s, l);
				if (k >= 0) dj->data[l] = s->data[k];
			}
		}
	}
	return dj;
}

/* extract column of dense matrix */
static mm_dense *
mm_real_dj_col (mm_dense *d, MM_INT j)
{
	mm_dense	*c = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, d->m, 1, d->m);
	if (!mm_real_is_symmetric (d)) {
		dcopy_ (&d->m, d->data + j * d->m, &ione, c->data, &ione);
	} else {
		MM_INT	h, n;
		h = (MM_INT) (d->n / 2);
		if (d->n % 2 == 0) h--;
		n = d->n - h;
		if (mm_real_is_upper (d)) {
			dcopy_ (&h, d->data + j * d->m, &ione, c->data, &ione);
			dcopy_ (&n, d->data + j + h * d->m, &d->m, c->data + h, &ione);
		} else {
			dcopy_ (&h, d->data + j, &d->m, c->data, &ione);
			dcopy_ (&n, d->data + j + h * d->m, &ione, c->data + h, &ione);
		}
	}
	return c;
}

/*** extract column of matrix ***/
mm_dense *
mm_real_xj_col (mm_real *x, MM_INT j)
{
	if (j < 0 || x->n <= j) error_and_exit ("mm_real_xj_col", "index out of range.", __FILE__, __LINE__);
//	if (mm_real_is_symmetric (x)) error_and_exit ("mm_real_xj_col", "x must be general.", __FILE__, __LINE__);
	return (mm_real_is_sparse (x)) ? mm_real_sj_col (x, j) : mm_real_dj_col (x, j);
}

/* extract row of sparse matrix */
static mm_dense *
mm_real_si_row (mm_sparse *s, MM_INT i)
{
	MM_INT	l;
	MM_INT	*si = s->i;
	MM_INT	*sp = s->p;
	MM_DBL	*sd = s->data;

	mm_real	*di = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, 1, s->n, s->n);
	mm_real_set_all (di, 0.);

#pragma omp parallel for
	for (MM_INT k = 0; k < s->n; k++) {
		l = find_row_element (i, s, k);
		if (l >= 0 && l < s->nnz) di->data[k] = sd[l];	// found
	}

	if (mm_real_is_symmetric (s)) {
		MM_INT	p, n;
		p = *(sp + i);
		n = *(sp + i + 1) - p;
		si = s->i + p;
		sd = s->data + p;
		//for (l = 0; l < n; l++, si++, sd++) di->data[*si] = *sd;
		for (l = 0; l < n; l++) di->data[si[l]] = sd[l];
	}

	return di;
}

#ifdef DEBUG
static mm_dense *
mm_real_si_row1 (mm_sparse *s, MM_INT i)
{
	MM_INT	*si = s->i;
	MM_INT	*sp = s->p;
	MM_DBL	*data = s->data;

	mm_real	*r = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, 1, s->n, s->n);
	mm_real_set_all (r, 0.);

#pragma omp parallel for
	for (MM_INT j = 0; j < s->n; j++) {
		for (MM_INT k = sp[j]; k < sp[j + 1]; k++) {
			MM_INT	sik = si[k];
			if (sik == i) {
				r->data[j] = data[k];
				break;
			} else if (sik > i) break;
		}
	}
	return r;
}
#endif

/* extract row of dense matrix */
static mm_dense *
mm_real_di_row (mm_dense *d, MM_INT i)
{
	mm_real	*r = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, 1, d->n, d->n);

	if (!mm_real_is_symmetric (d)) {
		dcopy_ (&d->n, d->data + i, &d->m, r->data, &ione);
	} else {
		MM_INT	h, n;
		h = (MM_INT) (d->n / 2);
		if (d->n % 2 == 0) h--;
		n = d->n - h;
		if (mm_real_is_upper (d)) {
			dcopy_ (&h, d->data + i * d->m, &ione, r->data, &ione);
			dcopy_ (&n, d->data + i + h * d->m, &d->m, r->data + h, &ione);
		} else {
			dcopy_ (&h, d->data + i, &d->m, r->data, &ione);
			dcopy_ (&n, d->data + i + h * d->m, &ione, r->data + h, &ione);
		}
	}
	return r;
}

/*** extract row of matrix ***/
mm_dense *
mm_real_xi_row (mm_real *x, MM_INT i)
{
	if (i < 0 || x->m <= i) error_and_exit ("mm_real_xi_row", "index out of range.", __FILE__, __LINE__);
//	if (mm_real_is_symmetric (x)) error_and_exit ("mm_real_xi_row", "x must be general.", __FILE__, __LINE__);
	return (mm_real_is_sparse (x)) ? mm_real_si_row (x, i) : mm_real_di_row (x, i);
}

/* transpose dense matrix */
static void
mm_real_transpose_dense (mm_dense *d)
{
	MM_INT	i, j;
	MM_INT	m = d->m;
	MM_INT	n = d->n;
	if (m > 1 && n > 1) {
		MM_DBL	*data = (MM_DBL *) malloc (m * n * sizeof (MM_DBL));
		for (j = 0; j < d->n; j++) {
			for (i = 0; i < d->m; i++) data[j + i * n] = d->data[i + j * m];
		}
		for (i = 0; i < d->nnz; i++) d->data[i] = data[i];
		free (data);
	}
	d->m = n;
	d->n = m;
	return;
}

/* transpose sparse matrix */
static void
mm_real_transpose_sparse (mm_sparse *s)
{
	MM_INT	i, j, k;
	MM_INT	p;
	MM_INT	m = s->m;
	MM_INT	n = s->n;
	MM_INT	nnz = s->nnz;

	MM_INT	*ti = (MM_INT *) malloc (nnz * sizeof (MM_INT));
	MM_INT	*tp = (MM_INT *) malloc ((m + 1) * sizeof (MM_INT));
	MM_DBL	*data = (MM_DBL *) malloc (nnz * sizeof (MM_DBL));

	p = 0;
	tp[0] = 0;
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			MM_INT	pj = s->p[j];
			MM_INT	l = s->p[j + 1] - pj;
			for (k = 0; k < l; k++) {
				MM_INT	ski = s->i[pj + k];
				if (ski == i) {
					ti[p] = j;
					data[p] = s->data[pj + k];
					p++;
					break;
				} else if (ski > i) break;
			}
		}
		tp[i + 1] = p;
	}
	for (k = 0; k < nnz; k++) s->i[k] = ti[k];
	dcopy_ (&s->nnz, data, &ione, s->data, &ione);
	s->p = (MM_INT *) realloc (s->p, (m + 1) * sizeof (MM_INT));
	for (k = 0; k <= m; k++) s->p[k] = tp[k];

	s->m = n;
	s->n = m;

	free (ti);
	free (tp);
	free (data);

	return;
}

/*** transpose matrix ***/
void
mm_real_transpose (mm_real *x)
{
	(mm_real_is_sparse (x)) ? mm_real_transpose_sparse (x) : mm_real_transpose_dense (x);
}


/* identity sparse matrix */
static mm_sparse *
mm_real_seye (MM_INT n)
{
	MM_INT		k;
	mm_sparse	*s = mm_real_new (MM_REAL_SPARSE, MM_REAL_GENERAL, n, n, n);
	MM_INT		*si = s->i;
	MM_INT		*sp = s->p;
	MM_DBL		*sd = s->data;
	for (k = 0; k < n; k++) {
		si[k] = k;
		sd[k] = 1.;
		sp[k + 1] = k + 1;
	}
	return s;
}

/* identity dense matrix */
static mm_dense *
mm_real_deye (MM_INT n)
{
	MM_INT		k;
	mm_dense	*d = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, n, n, n * n);
	MM_DBL		*dd = d->data;
	mm_real_set_all (d, 0.);
	for (k = 0; k < n; k++) dd[k + k * n] = 1.;
	return d;
}

/*** n x n identity matrix ***/
mm_real *
mm_real_eye (MMRealFormat format, MM_INT n)
{
	if (n <= 0) error_and_exit ("mm_real_eye", "invalid size.", __FILE__, __LINE__);
	return (format == MM_REAL_SPARSE) ? mm_real_seye (n) : mm_real_deye (n);
}

/* s = [s1; s2] */
static mm_sparse *
mm_real_vertcat_sparse (mm_sparse *s1, mm_sparse *s2)
{
	MM_INT		i, j, k;
	MM_INT		m = s1->m + s2->m;
	MM_INT		n = s1->n;
	MM_INT		nnz = s1->nnz + s2->nnz;
	mm_sparse	*s = mm_real_new (MM_REAL_SPARSE, MM_REAL_GENERAL, m, n, nnz);

	MM_INT		*si = s->i;
	MM_INT		*sp = s->p;
	MM_INT		*si1 = s1->i;
	MM_INT		*sp1 = s1->p;
	MM_INT		*si2 = s2->i;
	MM_INT		*sp2 = s2->p;

	k = 0;
	for (j = 0; j < n; j++) {
		MM_INT	n1 = sp1[j + 1] - sp1[j];
		MM_INT	n2 = sp2[j + 1] - sp2[j];
		MM_INT	pend;
		dcopy_ (&n1, s1->data + sp1[j], &ione, s->data + k, &ione);
		pend = sp1[j + 1];
		for (i = sp1[j]; i < pend; i++) si[k++] = si1[i];
		dcopy_ (&n2, s2->data + sp2[j], &ione, s->data + k, &ione);
		pend = sp2[j + 1];
		for (i = sp2[j]; i < pend; i++) si[k++] = si2[i] + s1->m;
		sp[j + 1] = k;
	}
	return s;
}

/* d = [d1; d2] */
static mm_dense *
mm_real_vertcat_dense (mm_dense *d1, mm_dense *d2)
{
	MM_INT		j;
	MM_INT		m = d1->m + d2->m;
	MM_INT		n = d1->n;
	MM_INT		nnz = d1->nnz + d2->nnz;
	mm_dense	*d = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, n, nnz);

	for (j = 0; j < n; j++) {
		dcopy_ (&d1->m, d1->data + j * d1->m, &ione, d->data + j * d->m, &ione);
		dcopy_ (&d2->m, d2->data + j * d2->m, &ione, d->data + j * d->m + d1->m, &ione);
	}
	return d;
}

/*** x = [x1; x2] ***/
mm_real *
mm_real_vertcat (mm_real *x1, mm_real *x2)
{
	if ((mm_real_is_sparse (x1) && mm_real_is_dense (x1)) || (mm_real_is_dense (x1) && mm_real_is_sparse (x1)))
		error_and_exit ("mm_real_vertcat", "format of matrix x1 and x2 are incompatible.", __FILE__, __LINE__);
	if (mm_real_is_symmetric (x1) || mm_real_is_symmetric (x2))
		error_and_exit ("mm_real_vertcat", "matrix must be general.", __FILE__, __LINE__);
	if (x1->n != x2->n) error_and_exit ("mm_real_vertcat", "matrix size is incompatible.", __FILE__, __LINE__);

	return (mm_real_is_sparse (x1)) ? mm_real_vertcat_sparse (x1, x2) : mm_real_vertcat_dense (x1, x2);
}

/* s = [s1, s2] */
static mm_sparse *
mm_real_horzcat_sparse (mm_sparse *s1, mm_sparse *s2)
{
	MM_INT		j, k;
	MM_INT		m = s1->m;
	MM_INT		n = s1->n + s2->n;
	MM_INT		nnz = s1->nnz + s2->nnz;
	mm_sparse	*s = mm_real_new (MM_REAL_SPARSE, MM_REAL_GENERAL, m, n, nnz);

	MM_INT		*si = s->i;
	MM_INT		*sp = s->p;
	MM_INT		nnz1 = s1->nnz;
	MM_INT		*si1 = s1->i;
	MM_INT		*sp1 = s1->p;
	MM_INT		nnz2 = s2->nnz;
	MM_INT		*si2 = s2->i;
	MM_INT		*sp2 = s2->p;

	for (k = 0; k < nnz1; k++) si[k] = si1[k];
	for (k = 0; k < nnz2; k++) si[k + nnz1] = si2[k];
	for (j = 0; j <= s1->n; j++) sp[j] = sp1[j];
	for (j = 0; j <= s2->n; j++) sp[j + s1->n] = sp2[j] + nnz1;
	if (nnz1 > 0) dcopy_ (&nnz1, s1->data, &ione, s->data, &ione);
	if (nnz2 > 0) dcopy_ (&nnz2, s2->data, &ione, s->data + nnz1, &ione);

	return s;
}

/* d = [d1, d2] */
static mm_dense *
mm_real_horzcat_dense (mm_dense *d1, mm_dense *d2)
{
	mm_dense	*d = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, d1->m, d1->n + d2->n, d1->nnz + d2->nnz);

	dcopy_ (&d1->nnz, d1->data, &ione, d->data, &ione);
	dcopy_ (&d2->nnz, d2->data, &ione, d->data + d1->nnz, &ione);

	return d;
}

/*** x = [x1, x2] ***/
mm_real *
mm_real_horzcat (mm_real *x1, mm_real *x2)
{
	if ((mm_real_is_sparse (x1) && mm_real_is_dense (x1)) || (mm_real_is_dense (x1) && mm_real_is_sparse (x1)))
		error_and_exit ("mm_real_horzcat", "format of matrix x1 and x2 are incompatible.", __FILE__, __LINE__);
	if (mm_real_is_symmetric (x1) || mm_real_is_symmetric (x2))
		error_and_exit ("mm_real_horzcat", "matrix must be general.", __FILE__, __LINE__);
	if (x1->m != x2->m) error_and_exit ("mm_real_horzcat", "matrix size is incompatible.", __FILE__, __LINE__);

	return (mm_real_is_sparse (x1)) ? mm_real_horzcat_sparse (x1, x2) : mm_real_horzcat_dense (x1, x2);
}

/*** x(:,j) += alpha ***/
void
mm_real_xj_add (mm_real *x, MM_INT j, MM_DBL alpha)
{
	MM_INT	k;
	MM_INT	n;
	MM_DBL	*data;
	if (mm_real_is_symmetric (x)) error_and_exit ("mm_real_xj_add_const", "matrix must be general.", __FILE__, __LINE__);
	if (j < 0 || x->n <= j) error_and_exit ("mm_real_xj_add_const", "index out of range.", __FILE__, __LINE__);

	if (mm_real_is_sparse (x)) {
		n = x->p[j + 1] - x->p[j];
		data = x->data + x->p[j];
	} else {
		n = x->m;
		data = x->data + j * x->m;
	}
	for (k = 0; k < n; k++) data[k] += alpha;

	return;
}

/*** x(:,j) *= alpha ***/
void
mm_real_xj_scale (mm_real *x, MM_INT j, MM_DBL alpha)
{
	MM_INT	n;
	MM_DBL	*data;
	if (mm_real_is_symmetric (x)) error_and_exit ("mm_real_xj_scale", "matrix must be general.", __FILE__, __LINE__);
	if (j < 0 || x->n <= j) error_and_exit ("mm_real_xj_scale", "index out of range.", __FILE__, __LINE__);

	if (mm_real_is_sparse (x)) {
		MM_INT	p = x->p[j];
		n = x->p[j + 1] - p;
		data = x->data + p;
	} else {
		n = x->m;
		data = x->data + j * x->m;
	}
	dscal_ (&n, &alpha, data, &ione);

	return;
}

/*** x *= alpha ***/
void
mm_real_scale (mm_real *x, MM_DBL alpha)
{
	MM_INT	n = x->nnz;
	dscal_ (&n, &alpha, x->data, &ione);
	return;
}

/* sum |s(:,j)| */
static MM_DBL
mm_real_sj_asum (mm_sparse *s, MM_INT j)
{
	MM_INT	p = s->p[j];
	MM_INT	n = s->p[j + 1] - p;
	MM_DBL	*sd = s->data;
	MM_DBL	asum = dasum_ (&n, sd + p, &ione);
	if (mm_real_is_symmetric (s)) {
		MM_INT	k;
		MM_INT	k0;
		MM_INT	k1;
		if (mm_real_is_upper (s)) {
			k0 = j + 1;
			k1 = s->n;
		} else {
			k0 = 0;
			k1 = j;
		}
		for (k = k0; k < k1; k++) {
			MM_INT	l = find_row_element (j, s, k);
			// if found
			if (l >= 0) asum += fabs (sd[l]);
		}
	}
	return asum;
}

/* sum |d(:,j)| */
static MM_DBL
mm_real_dj_asum (mm_dense *d, MM_INT j)
{
	MM_DBL	val = 0.;
	if (!mm_real_is_symmetric (d)) val = dasum_ (&d->m, d->data + j * d->m, &ione);
	else {
		MM_INT	n;
		if (mm_real_is_upper (d)) {
			n = j;
			val = dasum_ (&n, d->data + j * d->m, &ione);
			n = d->m - j;
			val += dasum_ (&n, d->data + j * d->m + j, &d->m);
		} else if (mm_real_is_lower (d)) {
			n = d->m - j;
			val = dasum_ (&n, d->data + j * d->m + j, &ione);
			n = j;
			val += dasum_ (&n, d->data + j, &d->m);
		}
	}
	return val;
}

/*** sum |x(:,j)| ***/
MM_DBL
mm_real_xj_asum (mm_real *x, MM_INT j)
{
	if (j < 0 || x->n <= j) error_and_exit ("mm_real_xj_asum", "index out of range.", __FILE__, __LINE__);
	return (mm_real_is_sparse (x)) ? mm_real_sj_asum (x, j) : mm_real_dj_asum (x, j);
}

/* sum s(:,j) */
static MM_DBL
mm_real_sj_sum (mm_sparse *s, MM_INT j)
{
	MM_INT	k;
	MM_INT	p = s->p[j];
	MM_INT	n = s->p[j + 1] - p;
	MM_DBL	*sd = s->data + p;
	MM_DBL	sum = 0.;
	for (k = 0; k < n; k++) sum += sd[k];
	if (mm_real_is_symmetric (s)) {
		MM_INT	k0;
		MM_INT	k1;
		if (mm_real_is_upper (s)) {
			k0 = j + 1;
			k1 = s->n;
		} else {
			k0 = 0;
			k1 = j;
		}
		sd = s->data;
		for (k = k0; k < k1; k++) {
			MM_INT	l = find_row_element (j, s, k);
			// if found
			if (l >= 0) sum += sd[l];
		}
	}
	return sum;
}

/* sum d(:,j) */
static MM_DBL
mm_real_dj_sum (mm_dense *d, MM_INT j)
{
	MM_INT	k;
	MM_DBL	*dd;
	MM_DBL	sum = 0.;
	if (!mm_real_is_symmetric (d)) {
		dd = d->data + j * d->m;
		for (k = 0; k < d->m; k++) sum += dd[k];
	} else {
		MM_INT	n;
		if (mm_real_is_upper (d)) {
			n = j;
			dd = d->data + j * d->m;
			for (k = 0; k < n; k++) sum += dd[k];
			n = d->m - j;
			dd = d->data + j + j * d->m;
			for (k = 0; k < n; k++) sum += dd[k * d->m];
		} else if (mm_real_is_lower (d)) {
			n = d->m - j;
			dd = d->data + j + j * d->m;
			for (k = 0; k < n; k++) sum += dd[k];
			n = j;
			dd = d->data + j;
			for (k = 0; k < n; k++) sum += dd[k * d->m];
		}
	}
	return sum;
}

/*** sum x(:,j) ***/
MM_DBL
mm_real_xj_sum (mm_real *x, MM_INT j)
{
	if (j < 0 || x->n <= j) error_and_exit ("mm_real_xj_sum", "index out of range.", __FILE__, __LINE__);
	return (mm_real_is_sparse (x)) ? mm_real_sj_sum (x, j) : mm_real_dj_sum (x, j);
}

/* sum_i s(i,j)^2 */
static MM_DBL
mm_real_sj_ssq (mm_sparse *s, MM_INT j)
{
	MM_INT	p = s->p[j];
	MM_INT	n = s->p[j + 1] - p;
	MM_DBL	*sd = s->data;
	MM_DBL	ssq = ddot_ (&n, sd + p, &ione, sd + p, &ione);
	if (mm_real_is_symmetric (s)) {
		MM_INT	k;
		MM_INT	k0;
		MM_INT	k1;
		if (mm_real_is_upper (s)) {
			k0 = j + 1;
			k1 = s->n;
		} else {
			k0 = 0;
			k1 = j;
		}
		for (k = k0; k < k1; k++) {
			MM_INT		l = find_row_element (j, s, k);
			// if found
			if (l >= 0) ssq += pow (sd[l], 2.);
		}
	}
	return ssq;
}

/* sum_i d(i,j)^2 */
static MM_DBL
mm_real_dj_ssq (mm_dense *d, MM_INT j)
{
	MM_DBL	ssq;
	if (!mm_real_is_symmetric (d)) ssq = ddot_ (&d->m, d->data + j * d->m, &ione, d->data + j * d->m, &ione);
	else {
		MM_INT	n;
		ssq = 0.;
		if (mm_real_is_upper (d)) {
			n = j;
			ssq = ddot_ (&n, d->data + j * d->m, &ione, d->data + j * d->m, &ione);
			n = d->m - j;
			ssq += ddot_ (&n, d->data + j * d->m + j, &d->m, d->data + j * d->m + j, &d->m);
		} else if (mm_real_is_lower (d)) {
			n = d->m - j;
			ssq = ddot_ (&n, d->data + j * d->m + j, &ione, d->data + j * d->m + j, &ione);
			n = j;
			ssq += ddot_ (&n, d->data + j, &d->m, d->data + j, &d->m);
		}
	}
	return ssq;
}

/*** sum_i x(i,j)^2 ***/
MM_DBL
mm_real_xj_ssq (mm_real *x, MM_INT j)
{
	if (j < 0 || x->n <= j) error_and_exit ("mm_real_xj_ssq", "index out of range.", __FILE__, __LINE__);
	return (mm_real_is_sparse (x)) ? mm_real_sj_ssq (x, j) : mm_real_dj_ssq (x, j);
}

/*** norm2 x(:,j) ***/
MM_DBL
mm_real_xj_nrm2 (mm_real *x, MM_INT j)
{
	MM_DBL	ssq;
	if (j < 0 || x->n <= j) error_and_exit ("mm_real_xj_nrm2", "index out of range.", __FILE__, __LINE__);
	ssq = mm_real_xj_ssq (x, j);
	return sqrt (ssq);
}

#ifdef DEBUG
static void
mm_real_s_dot_dk1 (bool trans, MM_DBL alpha, mm_sparse *s, mm_dense *d, MM_INT k, MM_DBL beta, mm_dense *z, MM_INT q)
{
	MM_INT	i, j;
	if (trans) {
//#pragma omp parallel for
		for (j = 0; j < s->n; j++) {
			mm_real	*sj = mm_real_sj_col (s, j);
			z->data[j + q * z->m] = alpha * ddot_ (&sj->m, sj->data, &ione, d->data + k * d->m, &ione) + beta * z->data[j + q * z->m];
			mm_real_free (sj);
		}
	} else {
//#pragma omp parallel for
		for (i = 0; i < s->m; i++) {
			mm_real	*si = mm_real_si_row (s, i);
			z->data[i + q * z->m] = alpha * ddot_ (&si->n, si->data, &ione, d->data + k * d->m, &ione) + beta * z->data[i + q * z->m];
			mm_real_free (si);
		}
	}
	return;
}
#endif

/* z = alpha * s * dy(:,k) + beta * z, where s is sparse matrix and y is dense general */
static void
mm_real_s_dot_dk (bool trans, MM_DBL alpha, mm_sparse *s, mm_dense *y, MM_INT k, MM_DBL beta, mm_dense *z, MM_INT q)
{
	MM_INT		j, l;
	MM_INT		p, m, n;

	MM_INT		*si = s->i;
	MM_DBL		*sd = s->data;
	MM_INT		*sp = s->p;
	MM_DBL		*yk = y->data + k * y->m;
	MM_DBL		*zl = z->data + q * z->m;
	MM_DBL		*zl0 = zl;

	m = (!trans) ? s->m : s->n;
	if (fabs (beta) > __DBL_EPSILON__) dscal_ (&m, &beta, zl, &ione);
	else mm_real_array_set_all (m, zl, 0.);
	if (trans) { // trans s
		if (!mm_real_is_symmetric (s)) { // s is not symmetric
			for (j = 0; j < s->n; j++, zl++) {
				//p = *(sp + j);
				//n = *(sp + j + 1) - p;
				p = *sp;
				n = *(++sp) - p;
				if (n <= 0) continue;
				//si = s->i + p;
				//sd = s->data + p;
				for (l = 0; l < n; l++, si++, sd++) *zl += alpha * (*sd) * (*(yk + *si));
			}
		} else if (mm_real_is_upper (s)) { // s is symmetric upper
			MM_INT	sil;
			MM_DBL	alpha_sdk;
			MM_DBL	*ykj = yk;
			for (j = 0; j < s->n; j++, zl++, ykj++) {
				p = *sp;
				n = *(++sp) - p;
				if (n <= 0) continue;
				//si = s->i + p;
				//sd = s->data + p;
				for (l = 0; l < n - 1; l++, si++, sd++) {
					sil = *si;
					alpha_sdk = alpha * (*sd);
					*zl += alpha_sdk * (*(yk + sil));
					*(zl0 + sil) += alpha_sdk * (*ykj);
				}
				sil = *si;
				alpha_sdk = alpha * (*sd);
				*zl += alpha_sdk * (*(yk + sil));
				// omit if diagonal element
				if (j != sil) *(zl0 + sil) += alpha_sdk * (*ykj);
				//////
				si++;
				sd++;
			}
		} else { // s is symmetric lower
			MM_INT	sil;
			MM_DBL	alpha_sdl;
			MM_DBL	*ykj = yk;
			for (j = 0; j < s->n; j++, zl++, ykj++) {
				p = *sp;
				n = *(++sp) - p;
				if (n <= 0) continue;
				//si = s->i + p;
				//sd = s->data + p;
				sil = *si;
				alpha_sdl = alpha * (*sd);
				*zl += alpha_sdl * (*(yk + sil));
				// omit if diagonal element
				if (j != sil) *(zl0 + sil) += alpha_sdl * (*ykj);
				si++;
				sd++;
				for (l = 1; l < n; l++, si++, sd++) {
					sil = *si;
					alpha_sdl = alpha * (*sd);
					*zl += alpha_sdl * (*(yk + sil));
					*(zl0 + sil) += alpha_sdl * (*ykj);
				}
			}
		}
	} else { // no trans s
		if (!mm_real_is_symmetric (s)) { // s is not symmetric
			MM_DBL	*ykj = yk;
			for (j = 0; j < s->n; j++, ykj++) {
				p = *sp;
				n = *(++sp) - p;
				if (n <= 0) continue;
				//si = s->i + p;
				//sd = s->data + p;
				for (l = 0; l < n; l++, si++, sd++) *(zl + *si) += alpha * (*sd) * (*ykj);
			}
		} else if (mm_real_is_upper (s)) { // s is symmetric upper
			MM_INT	sil;
			MM_DBL	alpha_sdl;
			MM_DBL	*ykj = yk;
			for (j = 0; j < s->n; j++, zl++, ykj++) {
				p = *sp;
				n = *(++sp) - p;
				if (n <= 0) continue;
				//si = s->i + p;
				//sd = s->data + p;
				for (l = 0; l < n - 1; l++, si++, sd++) {
					sil = *si;
					alpha_sdl = alpha * (*sd);
					*(zl0 + sil) += alpha_sdl * (*ykj);
					*zl += alpha_sdl * (*(yk + sil));
				}
				sil = *si;
				alpha_sdl = alpha * (*sd);
				*(zl0 + sil) += alpha_sdl * (*ykj);
				// omit if diagonal element
				if (j != sil) *zl += alpha_sdl * (*(yk + sil));
				//////
				si++;
				sd++;
			}
		} else { // s is symmetric lower
			MM_INT	sil;
			MM_DBL	alpha_sdl;
			MM_DBL	*ykj = yk;
			for (j = 0; j < s->n; j++, zl++, ykj++) {
				p = *sp;
				n = *(++sp) - p;
				if (n <= 0) continue;
				//si = s->i + p;
				//sd = s->data + p;
				sil = *si;
				alpha_sdl = alpha * (*sd);
				*(zl0 + sil) += alpha_sdl * (*ykj);
				// omit if diagonal element
				if (j != sil) *zl += alpha_sdl * (*(yk + sil));
				si = si + 1;
				sd = sd + 1;
				for (l = 1; l < n; l++, si++, sd++) {
					sil = *si;
					alpha_sdl = alpha * (*sd);
					*(zl0 + sil) += alpha_sdl * (*ykj);
					*zl += alpha_sdl * (*(yk + sil));
				}
			}
		}
	}
	return;
}

/* z = alpha * d * dy(:,k) + beta * z, where d is dense matrix and y is dense general */
static void
mm_real_d_dot_dk (bool trans, MM_DBL alpha, mm_dense *d, mm_dense *y, MM_INT k, MM_DBL beta, mm_dense *z, MM_INT l)
{
	MM_DBL	*yk = y->data + k * y->m;
	MM_DBL	*zl = z->data + l * z->m;
	char	_trans = 'T';
	char	_notrans = 'N';
	if (!mm_real_is_symmetric (d)) {
		// z = alpha * d * y + beta * z
		dgemv_ ((trans) ? &_trans : &_notrans, &d->m, &d->n, &alpha, d->data, &d->m, yk, &ione, &beta, zl, &ione);
	} else {
		char	uplo = (mm_real_is_upper (d)) ? 'U' : 'L';
		// z = alpha * d * y + beta * z
		dsymv_ (&uplo, &d->m, &alpha, d->data, &d->m, yk, &ione, &beta, zl, &ione);
	}
	return;
}

/* z = alpha * s * sy(:,k) + beta * z, where s is sparse matrix and y is sparse general */
static void
mm_real_s_dot_sk (bool trans, MM_DBL alpha, mm_sparse *s, mm_sparse *y, MM_INT k, MM_DBL beta, mm_dense *z, MM_INT l)
{
	mm_real	*sk = mm_real_sj_col (y, k);
	mm_real_s_dot_dk (trans, alpha, s, sk, 0, beta, z, l);
	mm_real_free (sk);
	return;
}

/* z = alpha * d * sy(:,k) + beta * z, where d is dense matrix and y is sparse general */
static void
mm_real_d_dot_sk (bool trans, MM_DBL alpha, mm_dense *d, mm_sparse *y, MM_INT k, MM_DBL beta, mm_dense *z, MM_INT l)
{
	mm_real	*sk = mm_real_sj_col (y, k);
	mm_real_d_dot_dk (trans, alpha, d, sk, 0, beta, z, l);
	mm_real_free (sk);
	return;
}

/*** z(:,k) = alpha * x * y(:,k), where x and y are sparse/dense matrix ***/
void
mm_real_x_dot_yk (bool trans, MM_DBL alpha, mm_real *x, mm_real *y, MM_INT k, MM_DBL beta, mm_dense *z)
{
	if (y->n < k) error_and_exit ("mm_real_x_dot_yk", "k exceeds num of col of y.", __FILE__, __LINE__);
	if (z->n < k) error_and_exit ("mm_real_x_dot_yk", "k exceeds num of col of z.", __FILE__, __LINE__);
	if ((trans && x->m != y->m) || (!trans && x->n != y->m))
		error_and_exit ("mm_real_x_dot_yk", "dimensions of x and y do not match.", __FILE__, __LINE__);
	if ((trans && x->n != z->m) || (!trans && x->m != z->m))
		error_and_exit ("mm_real_x_dot_yk", "dimensions of x and z do not match.", __FILE__, __LINE__);
	if (mm_real_is_symmetric (y)) error_and_exit ("mm_real_x_dot_yk", "y must be general.", __FILE__, __LINE__);
	if (mm_real_is_sparse (x)) {
		// x is sparse
		(mm_real_is_sparse (y)) ?\
			mm_real_s_dot_sk (trans, alpha, x, y, k, beta, z, k) : mm_real_s_dot_dk (trans, alpha, x, y, k, beta, z, k);
	} else {
		// x is dense
		(mm_real_is_sparse (y)) ?\
			mm_real_d_dot_sk (trans, alpha, x, y, k, beta, z, k) : mm_real_d_dot_dk (trans, alpha, x, y, k, beta, z, k);
	}
	return;
}

void
mm_real_x_dot_y (bool transx, bool transy, double alpha, mm_real *x, mm_real *y, double beta, mm_real *z)
{
	if (mm_real_is_symmetric (y)) error_and_exit ("mm_real_x_dot_y", "y must be general.", __FILE__, __LINE__);
	{
		MM_INT	mx = (transx) ? x->n : x->m;
		MM_INT	nx = (transx) ? x->m : x->n;
		MM_INT	my = (transy) ? y->n : y->m;
		MM_INT	ny = (transy) ? y->m : y->n;
		if (nx != my) error_and_exit ("mm_real_x_dot_y", "dimensions of x and y do not match.", __FILE__, __LINE__);
		if (z->m != mx) error_and_exit ("mm_real_x_dot_y", "dimensions of x and z do not match.", __FILE__, __LINE__);
		if (z->n != ny) error_and_exit ("mm_real_x_dot_y", "dimensions of y and z do not match.", __FILE__, __LINE__);
	}	
	if (!mm_real_is_dense (z)) error_and_exit ("mm_real_x_dot_y", "z must be dense.", __FILE__, __LINE__);
	if (mm_real_is_symmetric (z)) error_and_exit ("mm_real_x_dot_y", "z must be general.", __FILE__, __LINE__);

	if (mm_real_is_dense (x)) { // x is dense
		char	tx = (transx) ? 'T' : 'N';
		char	ty = (transy) ? 'T' : 'N';
		if (mm_real_is_dense (y)) { // y is dense
			// x is dense, y is dense
			MM_INT	m;
			MM_INT	n;
			MM_INT	k;
			MM_INT	lda;
			MM_INT	ldb;
			if (!mm_real_is_symmetric (x)) { // matrix x is asymmetric
				m = (transx) ? x->n : x->m;
				n = (transy) ? y->m : y->n;
				k = (transx) ? x->m : x->n;
				lda = x->m;
				ldb = (transy) ? n : k;
				dgemm_ (&tx, &ty, &m, &n, &k, &alpha, x->data, &lda, y->data, &ldb, &beta, z->data, &z->m);
			} else { // matrix x is symmetric
				char	side = 'L';
				char	uplo = (mm_real_is_upper (x)) ? 'U' : 'L';
				if (!transy) {
					// x * y
					dsymm_ (&side, &uplo, &z->m, &z->n, &alpha, x->data, &x->m, y->data, &y->m, &beta, z->data, &z->m);
				} else {
					// x * y'
					MM_INT	i;
#pragma omp parallel for
					for (i = 0; i < y->m; i++) {
						double	*d = z->data + i * z->m;
						mm_real	*c = mm_real_xi_row (y, i);
						//mm_real_transpose (c);
						// d = alpha * x * c + beta * d
						dsymv_ (&uplo, &x->m, &alpha, x->data, &x->m, c->data, &ione, &beta, d, &ione);
						mm_real_free (c);
					}
				}
			}
		} else { // y is sparse
			// x is dense, y is sparse
			if (transy) {
#pragma omp parallel for
				for (MM_INT j = 0; j < y->m; j++) {
					mm_real	*c = mm_real_xi_row (y, j);
					double	*d = z->data + j * z->m;
					if (!mm_real_is_symmetric (x)) {
						// d = alpha * x * c + beta * d
						dgemv_ (&tx, &x->m, &x->n, &alpha, x->data, &x->m, c->data, &ione, &beta, d, &ione);
					} else {
						char	uplo = (mm_real_is_upper (x)) ? 'U' : 'L';
						// d = alpha * x * c + beta * d
						dsymv_ (&uplo, &x->m, &alpha, x->data, &x->m, c->data, &ione, &beta, d, &ione);
					}
					mm_real_free (c);
				}
			} else {
#pragma omp parallel for
				for (MM_INT k = 0; k < y->n; k++) mm_real_x_dot_yk (transx, alpha, x, y, k, beta, z);
			}
		}
	} else { // x is sparse
		// x is sparse, y is dense / sparse
		if (transy) {
#pragma omp parallel for
			for (MM_INT j = 0; j < y->m; j++) {
				mm_real	*c = mm_real_xi_row (y, j);
				mm_real_transpose (c);
				mm_real	*d = mm_real_view_array (MM_REAL_DENSE, MM_REAL_GENERAL, z->m, 1, z->m, z->data + j * z->m);
				mm_real_x_dot_yk (transx, alpha, x, c, 0, beta, d);
				mm_real_free (d);
				mm_real_free (c);
			}
		} else {
#pragma omp parallel for
			for (MM_INT k = 0; k < y->n; k++) mm_real_x_dot_yk (transx, alpha, x, y, k, beta, z);
		}
	}
	return;
}


/*****************************************************************************/

/* s(:,j)' * d(:,k) */
static MM_DBL
mm_real_sj_trans_dot_dk (mm_sparse *s, MM_INT j, mm_dense *y, MM_INT k)
{
	MM_DBL	val = 0;

	MM_INT	*si;
	MM_DBL	*sd;
	MM_DBL	*yk = y->data + k * y->m;

	MM_INT	l;
	MM_INT	p = s->p[j];
	MM_INT	n = s->p[j + 1] - p;

	si = s->i + p;
	sd = s->data + p;
	for (l = 0; l < n; l++) val += sd[l] * yk[si[l]];

	if (mm_real_is_symmetric (s)) {
		MM_INT	l0;
		MM_INT	l1;
		if (mm_real_is_upper (s)) {
			l0 = j + 1;
			l1 = s->n;
		} else {
			l0 = 0;
			l1 = j;
		}
		sd = s->data;
		for (l = l0; l < l1; l++) {
			MM_INT	t = find_row_element (j, s, l);
			// if found
			if (t >= 0) val += sd[t] * yk[l];
		}
	}

	return val;
}

/* d(:,j)' * d(:,k) */
static MM_DBL
mm_real_dj_trans_dot_dk (mm_dense *d, MM_INT j, mm_dense *y, MM_INT k)
{
	MM_DBL	val = 0.;
	MM_DBL	*yk = y->data + k * y->m;
	if (!mm_real_is_symmetric (d)) val = ddot_ (&d->m, d->data + j * d->m, &ione, yk, &ione);
	else {
		MM_INT	n;
		if (mm_real_is_upper (d)) {
			n = j;
			val = ddot_ (&n, d->data + j * d->m, &ione, yk, &ione);
			n = d->m - j;
			val += ddot_ (&n, d->data + j * d->m + j, &d->m, yk + j, &ione);
		} else if (mm_real_is_lower (d)) {
			n = d->m - j;
			val = ddot_ (&n, d->data + j * d->m + j, &ione, yk + j, &ione);
			n = j;
			val += ddot_ (&n, d->data + j, &d->m, yk, &ione);
		}
	}
	return val;
}

/* s(:,j)' * s(:,k) */
static MM_DBL
mm_real_sj_trans_dot_sk (mm_sparse *s, MM_INT j, mm_sparse *y, MM_INT k)
{
	mm_real	*yk = mm_real_sj_col (y, k);
	MM_DBL	val = mm_real_sj_trans_dot_dk (s, j, yk, 0);
	mm_real_free (yk);
	return val;
}

/* d(:,j)' * s(:,k) */
static MM_DBL
mm_real_dj_trans_dot_sk (mm_dense *d, MM_INT j, mm_sparse *y, MM_INT k)
{
	mm_real	*yk = mm_real_sj_col (y, k);
	MM_DBL	val = mm_real_dj_trans_dot_dk (d, j, yk, 0);
	mm_real_free (yk);
	return val;
}

/*****************************************************************************/

/*** x(:,j)' * y(:,k) ***/
MM_DBL
mm_real_xj_trans_dot_yk (mm_real *x, MM_INT j, mm_dense *y, MM_INT k)
{
	MM_DBL	val = 0.;
	if (j < 0 || x->n <= j) error_and_exit ("mm_real_xj_trans_dot_yk", "first index out of range.", __FILE__, __LINE__);
	if (k < 0 || y->n <= k) error_and_exit ("mm_real_xj_trans_dot_yk", "second index out of range.", __FILE__, __LINE__);
	if (mm_real_is_symmetric (y)) error_and_exit ("mm_real_xj_trans_dot_yk", "y must be general.", __FILE__, __LINE__);
	if (x->m != y->m) error_and_exit ("mm_real_xj_trans_dot_yk", "matrix dimensions do not match.", __FILE__, __LINE__);

	if (mm_real_is_sparse (x)) { // x is sparse
		val = (mm_real_is_sparse (y)) ? mm_real_sj_trans_dot_sk (x, j, y, k) : mm_real_sj_trans_dot_dk (x, j, y, k);
	} else { // x is dense
		val = (mm_real_is_sparse (y)) ? mm_real_dj_trans_dot_sk (x, j, y, k) : mm_real_dj_trans_dot_dk (x, j, y, k);
	}
	return val;
}

/*** x(:,j)' * y ***/
mm_dense *
mm_real_xj_trans_dot_y (mm_real *x, MM_INT j, mm_dense *y)
{
	MM_INT		k;
	mm_dense	*z;
	if (j < 0 || x->n <= j) error_and_exit ("mm_real_xj_trans_dot_y", "first index out of range.", __FILE__, __LINE__);
	if (mm_real_is_symmetric (y)) error_and_exit ("mm_real_xj_trans_dot_y", "y must be general.", __FILE__, __LINE__);
	if (x->m != y->m) error_and_exit ("mm_real_xj_trans_dot_y", "matrix dimensions do not match.", __FILE__, __LINE__);
	z = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, 1, y->n, y->n);
#pragma omp parallel for
	for (k = 0; k < y->n; k++) z->data[k] = mm_real_xj_trans_dot_yk (x, j, y, k);
	return z;
}

/* y = alpha * s(:,j) + y */
static void
mm_real_asjpy (MM_DBL alpha, mm_sparse *s, MM_INT j, mm_dense *y)
{
	MM_INT	*si;
	MM_DBL	*sd;
	MM_DBL	*yd = y->data;

	MM_INT	k;
	MM_INT	p = s->p[j];
	MM_INT	n = s->p[j + 1] - p;
	si = s->i + p;
	sd = s->data + p;
	for (k = 0; k < n; k++) yd[si[k]] += alpha * sd[k];
	if (mm_real_is_symmetric (s)) {
		MM_INT	k0;
		MM_INT	k1;
		if (mm_real_is_upper (s)) {
			k0 = j + 1;
			k1 = s->n;
		} else {
			k0 = 0;
			k1 = j;
		}
		sd = s->data;
		for (k = k0; k < k1; k++) {
			MM_INT	l = find_row_element (j, s, k);
			// if found
			if (l >= 0) yd[k] += alpha * sd[l];
		}
	}
	return;
}

/* y = alpha * d(:,j) + y */
static void
mm_real_adjpy (MM_DBL alpha, mm_dense *d, MM_INT j, mm_dense *y)
{
	if (!mm_real_is_symmetric (d)) daxpy_ (&d->m, &alpha, d->data + j * d->m, &ione, y->data, &ione);
	else {
		MM_INT	n;
		if (mm_real_is_upper (d)) {
			n = j;
			daxpy_ (&n, &alpha, d->data + j * d->m, &ione, y->data, &ione);
			n = d->m - j;
			daxpy_ (&n, &alpha, d->data + j * d->m + j, &d->m, y->data + j, &ione);
		} else if (mm_real_is_lower (d)) {
			n = d->m - j;
			daxpy_ (&n, &alpha, d->data + j * d->m + j, &ione, y->data + j, &ione);
			n = j;
			daxpy_ (&n, &alpha, d->data + j, &d->m, y->data, &ione);
		}
	}
	return;
}

/*** y = alpha * x(:,j) + y ***/
void
mm_real_axjpy (MM_DBL alpha, mm_real *x, MM_INT j, mm_dense *y)
{
	if (j < 0 || x->n <= j) error_and_exit ("mm_real_axjpy", "index out of range.", __FILE__, __LINE__);
	if (!mm_real_is_dense (y)) error_and_exit ("mm_real_axjpy", "y must be dense.", __FILE__, __LINE__);
	if (mm_real_is_symmetric (y)) error_and_exit ("mm_real_axjpy", "y must be general.", __FILE__, __LINE__);
	if (y->n != 1) error_and_exit ("mm_real_axjpy", "y must be vector.", __FILE__, __LINE__);
	if (x->m != y->m) error_and_exit ("mm_real_axjpy", "vector and matrix dimensions do not match.", __FILE__, __LINE__);

	return (mm_real_is_sparse (x)) ? mm_real_asjpy (alpha, x, j, y) : mm_real_adjpy (alpha, x, j, y);
}

/* fread sparse */
static mm_sparse *
mm_real_fread_sparse (FILE *fp, MM_typecode typecode)
{
	MM_INT		k, l;
	MM_INT		m, n;
	MM_INT		nnz;
	MM_INT		*j;
	mm_sparse	*s;

	if (mm_read_mtx_crd_size (fp, &m, &n, &nnz) != 0) return NULL;
	s = mm_real_new (MM_REAL_SPARSE, MM_REAL_GENERAL, m, n, nnz);

	j = (MM_INT *) malloc (s->nnz * sizeof (MM_INT));
	if (mm_read_mtx_crd_data (fp, s->m, s->n, s->nnz, s->i, j, s->data, typecode) != 0) {
		free (j);
		mm_real_free (s);
		return NULL;
	}

	l = 0;
	for (k = 0; k < nnz; k++) {
		s->i[k]--;	// fortran -> c
		while (l < j[k]) s->p[l++] = k;
	}
	while (l <= n) s->p[l++] = k;

	if (mm_is_symmetric (typecode)) {
		mm_real_set_symmetric (s);
		for (k = 0; k < nnz; k++) {
			if (s->i[k] == j[k] - 1) continue;
			(s->i[k] < j[k] - 1) ? mm_real_set_upper (s) : mm_real_set_lower (s);
			break;
		}
	}
	free (j);

	return s;
}

/* fread dense */
static mm_dense *
mm_real_fread_dense (FILE *fp, MM_typecode typecode)
{
	MM_INT		k;
	MM_INT		m, n;
	MM_INT		ret;
	mm_dense	*d;

	if (mm_read_mtx_array_size (fp, &m, &n) != 0) return NULL;
	d = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, n, m * n);

	k = 0;
	do {
		ret = fscanf (fp, "%lf", &d->data[k]);
		if (ret > 0 && ++k >= d->nnz) break;
	} while (ret != EOF);

	return d;
}

/*** fread MatrixMarket format file ***/
mm_real *
mm_real_fread (FILE *fp)
{
	MM_typecode	typecode;
	mm_real		*x;
	if (mm_read_banner (fp, &typecode) != 0) error_and_exit ("mm_real_fread", "failed to read mm_real.", __FILE__, __LINE__);
	if (!is_type_supported (typecode)) {
		char	msg[128];
		sprintf (msg, "matrix type does not supported :[%s].", mm_typecode_to_str (typecode));
		error_and_exit ("mm_real_fread", msg, __FILE__, __LINE__);
	}
	x = (mm_is_sparse (typecode)) ? mm_real_fread_sparse (fp, typecode) : mm_real_fread_dense (fp, typecode);
	if (!x) error_and_exit ("mm_real_fread", "failed to read mm_real.", __FILE__, __LINE__);
	if (mm_real_is_symmetric (x) && x->m != x->n) error_and_exit ("mm_real_fread", "symmetric matrix must be square.", __FILE__, __LINE__);
	return x;
}

/* fwrite sparse */
static void
mm_real_fwrite_sparse (FILE *stream, mm_sparse *s, const char *format)
{
	MM_INT		j;
	mm_write_banner (stream, s->typecode);
	mm_write_mtx_crd_size (stream, s->m, s->n, s->nnz);
	for (j = 0; j < s->n; j++) {
		MM_INT	k = s->p[j];
		MM_INT	pend = s->p[j + 1];
		for (; k < pend; k++) {
#ifdef MMINT
			fprintf (stream, "%d %d ", s->i[k] + 1, j + 1);	// c -> fortran
#else
			fprintf (stream, "%ld %ld ", s->i[k] + 1, j + 1);	// c -> fortran
#endif
			fprintf (stream, format, s->data[k]);
			fprintf (stream, "\n");
		}
	}
	return;
}

/* fwrite dense */
static void
mm_real_fwrite_dense (FILE *stream, mm_dense *d, const char *format)
{
	MM_INT		k;
	mm_write_banner (stream, d->typecode);
	mm_write_mtx_array_size (stream, d->m, d->n);
	for (k = 0; k < d->nnz; k++) {
		fprintf (stream, format, d->data[k]);
		fprintf (stream, "\n");
	}
	return;
}

/*** fwrite in MatrixMarket format ***/
void
mm_real_fwrite (FILE *stream, mm_real *x, const char *format)
{
	return (mm_real_is_sparse (x)) ? mm_real_fwrite_sparse (stream, x, format) : mm_real_fwrite_dense (stream, x, format);
}

static void
check_fread_status (int stat, const char *fname, const char *file, int line)
{
	if (stat < 1) error_and_exit (fname, "fread failed", fname, line);
	return;
}

/* fread binary sparse */
static mm_dense *
mm_real_fread_binary_sparse (FILE *fp)
{
	int		ret;
	MM_INT	m, n;
	MM_INT	nnz;
	mm_real	*s;
	ret = fread (&m, sizeof (MM_INT), 1, fp);
	check_fread_status (ret, "mm_real_fread_binary_sparse", __FILE__, __LINE__);
	ret = fread (&n, sizeof (MM_INT), 1, fp);
	check_fread_status (ret, "mm_real_fread_binary_sparse", __FILE__, __LINE__);
	ret = fread (&nnz, sizeof (MM_INT), 1, fp);
	check_fread_status (ret, "mm_real_fread_binary_sparse", __FILE__, __LINE__);

	s = mm_real_new (MM_REAL_SPARSE, MM_REAL_GENERAL, m, n, nnz);
	for (MM_INT i = 0; i < nnz; i++) {
		ret = fread (&s->i[i], sizeof (MM_INT), 1, fp);
		check_fread_status (ret, "mm_real_fread_binary_sparse", __FILE__, __LINE__);
	}
	for (MM_INT i = 0; i < n + 1; i++) {
		ret = fread (&s->p[i], sizeof (MM_INT), 1, fp);
		check_fread_status (ret, "mm_real_fread_binary_sparse", __FILE__, __LINE__);
	}
	for (MM_INT i = 0; i < nnz; i++) {
		ret = fread (&s->data[i], sizeof (MM_DBL), 1, fp);
		check_fread_status (ret, "mm_real_fread_binary_sparse", __FILE__, __LINE__);
	}
	return s;
}

/* fread binary dense */
static mm_dense *
mm_real_fread_binary_dense (FILE *fp)
{
	int		ret;
	MM_INT	m, n;
	MM_INT	nnz;
	mm_real	*d;
	ret = fread (&m, sizeof (MM_INT), 1, fp);
	check_fread_status (ret, "mm_real_fread_binary_dense", __FILE__, __LINE__);
	ret = fread (&n, sizeof (MM_INT), 1, fp);
	check_fread_status (ret, "mm_real_fread_binary_dense", __FILE__, __LINE__);
	ret = fread (&nnz, sizeof (MM_INT), 1, fp);
	check_fread_status (ret, "mm_real_fread_binary_dense", __FILE__, __LINE__);
	d = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, n, nnz);
	ret = fread (d->data, sizeof (MM_DBL), nnz, fp);
	check_fread_status (ret, "mm_real_fread_binary_dense", __FILE__, __LINE__);
	return d;
}

/*** fread binary ***/
mm_real *
mm_real_fread_binary (FILE *fp)
{
	int		ret;
	char	typecode[4];

	if (!fp) error_and_exit ("mm_real_fread_binary", "file not open.", __FILE__, __LINE__);

	ret = fread (typecode, sizeof (char), 4, fp);
	check_fread_status (ret, "mm_real_fread_binary", __FILE__, __LINE__);
	return (mm_is_dense (typecode)) ? mm_real_fread_binary_dense (fp) : mm_real_fread_binary_sparse (fp);
}

/* fwrite binary sparse */
static void
mm_real_fwrite_binary_sparse (FILE *fp, mm_real *s)
{
	fwrite (s->typecode, sizeof (char), 4, fp); 
	fwrite (&s->m, sizeof (MM_INT), 1, fp);
	fwrite (&s->n, sizeof (MM_INT), 1, fp);
	fwrite (&s->nnz, sizeof (MM_INT), 1, fp);
	for (MM_INT i = 0; i < s->nnz; i++) fwrite (&s->i[i], sizeof (MM_INT), 1, fp);
	for (MM_INT i = 0; i < s->n + 1; i++) fwrite (&s->p[i], sizeof (MM_INT), 1, fp);
	for (MM_INT i = 0; i < s->nnz; i++) fwrite (&s->data[i], sizeof (MM_DBL), 1, fp);
}

/* fwrite binary dense */
static void
mm_real_fwrite_binary_dense (FILE *fp, mm_real *d)
{
	fwrite (d->typecode, sizeof (char), 4, fp); 
	fwrite (&d->m, sizeof (MM_INT), 1, fp);
	fwrite (&d->n, sizeof (MM_INT), 1, fp);
	fwrite (&d->nnz, sizeof (MM_INT), 1, fp);
	fwrite (d->data, sizeof (MM_DBL), d->nnz, fp);
}

/*** fwrite binary ***/
void
mm_real_fwrite_binary (FILE *fp, mm_real *x)
{
	if (!fp) error_and_exit ("mm_real_fwrite_binary", "file not open.", __FILE__, __LINE__);
	return (mm_real_is_dense (x)) ? mm_real_fwrite_binary_dense (fp, x) : mm_real_fwrite_binary_sparse (fp, x);
}

