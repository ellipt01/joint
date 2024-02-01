/*
 * mm_real.h
 *
 *  Created on: 2014/06/20
 *      Author: utsugi
 */

#ifndef MMREAL_H_
#define MMREAL_H_

#include <mmio.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

/* dense / sparse */
typedef enum {
	MM_REAL_SPARSE = 0,	// sparse matrix
	MM_REAL_DENSE  = 1	// dense matrix
} MMRealFormat;

enum {
	MM_GENERAL   = 1 << 0,
	MM_SYMMETRIC = 1 << 1
};

enum {
	MM_UPPER = 1 << 2,
	MM_LOWER = 1 << 3
};

/* symmetric */
typedef enum {
	MM_REAL_GENERAL = MM_GENERAL,	// asymmetric
	MM_REAL_SYMMETRIC_UPPER = MM_SYMMETRIC | MM_UPPER,	// symmetric upper triangular
	MM_REAL_SYMMETRIC_LOWER = MM_SYMMETRIC | MM_LOWER	// symmetric lower triangular
} MMRealSymm;

#define mm_real_is_sparse(a)		mm_is_sparse((a)->typecode)
#define mm_real_is_dense(a)			mm_is_dense((a)->typecode)
#define mm_real_is_symmetric(a)		(mm_is_symmetric((a)->typecode) && ((a)->symm & MM_SYMMETRIC))
#define mm_real_is_upper(a) 		((a)->symm & MM_UPPER)
#define mm_real_is_lower(a) 		((a)->symm & MM_LOWER)

// MatrixMarket format matrix
typedef struct s_mm_real	mm_real;
typedef struct s_mm_real	mm_dense;
typedef struct s_mm_real	mm_sparse;

/*** implementation of dense / sparse matrix
 * In the case of dense matrix, nnz = m * n and
 * i = NULL, p = NULL. ***/
struct s_mm_real {
	MM_typecode	typecode;	// type of matrix. see mmio.h

	MMRealSymm	symm;

	MM_INT			m;			// num of rows of matrix
	MM_INT			n;			// num of columns
	MM_INT			nnz;		// num of nonzero matrix elements

	MM_INT			*i;			// row index of each nonzero elements: size = nnz
	MM_INT			*p;			// p[0] = 0, p[j+1] = num of nonzeros in X(:,1:j): size = n + 1
	MM_DBL		*data;		// nonzero matrix elements: size = nnz

	bool		owner;
};

mm_real		*mm_real_new (MMRealFormat format, MMRealSymm symm, MM_INT m, MM_INT n, MM_INT nnz);
mm_real		*mm_real_view_array (MMRealFormat format, MMRealSymm symm, MM_INT m, MM_INT n, MM_INT nnz, MM_DBL *data);
void		mm_real_free (mm_real *mm);
bool		mm_real_realloc (mm_real *mm, MM_INT nnz);

void		mm_real_sort (mm_real *x);

void		mm_real_memcpy (mm_real *dest, mm_real *src);
mm_real		*mm_real_copy (mm_real *mm);
mm_dense	*mm_real_copy_sparse_to_dense (mm_sparse *s);
mm_sparse	*mm_real_copy_dense_to_sparse (mm_dense *x, MM_DBL threshold);
void		mm_real_set_all (mm_real *mm, MM_DBL val);

bool		mm_real_sparse_to_dense (mm_sparse *s);
bool		mm_real_dense_to_sparse (mm_dense *d, MM_DBL threshold);
bool		mm_real_symmetric_to_general (mm_real *x);
bool		mm_real_general_to_symmetric (const char uplo, mm_real *x);

MM_INT		mm_real_xj_iamax (mm_real *x, MM_INT j);
MM_INT		mm_real_iamax (mm_real *x);
mm_dense	*mm_real_xj_col (mm_real *x, MM_INT j);
mm_dense	*mm_real_xi_row (mm_real *x, MM_INT i);
void		mm_real_transpose (mm_real *x);

mm_real		*mm_real_eye (MMRealFormat type, MM_INT n);

mm_real		*mm_real_vertcat (mm_real *x1, mm_real *x2);
mm_real		*mm_real_horzcat (mm_real *x1, mm_real *x2);

void		mm_real_xj_add (mm_real *x, MM_INT j, MM_DBL alpha);
void		mm_real_xj_scale (mm_real *x, MM_INT j, MM_DBL alpha);

void		mm_real_scale (mm_real *x, MM_DBL alpha);

MM_DBL		mm_real_xj_asum (mm_real *x, MM_INT j);
MM_DBL		mm_real_xj_sum (mm_real *x, MM_INT j);
MM_DBL		mm_real_xj_ssq (mm_real *x, MM_INT j);
MM_DBL		mm_real_xj_nrm2 (mm_real *x, MM_INT j);

void		mm_real_x_dot_yk (
				bool trans, MM_DBL alpha, mm_real *x,
				mm_dense *y, MM_INT k, MM_DBL beta, mm_dense *z
			);

void		mm_real_x_dot_y (
				bool transx, bool transy, MM_DBL alpha, mm_real *x,
				mm_dense *y, MM_DBL beta, mm_dense *z
			);

MM_DBL		mm_real_xj_trans_dot_yk (mm_real *x, MM_INT j, mm_dense *y, MM_INT k);
mm_dense	*mm_real_xj_trans_dot_y (mm_real *x, MM_INT j, mm_dense *y);

void		mm_real_axjpy (MM_DBL alpha, mm_real *x, MM_INT j, mm_dense *y);
void		mm_real_axjpy_atomic (MM_DBL alpha, mm_real *x, MM_INT j, mm_dense *y);

mm_real		*mm_real_fread (FILE *fp);
void		mm_real_fwrite (FILE *stream, mm_real *x, const char *format);

mm_real		*mm_real_fread_binary (FILE *fp);
void		mm_real_fwrite_binary (FILE *fp, mm_real *x);

#ifdef __cplusplus
}
#endif

#endif /* MMREAL_H_ */
