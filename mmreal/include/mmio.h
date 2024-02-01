/*
*   Matrix Market I/O library for ANSI C
*
*   See http://math.nist.gov/MatrixMarket for details.
*
*
*/
#ifndef MMIO_H
#define MMIO_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

#ifdef MMINT
	#define	MM_INT	int
	#define	MM_SIZ	int
	#define MM_FMT	"%d"
#else
	#define	MM_INT	long
	#define	MM_SIZ	size_t
	#define MM_FMT	"%ld"
#endif

#ifdef MMFLOAT
	#define	MM_DBL	float
#else
	#define	MM_DBL	double
#endif

#define MM_MAX_LINE_LENGTH	1025
#define MatrixMarketBanner	"%%MatrixMarket"
#define MM_MAX_TOKEN_LENGTH	64

typedef char MM_typecode[4];

/* MM_typecode matcode -> const MM_typecode matcode, 2014.11.10, by M.Utsugi */
char	*mm_typecode_to_str(const MM_typecode matcode);

int		mm_read_banner(FILE *f, MM_typecode *matcode);
int		mm_read_mtx_crd_size(FILE *f, MM_INT *M, MM_INT *N, MM_INT *nz);
int		mm_read_mtx_array_size(FILE *f, MM_INT *M, MM_INT *N);

/* MM_typecode matcode -> const MM_typecode matcode, 2014.11.10, by M.Utsugi */
int		mm_write_banner(FILE *f, const MM_typecode matcode);
/* MM_INT M, N, nz -> const MM_INT M, N, nz, 2014.11.10, by M.Utsugi */
int		mm_write_mtx_crd_size(FILE *f, const MM_INT M, const MM_INT N, const MM_INT nz);
/* MM_INT M, N -> const MM_INT M, N, 2014.11.10, by M.Utsugi */
int		mm_write_mtx_array_size(FILE *f, const MM_INT M, const MM_INT N);


/********************* MM_typecode query fucntions ***************************/

#define mm_is_matrix(typecode)	((typecode)[0]=='M')

#define mm_is_sparse(typecode)	((typecode)[1]=='C')
#define mm_is_coordinate(typecode)((typecode)[1]=='C')
#define mm_is_dense(typecode)	((typecode)[1]=='A')
#define mm_is_array(typecode)	((typecode)[1]=='A')

#define mm_is_complex(typecode)	((typecode)[2]=='C')
#define mm_is_real(typecode)		((typecode)[2]=='R')
#define mm_is_pattern(typecode)	((typecode)[2]=='P')
#define mm_is_integer(typecode) ((typecode)[2]=='I')

#define mm_is_symmetric(typecode)((typecode)[3]=='S')
#define mm_is_general(typecode)	((typecode)[3]=='G')
#define mm_is_skew(typecode)	((typecode)[3]=='K')
#define mm_is_hermitian(typecode)((typecode)[3]=='H')

/* MM_typecode matcode -> const MM_typecode matcode, 2014.11.10, by M.Utsugi */
int		mm_is_valid(const MM_typecode matcode);		/* too complex for a macro */


/********************* MM_typecode modify fucntions ***************************/

#define mm_set_matrix(typecode)	((*typecode)[0]='M')
#define mm_set_coordinate(typecode)	((*typecode)[1]='C')
#define mm_set_array(typecode)	((*typecode)[1]='A')
#define mm_set_dense(typecode)	mm_set_array(typecode)
#define mm_set_sparse(typecode)	mm_set_coordinate(typecode)

#define mm_set_complex(typecode)((*typecode)[2]='C')
#define mm_set_real(typecode)	((*typecode)[2]='R')
#define mm_set_pattern(typecode)((*typecode)[2]='P')
#define mm_set_integer(typecode)((*typecode)[2]='I')


#define mm_set_symmetric(typecode)((*typecode)[3]='S')
#define mm_set_general(typecode)((*typecode)[3]='G')
#define mm_set_skew(typecode)	((*typecode)[3]='K')
#define mm_set_hermitian(typecode)((*typecode)[3]='H')

#define mm_clear_typecode(typecode) ((*typecode)[0]=(*typecode)[1]= \
									(*typecode)[2]=' ',(*typecode)[3]='G')

#define mm_initialize_typecode(typecode) mm_clear_typecode(typecode)


/********************* Matrix Market error codes ***************************/


#define MM_COULD_NOT_READ_FILE	11
#define MM_PREMATURE_EOF		12
#define MM_NOT_MTX				13
#define MM_NO_HEADER			14
#define MM_UNSUPPORTED_TYPE		15
#define MM_LINE_TOO_LONG		16
#define MM_COULD_NOT_WRITE_FILE	17


/******************** Matrix Market internal definitions ********************

   MM_matrix_typecode: 4-character sequence

				    ojbect 		sparse/   	data        storage
						  		dense     	type        scheme

   string position:	 [0]        [1]			[2]         [3]

   Matrix typecode:  M(atrix)  C(oord)		R(eal)   	G(eneral)
						        A(array)	C(omplex)   H(ermitian)
											P(attern)   S(ymmetric)
								    		I(nteger)	K(kew)

 ***********************************************************************/

#define MM_MTX_STR			"matrix"
#define MM_ARRAY_STR		"array"
#define MM_DENSE_STR		"array"
#define MM_COORDINATE_STR	"coordinate"
#define MM_SPARSE_STR		"coordinate"
#define MM_COMPLEX_STR		"complex"
#define MM_REAL_STR			"real"
#define MM_INT_STR			"integer"
#define MM_GENERAL_STR		"general"
#define MM_SYMM_STR			"symmetric"
#define MM_HERM_STR			"hermitian"
#define MM_SKEW_STR			"skew-symmetric"
#define MM_PATTERN_STR		"pattern"


/*  high level routines */

/* char fname[] -> const char fname[],
 * MM_INT M, N, nz, I[], J[] -> const MM_INT M, N, nz, I[], J[],
 * MM_DBL val[] -> const MM_DBL val[],
 * MM_typecode matcode -> const MM_typecode matcode, 2014.11.10, by M.Utsugi */
int		mm_write_mtx_crd(const char fname[], const MM_INT M, const MM_INT N, const MM_INT nz,
			const MM_INT I[], const MM_INT J[],
			const MM_DBL val[], const MM_typecode matcode);
int		mm_read_mtx_crd_data(FILE *f, MM_INT M, MM_INT N, MM_INT nz, MM_INT I[], MM_INT J[],
			MM_DBL val[], MM_typecode matcode);
int		mm_read_mtx_crd_entry(FILE *f, MM_INT *I, MM_INT *J, MM_DBL *real, MM_DBL *img,
			MM_typecode matcode);

int		mm_read_unsymmetric_sparse(const char *fname, MM_INT *M_, MM_INT *N_, MM_INT *nz_,
			MM_DBL **val_, MM_INT **I_, MM_INT **J_);

#ifdef __cplusplus
}
#endif

#endif	// MMIO_H
