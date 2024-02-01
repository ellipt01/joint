#ifndef BLAS_INTERFACE_HEADER
#define BLAS_INTERFACE_HEADER

#ifdef __cplusplus 
extern "C"{
#endif

//Structs
typedef struct complex_Tag{
  float r;
  float i;
} complex;

typedef struct doublecomplex_Tag{
  double r;
  double i;
} doublecomplex;


//MM_INT xrebla_(char *srname, MM_INT *info);

//Level1

//AXPY
void saxpy_(MM_INT *n, float         *alpha, float         *x, MM_INT *incx, float         *y, MM_INT *incy);
void daxpy_(MM_INT *n, double        *alpha, double        *x, MM_INT *incx, double        *y, MM_INT *incy);
void caxpy_(MM_INT *n, complex       *alpha, complex       *x, MM_INT *incx, complex       *y, MM_INT *incy);
void zaxpy_(MM_INT *n, doublecomplex *alpha, doublecomplex *x, MM_INT *incx, doublecomplex *y, MM_INT *incy);

//SUM
float   sasum_(MM_INT *n, float         *x, MM_INT *incx);
float  scasum_(MM_INT *n, complex       *x, MM_INT *incx);
double  dasum_(MM_INT *n, double        *x, MM_INT *incx);
double dzasum_(MM_INT *n, doublecomplex *x, MM_INT *incx);

//COPY
void scopy_(MM_INT *n, float  *x, MM_INT *incx, float  *y, MM_INT *incy);
void dcopy_(MM_INT *n, double *x, MM_INT *incx, double *y, MM_INT *incy);
void ccopy_(MM_INT *n, float  *x, MM_INT *incx, float  *y, MM_INT *incy);
void zcopy_(MM_INT *n, double *x, MM_INT *incx, double *y, MM_INT *incy);

//DOT
float  sdot_(MM_INT *n, float  *x, MM_INT *incx, float  *y, MM_INT *incy);
double ddot_(MM_INT *n, double *x, MM_INT *incx, double *y, MM_INT *incy);

//DOTC
complex       cdotc_(MM_INT *n, complex       *x, MM_INT *incx, complex       *y, MM_INT *incy);
doublecomplex zdotc_(MM_INT *n, doublecomplex *x, MM_INT *incx, doublecomplex *y, MM_INT *incy);

//DOTU
complex       cdotu_(MM_INT *n, complex       *x, MM_INT *incx, complex       *y, MM_INT *incy);
doublecomplex zdotu_(MM_INT *n, doublecomplex *x, MM_INT *incx, doublecomplex *y, MM_INT *incy);

//NRM2
float   snrm2_(MM_INT *n, float         *x, MM_INT *incx);
double  dnrm2_(MM_INT *n, double        *x, MM_INT *incx);
float  scnrm2_(MM_INT *n, complex       *x, MM_INT *incx);
double dznrm2_(MM_INT *n, doublecomplex *x, MM_INT *incx);

//ROT
void  srot_(MM_INT *n, float         *x, MM_INT *incx, float         *y, MM_INT *incy, float  *c, float  *s);
void  drot_(MM_INT *n, double        *x, MM_INT *incx, double        *y, MM_INT *incy, double *c, double *s);
void csrot_(MM_INT *n, complex       *x, MM_INT *incx, complex       *y, MM_INT *incy, float  *c, float  *s);
void zdrot_(MM_INT *n, doublecomplex *x, MM_INT *incx, doublecomplex *y, MM_INT *incy, double *c, double *s);

//ROTG
void srotg_(float         *a,float         *b, float  *c, float  *s);
void drotg_(double        *a,double        *b, double *c, double *s);
void crotg_(complex       *a,complex       *b, float  *c, float  *s);
void zrotg_(doublecomplex *a,doublecomplex *b, double *c, double *s);

//Stub
//ROTMG
//ROTM


//SCAL
void  sscal_(MM_INT *n,  float         *a, float          *x, MM_INT *incx);
void  dscal_(MM_INT *n,  double        *a, double         *x, MM_INT *incx);
void  cscal_(MM_INT *n,  complex       *a, complex        *x, MM_INT *incx);
void  zscal_(MM_INT *n,  doublecomplex *a, doublecomplex  *x, MM_INT *incx);
void csscal_(MM_INT *n,  float         *a, complex        *x, MM_INT *incx);
void zdscal_(MM_INT *n,  double        *a, doublecomplex  *x, MM_INT *incx);

//SWAP
void sswap_(MM_INT *n, float         *x, MM_INT *incx, float         *y, MM_INT *incy);
void dswap_(MM_INT *n, double        *x, MM_INT *incx, double        *y, MM_INT *incy);
void cswap_(MM_INT *n, complex       *x, MM_INT *incx, complex       *y, MM_INT *incy);
void zswap_(MM_INT *n, doublecomplex *x, MM_INT *incx, doublecomplex *y, MM_INT *incy);

//IAMAX
MM_INT isamax_(MM_INT *n, float         *x, MM_INT *incx);
MM_INT idamax_(MM_INT *n, double        *x, MM_INT *incx);
MM_INT icamax_(MM_INT *n, complex       *x, MM_INT *incx);
MM_INT izamax_(MM_INT *n, doublecomplex *x, MM_INT *incx);

//IAMIN
MM_INT isamin_(MM_INT *n, float         *x, MM_INT *incx);
MM_INT idamin_(MM_INT *n, double        *x, MM_INT *incx);
MM_INT icamin_(MM_INT *n, complex       *x, MM_INT *incx);
MM_INT izamin_(MM_INT *n, doublecomplex *x, MM_INT *incx);

//IMAX
MM_INT ismax_(MM_INT *n, float  *x, MM_INT *incx);
MM_INT idmax_(MM_INT *n, double *x, MM_INT *incx);

//IMIN
MM_INT ismin_(MM_INT *n, float  *x, MM_INT *incx);
MM_INT idmin_(MM_INT *n, double *x, MM_INT *incx);

//Level2

//GBMV
void sgbmv_(char *trans, MM_INT *m, MM_INT *n, MM_INT *kl, MM_INT *ku, 
            float         *alpha, float         *A, MM_INT *ldA, float         *x, MM_INT *incx,
            float         *beta , float         *y, MM_INT *incy);
void dgbmv_(char *trans, MM_INT *m, MM_INT *n, MM_INT *kl, MM_INT *ku, 
            double        *alpha, double        *A, MM_INT *ldA, double        *x, MM_INT *incx, 
            double        *beta , double        *y, MM_INT *incy);
void cgbmv_(char *trans, MM_INT *m, MM_INT *n, MM_INT *kl, MM_INT *ku, 
            complex       *alpha, complex       *A, MM_INT *ldA, complex       *x, MM_INT *incx, 
            complex       *beta , complex       *y, MM_INT *incy);
void zgbmv_(char *trans, MM_INT *m, MM_INT *n, MM_INT *kl, MM_INT *ku, 
            doublecomplex *alpha, doublecomplex *A, MM_INT *ldA, doublecomplex *x, MM_INT *incx, 
            doublecomplex *beta , doublecomplex *y, MM_INT *incy);

//GEMV
void sgemv_(char *trans, MM_INT *m, MM_INT *n, 
            float         *alpha, float         *A, MM_INT *ldA, float         *x, MM_INT *incx,
            float         *beta , float         *y, MM_INT *incy);
void dgemv_(char *trans, MM_INT *m, MM_INT *n, 
            double        *alpha, double        *A, MM_INT *ldA, double        *x, MM_INT *incx, 
            double        *beta , double        *y, MM_INT *incy);
void cgemv_(char *trans, MM_INT *m, MM_INT *n, 
            complex       *alpha, complex       *A, MM_INT *ldA, complex       *x, MM_INT *incx, 
            complex       *beta , complex       *y, MM_INT *incy);
void zgemv_(char *trans, MM_INT *m, MM_INT *n, 
            doublecomplex *alpha, doublecomplex *A, MM_INT *ldA, doublecomplex *x, MM_INT *incx, 
            doublecomplex *beta , doublecomplex *y, MM_INT *incy);

//GER
void sger_(MM_INT *m, MM_INT *n, float  *alpha,  float *x, MM_INT *incx,  float *y, MM_INT *incy,  float *A, MM_INT *ldA);
void dger_(MM_INT *m, MM_INT *n, double *alpha, double *x, MM_INT *incx, double *y, MM_INT *incy, double *A, MM_INT *ldA);

//GERC
void cgerc_(MM_INT *m, MM_INT *n, complex       *alpha, complex       *x, MM_INT *incx,
            complex       *y, MM_INT *incy, complex       *A, MM_INT *ldA);
void zgerc_(MM_INT *m, MM_INT *n, doublecomplex *alpha, doublecomplex *x, MM_INT *incx, 
            doublecomplex *y, MM_INT *incy, doublecomplex *A, MM_INT *ldA);

//GREU
void cgeru_(MM_INT *m, MM_INT *n, complex       *alpha, complex       *x, MM_INT *incx,
            complex       *y, MM_INT *incy, complex       *A, MM_INT *ldA);
void zgeru_(MM_INT *m, MM_INT *n, doublecomplex *alpha, doublecomplex *x, MM_INT *incx, 
            doublecomplex *y, MM_INT *incy, doublecomplex *A, MM_INT *ldA);

//HBMV
void chbmv_(char *uplo, MM_INT *n, MM_INT *k, complex       *alpha, complex       *A, MM_INT *ldA,
            complex       *x, MM_INT *incx, complex       *beta, complex       *y, MM_INT *incy);
void zhbmv_(char *uplo, MM_INT *n, MM_INT *k, doublecomplex *alpha, doublecomplex *A, MM_INT *ldA,
            doublecomplex *x, MM_INT *incx, doublecomplex *beta, doublecomplex *y, MM_INT *incy);

//HEMV
void chemv_(char *uplo, MM_INT *n, complex       *alpha, complex       *A, MM_INT *ldA,
            complex       *x, MM_INT *incx, complex       *beta, complex       *y, MM_INT *incy);
void zhemv_(char *uplo, MM_INT *n, doublecomplex *alpha, doublecomplex *A, MM_INT *ldA,
            doublecomplex *x, MM_INT *incx, doublecomplex *beta, doublecomplex *y, MM_INT *incy);

//HER
void cher_(char *uplo, MM_INT *n, float  *alpha, complex       *x, MM_INT *incx, complex       *A, MM_INT *ldA);
void zher_(char *uplo, MM_INT *n, double *alpha, doublecomplex *x, MM_INT *incx, doublecomplex *A, MM_INT *ldA);

//Stub
//HER2

//HPMV
void chpmv_(char *uplo, MM_INT *n, complex       *alpha, complex       *A,
            complex       *x, MM_INT *incx, complex       *beta, complex       *y, MM_INT *incy);
void zhpmv_(char *uplo, MM_INT *n, doublecomplex *alpha, doublecomplex *A,
            doublecomplex *x, MM_INT *incx, doublecomplex *beta, doublecomplex *y, MM_INT *incy);

//HPR
void chpr_ (char *uplo, MM_INT *n, float   *alpha, complex       *x, MM_INT *incx, complex       *A);
void zhpr_ (char *uplo, MM_INT *n, double  *alpha, doublecomplex *x, MM_INT *incx, doublecomplex *A);

//Stub
//HPR2

//SBMV
void ssbmv_(char *uplo, MM_INT *n, MM_INT *k, float  *alpha, float  *A, MM_INT *ldA,
            float  *x, MM_INT *incx, float  *beta, float  *y, MM_INT *incy);
void dsbmv_(char *uplo, MM_INT *n, MM_INT *k, double *alpha, double *A, MM_INT *ldA,
            double *x, MM_INT *incx, double *beta, double *y, MM_INT *incy);

//SPMV
void sspmv_(char *uplo, MM_INT *n, float  *alpha, float  *A, float  *x, MM_INT *incx, float  *beta, float  *y, MM_INT *incy);
void dspmv_(char *uplo, MM_INT *n, double *alpha, double *A, double *x, MM_INT *incx, double *beta, double *y, MM_INT *incy);

//SPR
void sspr_(char *uplo, MM_INT *n, float  *alpha, float  *x, MM_INT *incx, float  *A);
void dspr_(char *uplo, MM_INT *n, double *alpha, double *x, MM_INT *incx, double *A);

//Stub
//SPR2

//SYMV
void ssymv_(char *uplo, MM_INT *n, float  *alpha, float  *A, MM_INT *ldA,
            float  *x, MM_INT *incx, float  *beta, float  *y, MM_INT *incy);
void dsymv_(char *uplo, MM_INT *n, double *alpha, double *A, MM_INT *ldA,
            double *x, MM_INT *incx, double *beta, double *y, MM_INT *incy);

//SYR
void ssyr_(char *uplo, MM_INT *n, float  *alpha, float  *x, MM_INT *incx, float  *A, MM_INT *ldA);
void dsyr_(char *uplo, MM_INT *n, double *alpha, double *x, MM_INT *incx, double *A, MM_INT *ldA);

//Stub
//SYR2

//TBMV
void stbmv_(char *uplo, char *trans, char *diag, MM_INT *n, MM_INT *k, float         *A, MM_INT *ldA, float         *x, MM_INT *incx);
void dtbmv_(char *uplo, char *trans, char *diag, MM_INT *n, MM_INT *k, double        *A, MM_INT *ldA, double        *x, MM_INT *incx);
void ctbmv_(char *uplo, char *trans, char *diag, MM_INT *n, MM_INT *k, complex       *A, MM_INT *ldA, complex       *x, MM_INT *incx);
void ztbmv_(char *uplo, char *trans, char *diag, MM_INT *n, MM_INT *k, doublecomplex *A, MM_INT *ldA, doublecomplex *x, MM_INT *incx);

//TBSV
void stbsv_(char *uplo, char *trans, char *diag, MM_INT *n, MM_INT *k, float         *A, MM_INT *ldA, float         *x, MM_INT *incx);
void dtbsv_(char *uplo, char *trans, char *diag, MM_INT *n, MM_INT *k, double        *A, MM_INT *ldA, double        *x, MM_INT *incx);
void ctbsv_(char *uplo, char *trans, char *diag, MM_INT *n, MM_INT *k, complex       *A, MM_INT *ldA, complex       *x, MM_INT *incx);
void ztbsv_(char *uplo, char *trans, char *diag, MM_INT *n, MM_INT *k, doublecomplex *A, MM_INT *ldA, doublecomplex *x, MM_INT *incx);

//TPMV
void stpmv_(char *uplo, char *trans, char *diag, MM_INT *n, float         *A, float         *x, MM_INT *incx);
void dtpmv_(char *uplo, char *trans, char *diag, MM_INT *n, double        *A, double        *x, MM_INT *incx);
void ctpmv_(char *uplo, char *trans, char *diag, MM_INT *n, complex       *A, complex       *x, MM_INT *incx);
void ztpmv_(char *uplo, char *trans, char *diag, MM_INT *n, doublecomplex *A, doublecomplex *x, MM_INT *incx);

//TPSV
void stpsv_(char *uplo, char *trans, char *diag, MM_INT *n, float         *A, float         *x, MM_INT *incx);
void dtpsv_(char *uplo, char *trans, char *diag, MM_INT *n, double        *A, double        *x, MM_INT *incx);
void ctpsv_(char *uplo, char *trans, char *diag, MM_INT *n, complex       *A, complex       *x, MM_INT *incx);
void ztpsv_(char *uplo, char *trans, char *diag, MM_INT *n, doublecomplex *A, doublecomplex *x, MM_INT *incx);

//TRSV
void strsv_(char *uplo, char *trans, char *diag, MM_INT *n, float         *A, MM_INT *ldA, float         *x, MM_INT *incx);
void dtrsv_(char *uplo, char *trans, char *diag, MM_INT *n, double        *A, MM_INT *ldA, double        *x, MM_INT *incx);
void ctrsv_(char *uplo, char *trans, char *diag, MM_INT *n, complex       *A, MM_INT *ldA, complex       *x, MM_INT *incx);
void ztrsv_(char *uplo, char *trans, char *diag, MM_INT *n, doublecomplex *A, MM_INT *ldA, doublecomplex *x, MM_INT *incx);

//TRMV
void strmv_(char *uplo, char *trans, char *diag, MM_INT *n, float         *A, MM_INT *ldA, float         *x, MM_INT *incx);
void dtrmv_(char *uplo, char *trans, char *diag, MM_INT *n, double        *A, MM_INT *ldA, double        *x, MM_INT *incx);
void ctrmv_(char *uplo, char *trans, char *diag, MM_INT *n, complex       *A, MM_INT *ldA, complex       *x, MM_INT *incx);
void ztrmv_(char *uplo, char *trans, char *diag, MM_INT *n, doublecomplex *A, MM_INT *ldA, doublecomplex *x, MM_INT *incx);

//Level3

//GEMM
void sgemm_(char *transa, char *transb, MM_INT *m, MM_INT *n, MM_INT *k,
            float         *alpha, float         *A, MM_INT *ldA, float         *B, MM_INT *ldB, 
            float         *beta , float         *C, MM_INT *ldC);
void dgemm_(char *transa, char *transb, MM_INT *m, MM_INT *n, MM_INT *k, 
            double        *alpha, double        *A, MM_INT *ldA, double        *B, MM_INT *ldB, 
            double        *beta , double        *C, MM_INT *ldC);
void cgemm_(char *transa, char *transb, MM_INT *m, MM_INT *n, MM_INT *k, 
            complex       *alpha, complex       *A, MM_INT *ldA, complex       *B, MM_INT *ldB, 
            complex       *beta , complex       *C, MM_INT *ldC);
void zgemm_(char *transa, char *transb, MM_INT *m, MM_INT *n, MM_INT *k,
            doublecomplex *alpha, doublecomplex *A, MM_INT *ldA, doublecomplex *B, MM_INT *ldB, 
            doublecomplex *beta , doublecomplex *C, MM_INT *ldC);

//HEMM
void chemm_(char *side, char *uplo, MM_INT *m, MM_INT *n, complex       *alpha, complex       *A, MM_INT *ldA,
            complex       *B, MM_INT *ldB, complex       *beta, complex       *C, MM_INT *ldC);
void zhemm_(char *side, char *uplo, MM_INT *m, MM_INT *n, doublecomplex *alpha, doublecomplex *A, MM_INT *ldA,
            doublecomplex *B, MM_INT *ldB, doublecomplex *beta, doublecomplex *C, MM_INT *ldC);

//HERK
void cherk_(char *uplo, char *trans, MM_INT *n, MM_INT *k, float  *alpha, complex       *A, MM_INT *ldA,
            float  *beta , complex       *C, MM_INT *ldC);
void zherk_(char *uplo, char *trans, MM_INT *n, MM_INT *k, double *aplha, doublecomplex *A, MM_INT *ldA,
            double *beta , doublecomplex *C, MM_INT *ldC);

//HERK2
void cher2k_(char *uplo, char *trans, MM_INT *n, MM_INT *k, complex       *alpha, complex       *A, MM_INT *ldA,
             complex       *B, MM_INT *ldB,float  *beta, complex       *C, MM_INT *ldC);
void zher2k_(char *uplo, char *trans, MM_INT *n, MM_INT *k, doublecomplex *alpha, doublecomplex *A, MM_INT *ldA,
             doublecomplex *B, MM_INT *ldB,double *beta, doublecomplex *C, MM_INT *ldC);

//SYMM
void ssymm_(char *side, char *uplo, MM_INT *m, MM_INT *n, 
            float         *alpha, float         *A, MM_INT *ldA, float         *B, MM_INT *ldB, 
            float         *beta , float         *C, MM_INT *ldC);
void dsymm_(char *side, char *uplo, MM_INT *m, MM_INT *n, 
            double        *alpha, double        *A, MM_INT *ldA, double        *B, MM_INT *ldB, 
            double        *beta , double        *C, MM_INT *ldC);
void csymm_(char *side, char *uplo, MM_INT *m, MM_INT *n, 
            complex       *alpha, complex       *A, MM_INT *ldA, complex       *B, MM_INT *ldB, 
            complex       *beta , complex       *C, MM_INT *ldC);
void zsymm_(char *side, char *uplo, MM_INT *m, MM_INT *n, 
            doublecomplex *alpha, doublecomplex *A, MM_INT *ldA, doublecomplex *B, MM_INT *ldB, 
            doublecomplex *beta , doublecomplex *C, MM_INT *ldC);

//SYRK
void ssyrk_(char *uplo, char *trans, MM_INT *n, MM_INT *k, float  *alpha, float         *A, MM_INT *ldA,
            float  *beta , float         *C, MM_INT *ldC);
void dsyrk_(char *uplo, char *trans, MM_INT *n, MM_INT *k, double *alpha, double        *A, MM_INT *ldA,
            double *beta , double        *C, MM_INT *ldC);
void csyrk_(char *uplo, char *trans, MM_INT *n, MM_INT *k, float  *alpha, complex       *A, MM_INT *ldA,
            float  *beta , complex       *C, MM_INT *ldC);
void zsyrk_(char *uplo, char *trans, MM_INT *n, MM_INT *k, double *aplha, doublecomplex *A, MM_INT *ldA,
            double *beta , doublecomplex *C, MM_INT *ldC);

//SYR2K
void ssyr2k_(char *uplo, char *trans, MM_INT *n, MM_INT *k, float  *alpha, float         *A, MM_INT *ldA, float         *B, MM_INT *ldB,
            float  *beta , float         *C, MM_INT *ldC);
void dsyr2k_(char *uplo, char *trans, MM_INT *n, MM_INT *k, double *alpha, double        *A, MM_INT *ldA, double        *B, MM_INT *ldB,
            double *beta , double        *C, MM_INT *ldC);
void csyr2k_(char *uplo, char *trans, MM_INT *n, MM_INT *k, float  *alpha, complex       *A, MM_INT *ldA, complex       *B, MM_INT *ldB,
            float  *beta , complex       *C, MM_INT *ldC);
void zsyr2k_(char *uplo, char *trans, MM_INT *n, MM_INT *k, double *aplha, doublecomplex *A, MM_INT *ldA, doublecomplex *B, MM_INT *ldB,
            double *beta , doublecomplex *C, MM_INT *ldC);

//TRMM
void strmm_(char *side, char *uplo, char *trans, char *diag, MM_INT *m, MM_INT *n, 
            float         *alpha, float         *A, MM_INT *ldA, float         *B, MM_INT *ldB);
void dtrmm_(char *side, char *uplo, char *trans, char *diag, MM_INT *m, MM_INT *n, 
            double        *alpha, double        *A, MM_INT *ldA, double        *B, MM_INT *ldB);
void ctrmm_(char *side, char *uplo, char *trans, char *diag, MM_INT *m, MM_INT *n, 
            complex       *alpha, complex       *A, MM_INT *ldA, complex       *B, MM_INT *ldB);
void ztrmm_(char *side, char *uplo, char *trans, char *diag, MM_INT *m, MM_INT *n, 
            doublecomplex *alpha, doublecomplex *A, MM_INT *ldA, doublecomplex *B, MM_INT *ldB);

//TRSM
void strsm_(char *side, char *uplo, char *trans, char *diag, MM_INT *m, MM_INT *n, 
            float         *alpha, float         *A, MM_INT *ldA, float         *B, MM_INT *ldB);
void dtrsm_(char *side, char *uplo, char *trans, char *diag, MM_INT *m, MM_INT *n, 
            double        *alpha, double        *A, MM_INT *ldA, double        *B, MM_INT *ldB);
void ctrsm_(char *side, char *uplo, char *trans, char *diag, MM_INT *m, MM_INT *n, 
            complex       *alpha, complex       *A, MM_INT *ldA, complex       *B, MM_INT *ldB);
void ztrsm_(char *side, char *uplo, char *trans, char *diag, MM_INT *m, MM_INT *n, 
            doublecomplex *alpha, doublecomplex *A, MM_INT *ldA, doublecomplex *B, MM_INT *ldB);

#ifdef __cplusplus 
}
#endif

#endif

