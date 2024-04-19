#include "bl_config.h"
#include "bl_dgemm_kernel.h"

#define a(i, j, ld) a[(i) * (ld) + (j)]
#define b(i, j, ld) b[(i) * (ld) + (j)]
#define c(i, j, ld) c[(i) * (ld) + (j)]

//
// C-based micorkernel
//
//#define NOPACK
#ifdef NOPACK
void bl_dgemm_ukr(int k, int m, int n, double *a, double *b, double *c,
                  unsigned long long ldc, aux_t *data) {
  int l, j, i;

  for (l = 0; l < k; ++l) {
    for (j = 0; j < n; ++j) {
      for (i = 0; i < m; ++i) {
        // ldc is used here because a[] and b[] are not packed by the
        // starter code
        // cse260 - you can modify the leading indice to DGEMM_NR and DGEMM_MR
        // as appropriate
        //
        c(i, j, ldc) += a(i, l, ldc) * b(l, j, ldc);
      }
    }
  }
}
#else
void bl_dgemm_ukr(int k, int m, int n, double *a, double *b, double *c,
                  unsigned long long ldc, aux_t *data) {
  int l, j, i;
  int q;
  ///double *cp;
  ///double *ap;
  ///double *bp;
  size_t ds = (int)sizeof(double)*8;
  int shift = ds - 1;
  for (l = 0; l < k; ++l) {
    for (i = 0; i < m; ++i) {
      //cp = &c[i*ldc];
      for (j = 0; j < n; j=j+1) {
        //c(i, j, ldc) += a(i, l, k) * b(l, j, n);
        c(i, j, ldc) += a(l, i, m) * b(l, j, n);
	//uint64_t m1,m2,m3,m4,m5,m6,m7,m8;
	//m1 = n-j-1;
	//m2 = n-j-2;
	//m3 = n-j-3;
	//m4 = n-j-4;
	//m5 = n-j-5;
	//m6 = n-j-6;
	//m7 = n-j-7;
	//m8 = n-j-8;
	////*(cp+0) += a(i, l, k) * b(l, j+0, n) * ((~(*(uint64_t*)&m1 >> ds-1))&1);
	////*(cp+1) += a(i, l, k) * b(l, j+1, n) * ((~(*(uint64_t*)&m2 >> ds-1))&1);
	////*(cp+2) += a(i, l, k) * b(l, j+2, n) * ((~(*(uint64_t*)&m3 >> ds-1))&1);
	////*(cp+3) += a(i, l, k) * b(l, j+3, n) * ((~(*(uint64_t*)&m4 >> ds-1))&1);
	//bp = &b[l*n+j];
	//ap= &a[i*k+l];
	//*(cp+0) += *ap * *(bp+0) * ((~(m1 >> shift ))&1);
	//*(cp+1) += *ap * *(bp+1) * ((~(m2 >> shift ))&1);
	//*(cp+2) += *ap * *(bp+2) * ((~(m3 >> shift ))&1);
	//*(cp+3) += *ap * *(bp+3) * ((~(m4 >> shift ))&1);
	//*(cp+4) += *ap * *(bp+4) * ((~(m5 >> shift ))&1);
	//*(cp+5) += *ap * *(bp+5) * ((~(m6 >> shift ))&1);
	//*(cp+6) += *ap * *(bp+6) * ((~(m7 >> shift ))&1);
	//*(cp+7) += *ap * *(bp+7) * ((~(m8 >> shift ))&1);
	//cp += 8;
	//*(cp+2) += a(i, l, k) * b(l, 2, n);
	//*(cp+3) += a(i, l, k) * b(l, 3, n);
	//*(cp+4) += a(i, l, k) * b(l, 4, n);
      }
    }
  }
}
#endif

// cse260
// you can put your optimized kernels here
// - put the function prototypes in bl_dgemm_kernel.h
// - define BL_MICRO_KERNEL appropriately in bl_config.h
//
