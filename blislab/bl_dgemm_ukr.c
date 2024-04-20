#include "bl_config.h"
#include "bl_dgemm_kernel.h"

#define a(i, j, ld) a[(i) * (ld) + (j)]
#define b(i, j, ld) b[(i) * (ld) + (j)]
#define c(i, j, ld) c[(i) * (ld) + (j)]

//
// C-based micorkernel
//
//#define NOPACK
#define ARM_SVE_KERNEL
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
	#ifndef ARM_SVE_KERNEL
		void bl_dgemm_ukr(int k, int m, int n, double *a, double *b, double *c,
		                  unsigned long long ldc, aux_t *data) {
		  int l, j, i;
		  int q;
		  size_t ds = (int)sizeof(double)*8;
		  int shift = ds - 1;
		  for (l = 0; l < k; ++l) {
		    for (i = 0; i < m; ++i) {
		      for (j = 0; j < n; j=j+1) {
		        c(i, j, ldc) += a(l, i, m) * b(l, j, n);
		      }
		    }
		  }
		}
	#else
		void bl_dgemm_ukr(int k, int m, int n, double *a, double *b, double *c,
		                  unsigned long long ldc, aux_t *data) {
		  //for (l = 0; l < k; ++l) {
		  //  for (i = 0; i < m; ++i) {
		  //    for (j = 0; j < n; j=j+1) {
		  //      c(i, j, ldc) += a(l, i, m) * b(l, j, n);
		  //    }
		  //  }
		  //}
		  register svfloat64_t ax;
		  register svfloat64_t bx;
		  //register svfloat64_t cx;
		  register svfloat64_t clx;
		  svbool_t npred = svwhilelt_b64_u64(0,n);
		  int i,j,l;
		  for (l = 0; l < k; ++l) {
		    bx = svld1_f64(svptrue_b64(), b + l*n );
		    for (i = 0; i < m; ++i) {
		      clx = svld1_f64(npred,c + i*ldc);	
		      register float64_t aval = *(a + l*m + i);
		      ax = svdup_f64(aval);
		      clx = svmla_f64_m(npred,clx, bx, ax);
		    svst1_f64(npred, c + i*ldc, clx);
		    }
		  }
		}
	#endif
#endif

// cse260
// you can put your optimized kernels here
// - put the function prototypes in bl_dgemm_kernel.h
// - define BL_MICRO_KERNEL appropriately in bl_config.h
//
