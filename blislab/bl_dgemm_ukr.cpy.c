#include "bl_config.h"
#include "bl_dgemm_kernel.h"

#define a(i, j, ld) a[(i) * (ld) + (j)]
#define b(i, j, ld) b[(i) * (ld) + (j)]
#define c(i, j, ld) c[(i) * (ld) + (j)]
#define SIZEINT ((sizeof(int))-(1))
#define MASK(i,m) ((~((m - i - 1) >> (SIZEINT))) & (1))
//#define MASK0 MASK(0,DGEMM_MR)
//#define MASK1 MASK(1,DGEMM_MR)
//#define MASK2 MASK(2,DGEMM_MR)
//#define MASK3 MASK(3,DGEMM_MR)
//#define MASK4 MASK(4,DGEMM_MR)
//#define MASK5 MASK(5,DGEMM_MR)
//#define MASK6 MASK(6,DGEMM_MR)
//#define MASK7 MASK(7,DGEMM_MR)

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
	printMat(DGEMM_MR,DGEMM_KC,"A:",a);
	printMat(DGEMM_KC,DGEMM_NR,"B:",b);
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
	printMat(DGEMM_MR,DGEMM_NR,"C:",c);
}
#else
void bl_dgemm_ukr(int k, int m, int n, double *restrict a, double *restrict b, double *restrict c,
		unsigned long long ldc, aux_t *restrict data) {
	svbool_t npred = svwhilelt_b64_u64(0,n);
	int i,j,l;
	register svfloat64_t bx;
	//register svfloat64_t ax0; register svfloat64_t clx0; register float64_t aval0;
	//register svfloat64_t ax1; register svfloat64_t clx1; register float64_t aval1;
	//register svfloat64_t ax2; register svfloat64_t clx2; register float64_t aval2;
	//register svfloat64_t ax3; register svfloat64_t clx3; register float64_t aval3;
	//register svfloat64_t ax4; register svfloat64_t clx4; register float64_t aval4;
	//register svfloat64_t ax5; register svfloat64_t clx5; register float64_t aval5;
	//register svfloat64_t ax6; register svfloat64_t clx6; register float64_t aval6;
	//register svfloat64_t ax7; register svfloat64_t clx7; register float64_t aval7;


	double *restrict ap;
	for (l = 0; l < k; ++l) {
		bx = svld1_f64(svptrue_b64(), b + l*n );
		ap = a + l*m;
		if( m>>3 == 1) // Only work for MR == 8
		{
	register svfloat64_t ax0; register svfloat64_t clx0; register float64_t aval0;
			clx0 = svld1_f64(npred,c + 0*ldc); aval0 = (*(ap + 0)); ax0 = svdup_f64(aval0); clx0 = svmla_f64_m(npred,clx0, bx, ax0); svst1_f64(npred, c + 0*ldc, clx0);
	register svfloat64_t ax1; register svfloat64_t clx1; register float64_t aval1;
			clx1 = svld1_f64(npred,c + 1*ldc); aval1 = (*(ap + 1)); ax1 = svdup_f64(aval1); clx1 = svmla_f64_m(npred,clx1, bx, ax1); svst1_f64(npred, c + 1*ldc, clx1);
	register svfloat64_t ax2; register svfloat64_t clx2; register float64_t aval2;
			clx2 = svld1_f64(npred,c + 2*ldc); aval2 = (*(ap + 2)); ax2 = svdup_f64(aval2); clx2 = svmla_f64_m(npred,clx2, bx, ax2); svst1_f64(npred, c + 2*ldc, clx2);
	register svfloat64_t ax3; register svfloat64_t clx3; register float64_t aval3;
			clx3 = svld1_f64(npred,c + 3*ldc); aval3 = (*(ap + 3)); ax3 = svdup_f64(aval3); clx3 = svmla_f64_m(npred,clx3, bx, ax3); svst1_f64(npred, c + 3*ldc, clx3);
	register svfloat64_t ax4; register svfloat64_t clx4; register float64_t aval4;
			clx4 = svld1_f64(npred,c + 4*ldc); aval4 = (*(ap + 4)); ax4 = svdup_f64(aval4); clx4 = svmla_f64_m(npred,clx4, bx, ax4); svst1_f64(npred, c + 4*ldc, clx4);
	register svfloat64_t ax5; register svfloat64_t clx5; register float64_t aval5;
			clx5 = svld1_f64(npred,c + 5*ldc); aval5 = (*(ap + 5)); ax5 = svdup_f64(aval5); clx5 = svmla_f64_m(npred,clx5, bx, ax5); svst1_f64(npred, c + 5*ldc, clx5);
	register svfloat64_t ax6; register svfloat64_t clx6; register float64_t aval6;
			clx6 = svld1_f64(npred,c + 6*ldc); aval6 = (*(ap + 6)); ax6 = svdup_f64(aval6); clx6 = svmla_f64_m(npred,clx6, bx, ax6); svst1_f64(npred, c + 6*ldc, clx6);
	register svfloat64_t ax7; register svfloat64_t clx7; register float64_t aval7;
			clx7 = svld1_f64(npred,c + 7*ldc); aval7 = (*(ap + 7)); ax7 = svdup_f64(aval7); clx7 = svmla_f64_m(npred,clx7, bx, ax7); svst1_f64(npred, c + 7*ldc, clx7);
		}
		else
		{
			for (i = 0; i < m; ++i) {
	register svfloat64_t ax; register svfloat64_t clx; register float64_t aval;
				clx = svld1_f64(npred,c + i*ldc);	
				//aval = *(a + l*m + i);
				aval = *(ap + i);
				ax = svdup_f64(aval);
				clx = svmla_f64_m(npred,clx, bx, ax);
				svst1_f64(npred, c + i*ldc, clx);
			}

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
