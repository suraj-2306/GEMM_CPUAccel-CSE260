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
void inline bl_dgemm_ukr(int k, int m, int n, double *restrict a, double *restrict b, double *restrict c,
		unsigned long long ldc, aux_t *restrict data) {
	//svbool_t npred = svwhilelt_b64_u64(0,n);
	//int i,j,l;
	//register svfloat64_t bx;
	//register svfloat64_t ax0; register svfloat64_t clx0; register float64_t aval0;
	//register svfloat64_t ax1; register svfloat64_t clx1; register float64_t aval1;
	//register svfloat64_t ax2; register svfloat64_t clx2; register float64_t aval2;
	//register svfloat64_t ax3; register svfloat64_t clx3; register float64_t aval3;
	//register svfloat64_t ax4; register svfloat64_t clx4; register float64_t aval4;
	//register svfloat64_t ax5; register svfloat64_t clx5; register float64_t aval5;
	//register svfloat64_t ax6; register svfloat64_t clx6; register float64_t aval6;
	//register svfloat64_t ax7; register svfloat64_t clx7; register float64_t aval7;

	//register svfloat64_t ax; register svfloat64_t clx; register float64_t aval;

	//double *restrict ap;
	//clx0 = svld1_f64(npred,c + 0*ldc);  
	//clx1 = svld1_f64(npred,c + 1*ldc);  
	//clx2 = svld1_f64(npred,c + 2*ldc);  
	//clx3 = svld1_f64(npred,c + 3*ldc);  
	//clx4 = svld1_f64(npred,c + 4*ldc);  
	//clx5 = svld1_f64(npred,c + 5*ldc);  
	//clx6 = svld1_f64(npred,c + 6*ldc);  
	//clx7 = svld1_f64(npred,c + 7*ldc);  

	//for (l = 0; l < k; l+=1) {
	//	bx = svld1_f64(svptrue_b64(), b + l*n );
	//	ap = a + l*m;
	//	if( m>>3 == 1) // Only work for MR == 8
	//	{
	//		aval0 = (*(ap + 0)); ax0 = svdup_f64(aval0); clx0 = svmla_f64_m(npred,clx0, bx, ax0);
	//		aval1 = (*(ap + 1)); ax1 = svdup_f64(aval1); clx1 = svmla_f64_m(npred,clx1, bx, ax1);
	//		aval2 = (*(ap + 2)); ax2 = svdup_f64(aval2); clx2 = svmla_f64_m(npred,clx2, bx, ax2);
	//		aval3 = (*(ap + 3)); ax3 = svdup_f64(aval3); clx3 = svmla_f64_m(npred,clx3, bx, ax3);
	//		aval4 = (*(ap + 4)); ax4 = svdup_f64(aval4); clx4 = svmla_f64_m(npred,clx4, bx, ax4);
	//		aval5 = (*(ap + 5)); ax5 = svdup_f64(aval5); clx5 = svmla_f64_m(npred,clx5, bx, ax5);
	//		aval6 = (*(ap + 6)); ax6 = svdup_f64(aval6); clx6 = svmla_f64_m(npred,clx6, bx, ax6);
	//		aval7 = (*(ap + 7)); ax7 = svdup_f64(aval7); clx7 = svmla_f64_m(npred,clx7, bx, ax7);

	//	}
	//}

	//svst1_f64(npred, c + 0*ldc , clx0);
	//svst1_f64(npred, c + 1*ldc , clx1);
	//svst1_f64(npred, c + 2*ldc , clx2);
	//svst1_f64(npred, c + 3*ldc , clx3);
	//svst1_f64(npred, c + 4*ldc , clx4);
	//svst1_f64(npred, c + 5*ldc , clx5);
	//svst1_f64(npred, c + 6*ldc , clx6);
	//svst1_f64(npred, c + 7*ldc , clx7);

	//for (l = 0; l < k; l+=1) {
	//	bx = svld1_f64(svptrue_b64(), b + l*n );
	//	ap = a + l*m;

	//	if( m>>3 == 0) // Only work for MR == 8
	//	{
	//		for (i = 0; i < m; ++i) {
	//			clx = svld1_f64(npred,c + i*ldc);	
	//			aval = *(ap + i);
	//			ax = svdup_f64(aval);
	//			clx = svmla_f64_m(npred,clx, bx, ax);
	//			svst1_f64(npred, c + i*ldc, clx);
	//		}

	//	}
	//}
	svbool_t npred = svwhilelt_b64_u64(0,n);
	int i,j,l;
	register svfloat64_t bx;
	register svfloat64_t ax0; register svfloat64_t clx0; register float64_t aval0;
	register svfloat64_t ax1; register svfloat64_t clx1; register float64_t aval1;
	register svfloat64_t ax2; register svfloat64_t clx2; register float64_t aval2;
	register svfloat64_t ax3; register svfloat64_t clx3; register float64_t aval3;
	register svfloat64_t ax4; register svfloat64_t clx4; register float64_t aval4;
	register svfloat64_t ax5; register svfloat64_t clx5; register float64_t aval5;
	register svfloat64_t ax6; register svfloat64_t clx6; register float64_t aval6;
	register svfloat64_t ax7; register svfloat64_t clx7; register float64_t aval7;
	register svfloat64_t ax8; register svfloat64_t clx8; register float64_t aval8;
	register svfloat64_t ax9; register svfloat64_t clx9; register float64_t aval9;
	register svfloat64_t ax10; register svfloat64_t clx10; register float64_t aval10;
	register svfloat64_t ax11; register svfloat64_t clx11; register float64_t aval11;
	register svfloat64_t ax12; register svfloat64_t clx12; register float64_t aval12;
	register svfloat64_t ax13; register svfloat64_t clx13; register float64_t aval13;
	register svfloat64_t ax14; register svfloat64_t clx14; register float64_t aval14;
	register svfloat64_t ax15; register svfloat64_t clx15; register float64_t aval15;
	register svfloat64_t ax; register svfloat64_t clx; register float64_t aval; //Residual
	double *restrict ap;
	clx0 = svld1_f64(npred,c + 0*ldc);
	clx1 = svld1_f64(npred,c + 1*ldc);
	clx2 = svld1_f64(npred,c + 2*ldc);
	clx3 = svld1_f64(npred,c + 3*ldc);
	clx4 = svld1_f64(npred,c + 4*ldc);
	clx5 = svld1_f64(npred,c + 5*ldc);
	clx6 = svld1_f64(npred,c + 6*ldc);
	clx7 = svld1_f64(npred,c + 7*ldc);
	clx8 = svld1_f64(npred,c + 8*ldc);
	clx9 = svld1_f64(npred,c + 9*ldc);
	clx10 = svld1_f64(npred,c + 10*ldc);
	clx11 = svld1_f64(npred,c + 11*ldc);
	clx12 = svld1_f64(npred,c + 12*ldc);
	clx13 = svld1_f64(npred,c + 13*ldc);
	clx14 = svld1_f64(npred,c + 14*ldc);
	clx15 = svld1_f64(npred,c + 15*ldc);
	for (l = 0; l < k; l+=1) {
		bx = svld1_f64(svptrue_b64(), b + l*n );
		ap = a + l*m;
		if( m>>4 == 1) // Only work for MR == 8
		{
			aval0 = (*(ap + 0)); ax0 = svdup_f64(aval0); clx0 = svmla_f64_m(npred,clx0, bx, ax0);
			aval1 = (*(ap + 1)); ax1 = svdup_f64(aval1); clx1 = svmla_f64_m(npred,clx1, bx, ax1);
			aval2 = (*(ap + 2)); ax2 = svdup_f64(aval2); clx2 = svmla_f64_m(npred,clx2, bx, ax2);
			aval3 = (*(ap + 3)); ax3 = svdup_f64(aval3); clx3 = svmla_f64_m(npred,clx3, bx, ax3);
			aval4 = (*(ap + 4)); ax4 = svdup_f64(aval4); clx4 = svmla_f64_m(npred,clx4, bx, ax4);
			aval5 = (*(ap + 5)); ax5 = svdup_f64(aval5); clx5 = svmla_f64_m(npred,clx5, bx, ax5);
			aval6 = (*(ap + 6)); ax6 = svdup_f64(aval6); clx6 = svmla_f64_m(npred,clx6, bx, ax6);
			aval7 = (*(ap + 7)); ax7 = svdup_f64(aval7); clx7 = svmla_f64_m(npred,clx7, bx, ax7);
			aval8 = (*(ap + 8)); ax8 = svdup_f64(aval8); clx8 = svmla_f64_m(npred,clx8, bx, ax8);
			aval9 = (*(ap + 9)); ax9 = svdup_f64(aval9); clx9 = svmla_f64_m(npred,clx9, bx, ax9);
			aval10 = (*(ap + 10)); ax10 = svdup_f64(aval10); clx10 = svmla_f64_m(npred,clx10, bx, ax10);
			aval11 = (*(ap + 11)); ax11 = svdup_f64(aval11); clx11 = svmla_f64_m(npred,clx11, bx, ax11);
			aval12 = (*(ap + 12)); ax12 = svdup_f64(aval12); clx12 = svmla_f64_m(npred,clx12, bx, ax12);
			aval13 = (*(ap + 13)); ax13 = svdup_f64(aval13); clx13 = svmla_f64_m(npred,clx13, bx, ax13);
			aval14 = (*(ap + 14)); ax14 = svdup_f64(aval14); clx14 = svmla_f64_m(npred,clx14, bx, ax14);
			aval15 = (*(ap + 15)); ax15 = svdup_f64(aval15); clx15 = svmla_f64_m(npred,clx15, bx, ax15);
		}
	}
	svst1_f64(npred, c + 0*ldc , clx0);
	svst1_f64(npred, c + 1*ldc , clx1);
	svst1_f64(npred, c + 2*ldc , clx2);
	svst1_f64(npred, c + 3*ldc , clx3);
	svst1_f64(npred, c + 4*ldc , clx4);
	svst1_f64(npred, c + 5*ldc , clx5);
	svst1_f64(npred, c + 6*ldc , clx6);
	svst1_f64(npred, c + 7*ldc , clx7);
	svst1_f64(npred, c + 8*ldc , clx8);
	svst1_f64(npred, c + 9*ldc , clx9);
	svst1_f64(npred, c + 10*ldc , clx10);
	svst1_f64(npred, c + 11*ldc , clx11);
	svst1_f64(npred, c + 12*ldc , clx12);
	svst1_f64(npred, c + 13*ldc , clx13);
	svst1_f64(npred, c + 14*ldc , clx14);
	svst1_f64(npred, c + 15*ldc , clx15);
	for (l = 0; l < k; l+=1) {
		bx = svld1_f64(svptrue_b64(), b + l*n );
		ap = a + l*m;

		if( m>>4 == 0) // Only work for MR == 8
		{
			for (i = 0; i < m; ++i) {
				clx = svld1_f64(npred,c + i*ldc);
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
