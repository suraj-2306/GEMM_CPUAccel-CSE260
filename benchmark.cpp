/*
 *  Driver code for Matrix Multplication
 *  Provided by Jim Demmel at UC Berkeley
 *  Modified by Scott B. Baden at UC San Diego to
 *      Enable user to select one problem size only via the -n option
 *      Support CBLAS interface
 */

#include <stdlib.h> // For: exit, random, malloc, free, NULL, EXIT_FAILURE
#include <stdio.h>  // For: perror
#include <string.h> // For: memset

#include <float.h>  // For: DBL_EPSILON
#include <math.h>   // For: fabs
#include "blislab/bl_dgemm.h" // For: bl_malloc_aligned

#ifdef USE_MKL
#include "mkl.h"
#else
#include "cblas.h"
#endif

extern "C" {
  void square_dgemm(int, double*, double*, double*);
  double wall_time();
}


void cmdLine(int argc, char *argv[], int& n, int& noCheck, int& identDebug, int& genDATA);

/* reference_dgemm wraps a call to the BLAS-3 routine DGEMM, */
/* via the standard FORTRAN interface - hence the reference semantics. */ 
void reference_dgemm (int N, double Alpha, double* A, double* B, double* C)
{
  const double Beta  = 1.0;
  const int M = N, K=N;
  const int LDA = N, LDB = N, LDC = N;
  const enum CBLAS_TRANSPOSE transA = CblasNoTrans;
  const enum CBLAS_TRANSPOSE transB = CblasNoTrans;
  /* Don't change this call */
  cblas_dgemm( CblasRowMajor, transA, transB, M, N, K,
               Alpha, A, LDA, B, LDB, Beta, C, LDC );
}   

/* Your function must have the following signature: */
extern const char* dgemm_desc;
extern void square_dgemm (int, double*, double*, double*);

extern double wall_time();


#include "debugMat.h"

void Fail (const char* message)
{
  perror (message);
  exit (EXIT_FAILURE);
}


void fill (double* p, int n)
{
  long int Rmax   = RAND_MAX;
  long int Rmax_2 = Rmax >> 1;
  long int RM     =  Rmax_2 + 1;
  for (int i = 0; i < n; ++i){
    long int r = rand();   // Uniformly distributed ints over [0,RAND_MAX]
                             // Typical value of RAND_MAX: 2^31 - 1
    long int R = r - RM;
    p[i] = (double) R / (double) RM; // Uniformly distributed over [-1, 1]
  }
}

void absolute_value (double *p, int n)
{
  for (int i = 0; i < n; ++i)
    p[i] = fabs (p[i]);
}

/* The benchmarking program */
int main (int argc, char **argv)
{
  
  /* We can pick just one size with the -n flag */
  int n0;
  int noCheck;
  int identDebug;
  int genDATA;
  cmdLine(argc, argv, n0, noCheck, identDebug, genDATA);
  
  if (!genDATA)
    printf ("Description:\t%s\n\n", dgemm_desc);

  /* Test sizes should highlight performance d`ips at multiples of certain powers-of-two */

  int test_sizes[] = 

  /* Multiples-of-32, +/- 1. Currently uncommented. */
  //{31,32,33,63,64,65,95,96,97,127,128,129,159,160,161,191,192,193,223,224,
  // 225,255,256,257,287,288,289,319,320,321,351,352,353,383,384,385,415,416,
  // 417,447,448,449,479,480,481,511,512,513,543,544,545,575,576,577,607,608,
  // 609,639,640,641,671,672,673,703,704,705,735,736,737,767,768,769,799,800,
  // 801,831,832,833,863,864,865,895,896,897,927,928,929,959,960,961,991,992,
  // 993,1023,1024,1025}; 

  //{513,543,544,545,575,576,577,607,608,
  // 609,639,640,641,671,672,673,703,704,705,735,736,737,767,768,769,799,800,
  // 801,831,832,833,863,864,865,895,896,897,927,928,929,959,960,961,991,992,
  // 993,1023,1024,1025}; 

   //{512, 514, 516, 518, 520, 522, 524, 526, 528, 530, 532, 534, 536, 538, 540, 542, 544, 546, 548, 550, 552, 554, 556, 558, 560, 562, 564, 566, 568, 570, 572, 574, 576, 578, 580, 
   //582, 584, 586, 588, 590, 592, 594, 596, 598, 600, 602, 604, 606, 608, 610, 612, 614, 616, 618, 620, 622, 624, 626, 628, 630, 632, 634, 636, 638, 640, 642, 644, 646, 648, 650, 
   //652, 654, 656, 658, 660, 662, 664, 666, 668, 670, 672, 674, 676, 678, 680, 682, 684, 686, 688, 690, 692, 694, 696, 698, 700, 702, 704, 706, 708, 710, 712, 714, 716, 718, 720, 
   //722, 724, 726, 728, 730, 732, 734, 736, 738, 740, 742, 744, 746, 748, 750, 752, 754, 756, 758, 760, 762, 764, 766, 768, 770, 772, 774, 776, 778, 780, 782, 784, 786, 788, 790, 
   //792, 794, 796, 798, 800, 802, 804, 806, 808, 810, 812, 814, 816, 818, 820, 822, 824, 826, 828, 830, 832, 834, 836, 838, 840, 842, 844, 846, 848, 850, 852, 854, 856, 858, 860, 
   //862, 864, 866, 868, 870, 872, 874, 876, 878, 880, 882, 884, 886, 888, 890, 892, 894, 896, 898, 900, 902, 904, 906, 908, 910, 912, 914, 916, 918, 920, 922, 924, 926, 928, 930, 
   //932, 934, 936, 938, 940, 942, 944, 946, 948, 950, 952, 954, 956, 958, 960, 962, 964, 966, 968, 970, 972, 974, 976, 978, 980, 982, 984, 986, 988, 990, 992, 994, 996, 998, 1000, 
   //1002, 1004, 1006, 1008, 1010, 1012, 1014, 1016, 1018, 1020, 1022, 1024};

   {32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110,
   112, 114, 116, 118, 120, 122, 124, 126, 128, 130, 132, 134, 136, 138, 140, 142, 144, 146, 148, 150, 152, 154, 156, 158, 160, 162, 164, 166, 168, 170, 172, 174, 176, 178, 180, 
   182, 184, 186, 188, 190, 192, 194, 196, 198, 200, 202, 204, 206, 208, 210, 212, 214, 216, 218, 220, 222, 224, 226, 228, 230, 232, 234, 236, 238, 240, 242, 244, 246, 248, 250, 
   252, 254, 256, 258, 260, 262, 264, 266, 268, 270, 272, 274, 276, 278, 280, 282, 284, 286, 288, 290, 292, 294, 296, 298, 300, 302, 304, 306, 308, 310, 312, 314, 316, 318, 320, 
   322, 324, 326, 328, 330, 332, 334, 336, 338, 340, 342, 344, 346, 348, 350, 352, 354, 356, 358, 360, 362, 364, 366, 368, 370, 372, 374, 376, 378, 380, 382, 384, 386, 388, 390, 
   392, 394, 396, 398, 400, 402, 404, 406, 408, 410, 412, 414, 416, 418, 420, 422, 424, 426, 428, 430, 432, 434, 436, 438, 440, 442, 444, 446, 448, 450, 452, 454, 456, 458, 460, 
   462, 464, 466, 468, 470, 472, 474, 476, 478, 480, 482, 484, 486, 488, 490, 492, 494, 496, 498, 500, 502, 504, 506, 508, 510, 512, 514, 516, 518, 520, 522, 524, 526, 528, 530, 
   532, 534, 536, 538, 540, 542, 544, 546, 548, 550, 552, 554, 556, 558, 560, 562, 564, 566, 568, 570, 572, 574, 576, 578, 580, 582, 584, 586, 588, 590, 592, 594, 596, 598, 600, 
   602, 604, 606, 608, 610, 612, 614, 616, 618, 620, 622, 624, 626, 628, 630, 632, 634, 636, 638, 640, 642, 644, 646, 648, 650, 652, 654, 656, 658, 660, 662, 664, 666, 668, 670, 
   672, 674, 676, 678, 680, 682, 684, 686, 688, 690, 692, 694, 696, 698, 700, 702, 704, 706, 708, 710, 712, 714, 716, 718, 720, 722, 724, 726, 728, 730, 732, 734, 736, 738, 740, 
   742, 744, 746, 748, 750, 752, 754, 756, 758, 760, 762, 764, 766, 768, 770, 772, 774, 776, 778, 780, 782, 784, 786, 788, 790, 792, 794, 796, 798, 800, 802, 804, 806, 808, 810, 
   812, 814, 816, 818, 820, 822, 824, 826, 828, 830, 832, 834, 836, 838, 840, 842, 844, 846, 848, 850, 852, 854, 856, 858, 860, 862, 864, 866, 868, 870, 872, 874, 876, 878, 880, 
   882, 884, 886, 888, 890, 892, 894, 896, 898, 900, 902, 904, 906, 908, 910, 912, 914, 916, 918, 920, 922, 924, 926, 928, 930, 932, 934, 936, 938, 940, 942, 944, 946, 948, 950, 
   952, 954, 956, 958, 960, 962, 964, 966, 968, 970, 972, 974, 976, 978, 980, 982, 984, 986, 988, 990, 992, 994, 996, 998, 1000, 1002, 1004, 1006, 1008, 1010, 1012, 1014, 1016, 1018, 1020, 1022, 1024};

   //{32,64,128,256,511,512,513,1023,1024,1025,2047,2048};

  /* Multiples-of-32, +/- 1. Currently uncommented. Large  only */
  /*  {511,512,513,543,544,545,575,576,577,607,608,609,639,640,641,671,672,673,703,704,705,735,736,737,767,768,769,799,800,801,831,832,833,863,864,865,895,896,897,927,928,929,959,960,961,991,992,993,1023,1024,1025};  */
/*
  {31,32,33,63,64,65,95,96,97,127,128,129,159,160,161,191,192,193,223,224,225,255,256,257,287,288,289,319,320,321,351,352,353,383,384,385,415,416,417,447,448,449,479,480,481,511,512,513,543,544,545,575,576,577,607,608,609,639,640,641,671,672,673,703,704,705,735,736,737,767,768,769};
  */

  /* A representative subset of the first list. Currently commented. */ 
  /*
  { 31, 32, 96, 97, 127, 128, 129, 191, 192, 229, 255, 256, 257,
    319, 320, 321, 417, 479, 480, 511, 512, 639, 640, 767, 768, 769 };
    */

  /* N= 31 through 768+-1 in steps of 16 */
  /*
  { 31, 32, 33, 39, 40, 41, 47, 48, 49, 55, 56, 57, 63, 64, 65, 71, 72, 73, 79, 80, 81, 87, 88, 89, 95, 96, 97, 103, 104, 105, 111, 112, 113, 119, 120, 121, 127, 128, 129, 135, 136, 137, 143, 144, 145, 151, 152, 153, 159, 160, 161, 167, 168, 169, 175, 176, 177, 183, 184, 185, 191, 192, 193, 199, 200, 201, 207, 208, 209, 215, 216, 217, 223, 224, 225, 231, 232, 233, 239, 240, 241, 247, 248, 249, 255, 256, 257, 263, 264, 265, 271, 272, 273, 279, 280, 281, 287, 288, 289, 295, 296, 297, 303, 304, 305, 311, 312, 313, 319, 320, 321, 327, 328, 329, 335, 336, 337, 343, 344, 345, 351, 352, 353, 359, 360, 361, 367, 368, 369, 375, 376, 377, 383, 384, 385, 391, 392, 393, 399, 400, 401, 407, 408, 409, 415, 416, 417, 423, 424, 425, 431, 432, 433, 439, 440, 441, 447, 448, 449, 455, 456, 457, 463, 464, 465, 471, 472, 473, 479, 480, 481, 487, 488, 489, 495, 496, 497, 503, 504, 505, 511, 512, 513, 519, 520, 521, 527, 528, 529, 535, 536, 537, 543, 544, 545, 551, 552, 553, 559, 560, 561, 567, 568, 569, 575, 576, 577, 583, 584, 585, 591, 592, 593, 599, 600, 601, 607, 608, 609, 615, 616, 617, 623, 624, 625, 631, 632, 633, 639, 640, 641, 647, 648, 649, 655, 656, 657, 663, 664, 665, 671, 672, 673, 679, 680, 681, 687, 688, 689, 695, 696, 697, 703, 704, 705, 711, 712, 713, 719, 720, 721, 727, 728, 729, 735, 736, 737, 743, 744, 745, 751, 752, 753, 759, 760, 761, 767, 768, 769 }; */

  int nsizes = sizeof(test_sizes)/sizeof(test_sizes[0]);

  /* assume last size is also the largest size */
  int nmax = test_sizes[nsizes-1];

  if (n0){
    nmax = n0;
    test_sizes[0] = n0;
  }

  /* allocate memory for all problems */
  double *A;
  if ((A = bl_malloc_aligned( nmax, nmax, sizeof(double))) == NULL){
    Fail ("Failed to allocate matrix A");
  }
  double *B;
  if ((B = bl_malloc_aligned( nmax, nmax, sizeof(double))) == NULL){
    Fail ("Failed to allocate matrix B");
  }
  double *C;
  if ((C = bl_malloc_aligned( nmax, nmax, sizeof(double))) == NULL){
    Fail ("Failed to allocated matrix C");
  }
  int sizes = sizeof(test_sizes)/sizeof(test_sizes[0]);
  if (n0)
    sizes = 1;

  /* For each test size */
  double sumLns  = 0.0;  // sum of log  performance
  for (int isize = 0; isize < sizes; ++isize)
  {
    /* Create and fill 3 random matrices A,B,C*/
    int n = test_sizes[isize];


    if (identDebug){
      identMat(n, B);
      seqMat(n, n, A);
      setMat(n, n, C, 0.0);
      printMat(n, n, "A:", A);
      printMat(n, n, "B:", B);
      square_dgemm(n, A, B, C);
      printMat(n, n, "C=", C);
      return -1;
    }else{
      fill (A, n*n);
      fill (B, n*n);
      fill (C, n*n);
    }

    /* Measure performance (in Gflops/s). */

    /* Time a "sufficiently long" sequence of calls to reduce noise */
    double Gflops_s, seconds = -1.0;
    double timeout = 0.1; // "sufficiently long" := at least 1/10 second.

    for (int n_iterations = 1; seconds < timeout; n_iterations *= 2) 
    {
      /* Warm-up */
      square_dgemm (n, A, B, C);

      /* Benchmark n_iterations runs of square_dgemm */
      seconds = -wall_time();
      for (int it = 0; it < n_iterations; ++it)
	square_dgemm (n, A, B, C);
      seconds += wall_time();

      /*  compute Mflop/s rate */
      Gflops_s = 2.e-9 * n_iterations * n * n * n / seconds;
    }


    if(genDATA){
      printf ("%d\t%.3g\n", n, Gflops_s);
    }
    else{
      printf ("Size: %d\tGflop/s: %.3g\n", n, Gflops_s);
    }
    
    sumLns +=  log(Gflops_s);

    // verify result
    //
    if (!noCheck){
      /* Ensure that error does not exceed the theoretical error bound. */

      /* C := A * B, computed with square_dgemm */
      memset (C, 0, n * n * sizeof(double));
      square_dgemm (n, A, B, C);
      
      /* Do not explicitly check that A and B were unmodified on square_dgemm exit
       *  - if they were, the following will most likely detect it:   
       * C := C - A * B, computed with reference_dgemm */
      reference_dgemm(n, -1., A, B, C);
      
      /* A := |A|, B := |B|, C := |C| */
      absolute_value (A, n * n);
      absolute_value (B, n * n);
      absolute_value (C, n * n);
      
      /* C := |C| - 3 * e_mach * n * |A| * |B|, computed with reference_dgemm */ 
      reference_dgemm (n, -3.*DBL_EPSILON*n, A, B, C);
      
      /* If any element in C is positive, then something went wrong in square_dgemm */
      for (int i = 0; i < n * n; ++i)
	if (C[i] > 0){
	  fprintf(stderr, "error at index C[%d] = %5.2f\n", i, C[i]);
	  Fail("*** FAILURE *** Error in matrix multiply exceeds componentwise error bounds.\n" );
	}
    }
  }
  if(!genDATA)
    printf("GeoMean  = %4.2f\n",  exp(1/(double)sizes * sumLns));

  free (A);
  free (B);
  free (C);
  return 0;
}
