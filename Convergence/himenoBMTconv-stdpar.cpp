/********************************************************************

 This benchmark test program is measuring a cpu performance
 of floating point operation by a Poisson equation solver.

 If you have any question, please ask me via email.
 written by Ryutaro HIMENO, November 26, 2001.
 Version 3.0
 ----------------------------------------------
 Ryutaro Himeno, Dr. of Eng.
 Head of Computer Information Division,
 RIKEN (The Institute of Pysical and Chemical Research)
 Email : himeno@postman.riken.go.jp
 ---------------------------------------------------------------
 You can adjust the size of this benchmark code to fit your target
 computer. In that case, please chose following sets of
 [mimax][mjmax][mkmax]:
 small : 33,33,65
 small : 65,65,129
 midium: 129,129,257
 large : 257,257,513
 ext.large: 513,513,1025
 This program is to measure a computer performance in MFLOPS
 by using a kernel which appears in a linear solver of pressure
 Poisson eq. which appears in an incompressible Navier-Stokes solver.
 A point-Jacobi method is employed in this solver as this method can
 be easyly vectrized and be parallelized.
 ------------------
 Finite-difference method, curvilinear coodinate system
 Vectorizable and parallelizable on each grid point
 No. of grid points : imax x jmax x kmax including boundaries
 ------------------
 A,B,C:coefficient matrix, wrk1: source term of Poisson equation
 wrk2 : working area, OMEGA : relaxation parameter
 BND:control variable for boundaries and objects ( = 0 or 1)
 P: pressure
********************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>

#include <algorithm>
#include <numeric>
#include <execution>
#include <ranges>
#include <memory>
#include <string>

#define MR(mt, n, r, c, d) \
  mt->m[(n) * mt->mrows * mt->mcols * mt->mdeps + (r) * mt->mcols * mt->mdeps + (c) * mt->mdeps + (d)]

/* Macro for referencing a pointer to a member of a Mat structure */
#define MP(mt, n, r, c, d) \
  mt.m[(n) * mt.mrows * mt.mcols * mt.mdeps + (r) * mt.mcols * mt.mdeps + (c) * mt.mdeps + (d)]

#define OMEGA 0.8f
#define EPS 1e-5
#define MAX_ITER 200000

struct Mat {
  float* m;
  int mnums;
  int mrows;
  int mcols;
  int mdeps;
};

/* prototypes */
typedef struct Mat Matrix;

int newMat(Matrix* Mat, int mnums, int mrows, int mcols, int mdeps);
void clearMat(Matrix* Mat);
void set_param(int i[], char* size);
void mat_set(Matrix* Mat, int l, float z);
void mat_set_init(Matrix* Mat);
float jacobi(int* n, int nmax, Matrix* M1, Matrix* M2, Matrix* M3, Matrix* M4, Matrix* M5, Matrix* M6, Matrix* M7);
float jacobi_for_each_n(int* n, int nmax, Matrix* M1, Matrix* M2, Matrix* M3, Matrix* M4, Matrix* M5, Matrix* M6, Matrix* M7);
double fflop(int, int, int);
double mflops(int, double, double);
double gflops(int, double, double);
double second();

int main(int argc, char* argv[]) {
  int nn, nnmax;
  int imax, jmax, kmax, mimax, mjmax, mkmax, msize[3];
  float gosa;
  double cpu0, cpu1, cpu, flop;
  char size[10];
  Matrix a, b, c, p, bnd, wrk1, wrk2;

  if (argc >= 2) {
    strcpy(size, argv[1]);
  } else {
    printf("For example: \n");
    printf(" Grid-size= XS (32x32x64)\n");
    printf("\t    S  (64x64x128)\n");
    printf("\t    M  (128x128x256)\n");
    printf("\t    L  (256x256x512)\n");
    printf("\t    XL (512x512x1024)\n\n");
    printf("Grid-size = ");
    scanf("%s", size);
    printf("\n");
  }

  set_param(msize, size);

  mimax = msize[0];
  mjmax = msize[1];
  mkmax = msize[2];
  imax = mimax - 1;
  jmax = mjmax - 1;
  kmax = mkmax - 1;

  printf("mimax = %d mjmax = %d mkmax = %d\n", mimax, mjmax, mkmax);
  printf("imax = %d jmax = %d kmax =%d\n", imax, jmax, kmax);

  /*
   *    Initializing matrixes
   */
  newMat(&p, 1, mimax, mjmax, mkmax);
  newMat(&bnd, 1, mimax, mjmax, mkmax);
  newMat(&wrk1, 1, mimax, mjmax, mkmax);
  newMat(&wrk2, 1, mimax, mjmax, mkmax);
  newMat(&a, 4, mimax, mjmax, mkmax);
  newMat(&b, 3, mimax, mjmax, mkmax);
  newMat(&c, 3, mimax, mjmax, mkmax);

  mat_set_init(&p);
  mat_set(&bnd, 0, 1.0);
  mat_set(&wrk1, 0, 0.0);
  mat_set(&wrk2, 0, 0.0);
  mat_set(&a, 0, 1.0);
  mat_set(&a, 1, 1.0);
  mat_set(&a, 2, 1.0);
  mat_set(&a, 3, 1.0 / 6.0);
  mat_set(&b, 0, 0.0);
  mat_set(&b, 1, 0.0);
  mat_set(&b, 2, 0.0);
  mat_set(&c, 0, 1.0);
  mat_set(&c, 1, 1.0);
  mat_set(&c, 2, 1.0);

  /* Start rehearsal measurement */
  flop = fflop(imax, jmax, kmax);
  nnmax = 3;
  nn = 0;
  printf(" Start rehearsal measurement process.\n");
  printf(" Measure the performance in %d times.\n\n", nnmax);

  cpu0 = second();
  gosa = jacobi(&nn, nnmax, &a, &b, &c, &p, &bnd, &wrk1, &wrk2);
  // gosa = jacobi_for_each_n(&nn, nnmax, &a, &b, &c, &p, &bnd, &wrk1, &wrk2);
  cpu1 = second();
  cpu = cpu1 - cpu0;

  printf(" GFLOPS: %f time(s): %f error: %e\n\n", gflops(nn, cpu, flop), cpu, gosa);

  /* Start the actual mesurement */
  nnmax = MAX_ITER;
  nn = 0;

  printf(" Now, start the actual measurement process.\n");
  printf(" Maximum nunber of iterations is %d.\n", nnmax);

  cpu0 = second();
  gosa = jacobi(&nn, nnmax, &a, &b, &c, &p, &bnd, &wrk1, &wrk2);
  // gosa = jacobi_for_each_n(&nn, nnmax, &a, &b, &c, &p, &bnd, &wrk1, &wrk2);
  cpu1 = second();
  cpu = cpu1 - cpu0;

  printf(" Loop executed for %d times\n", nn);
  printf(" Gosa : %e \n", gosa);
  printf(" GFLOPS measured : %f\tcpu : %f\n", gflops(nn, cpu, flop), cpu);
  printf(" Score based on Pentium III 600MHz using Fortran 77: %f\n", mflops(nn, cpu, flop) / 82.84);

  /* Matrix free */
  clearMat(&p);
  clearMat(&bnd);
  clearMat(&wrk1);
  clearMat(&wrk2);
  clearMat(&a);
  clearMat(&b);
  clearMat(&c);

  return (0);
}

double fflop(int mx, int my, int mz) {
  return ((double)(mz - 2) * (double)(my - 2) * (double)(mx - 2) * 34.0);
}

double mflops(int nn, double cpu, double flop) {
  return (flop / cpu * 1.e-6 * (double)nn);
}

double gflops(int nn, double cpu, double flop) {
  return (flop / cpu * 1.e-9 * (double)nn);
}

void set_param(int is[], char* size) {
  if (!strcmp(size, "XS") || !strcmp(size, "xs")) {
    is[0] = 32;
    is[1] = 32;
    is[2] = 64;
    return;
  }
  if (!strcmp(size, "S") || !strcmp(size, "s")) {
    is[0] = 64;
    is[1] = 64;
    is[2] = 128;
    return;
  }
  if (!strcmp(size, "M") || !strcmp(size, "m")) {
    is[0] = 128;
    is[1] = 128;
    is[2] = 256;
    return;
  }
  if (!strcmp(size, "L") || !strcmp(size, "l")) {
    is[0] = 256;
    is[1] = 256;
    is[2] = 512;
    return;
  }
  if (!strcmp(size, "XL") || !strcmp(size, "xl")) {
    is[0] = 512;
    is[1] = 512;
    is[2] = 1024;
    return;
  } else {
    printf("Invalid input character !!\n");
    exit(6);
  }
}

int newMat(Matrix* Mat, int mnums, int mrows, int mcols, int mdeps) {
  Mat->mnums = mnums;
  Mat->mrows = mrows;
  Mat->mcols = mcols;
  Mat->mdeps = mdeps;
  Mat->m = nullptr;
  // Mat->m = (float*)malloc(mnums * mrows * mcols * mdeps * sizeof(float));
  Mat->m = new float[mnums * mrows * mcols * mdeps];

  return (Mat->m != nullptr) ? 1 : 0;
}

void clearMat(Matrix* Mat) {
  // if (Mat->m) free(Mat->m);
  if (Mat->m) delete[] Mat->m;
  Mat->m = nullptr;
  Mat->mnums = 0;
  Mat->mcols = 0;
  Mat->mrows = 0;
  Mat->mdeps = 0;
}

void mat_set(Matrix* Mat, int l, float val) {
  int i, j, k;

  for (i = 0; i < Mat->mrows; i++)
    for (j = 0; j < Mat->mcols; j++)
      for (k = 0; k < Mat->mdeps; k++) MR(Mat, l, i, j, k) = val;
}

void mat_set_init(Matrix* Mat) {
  int i, j, k;

  for (i = 0; i < Mat->mrows; i++)
    for (j = 0; j < Mat->mcols; j++)
      for (k = 0; k < Mat->mdeps; k++)
        MR(Mat, 0, i, j, k) =
            (float)(i * i) / (float)((Mat->mrows - 1) * (Mat->mrows - 1));
}

/**
 * Jacobi method changed to Standard Parallelism implementation. 
 * The transform_reduce() function implements simultaneous 
 * computation of the asymptotic expression and error aggregation.
*/
float jacobi(int* nn, int nnmax, Matrix* _a, Matrix* _b, Matrix* _c, Matrix* _p, Matrix* _bnd, Matrix* _wrk1, Matrix* _wrk2) {
  /* -stdpar=gpu to deal with structure pointer members that cannot be changed */
  Matrix a = *_a;
  Matrix b = *_b;
  Matrix c = *_c;
  Matrix p = *_p;
  Matrix bnd = *_bnd;
  Matrix wrk1 = *_wrk1;
  Matrix wrk2 = *_wrk2;

  float gosa;
  int n;
  int imax = p.mrows - 2;
  int jmax = p.mcols - 2;
  int kmax = p.mdeps - 2;
  int gridSize = imax * jmax * kmax;

  /* View of iterator for transform_reduce */
  std::ranges::iota_view idx_iota{0, gridSize};

  for (n = 1; n < nnmax; n++) {
    gosa = std::transform_reduce(std::execution::par_unseq, 
      idx_iota.begin(), idx_iota.end(), 0.0f, std::plus<float>(), 
      [=](int idx) {
        int i = idx / (jmax * kmax) + 1,
            j = (idx % (jmax * kmax)) / kmax + 1,
            k = idx % kmax + 1;

        float s0 = MP(a, 0, i, j, k) * MP(p, 0, i + 1, j, k) +
                  MP(a, 1, i, j, k) * MP(p, 0, i, j + 1, k) +
                  MP(a, 2, i, j, k) * MP(p, 0, i, j, k + 1) +
                  MP(b, 0, i, j, k) * (MP(p, 0, i + 1, j + 1, k) - MP(p, 0, i + 1, j - 1, k) -
                                        MP(p, 0, i - 1, j + 1, k) + MP(p, 0, i - 1, j - 1, k)) +
                  MP(b, 1, i, j, k) * (MP(p, 0, i, j + 1, k + 1) - MP(p, 0, i, j - 1, k + 1) -
                                        MP(p, 0, i, j + 1, k - 1) + MP(p, 0, i, j - 1, k - 1)) +
                  MP(b, 2, i, j, k) * (MP(p, 0, i + 1, j, k + 1) - MP(p, 0, i - 1, j, k + 1) -
                                        MP(p, 0, i + 1, j, k - 1) + MP(p, 0, i - 1, j, k - 1)) +
                  MP(c, 0, i, j, k) * MP(p, 0, i - 1, j, k) +
                  MP(c, 1, i, j, k) * MP(p, 0, i, j - 1, k) +
                  MP(c, 2, i, j, k) * MP(p, 0, i, j, k - 1) + MP(wrk1, 0, i, j, k);

        float ss = (s0 * MP(a, 3, i, j, k) - MP(p, 0, i, j, k)) * MP(bnd, 0, i, j, k);
        MP(wrk2, 0, i, j, k) = MP(p, 0, i, j, k) + OMEGA * ss;

        return ss * ss;
      });

    std::for_each_n(std::execution::par_unseq, std::ranges::views::iota(0).begin(), gridSize, 
      [=](int idx) {
        int i = idx / (jmax * kmax) + 1,
            j = (idx % (jmax * kmax)) / kmax + 1,
            k = idx % kmax + 1;
        MP(p, 0, i, j, k) = MP(wrk2, 0, i, j, k);
      });

    /* Convergence check */
    if (gosa < EPS) break;

  } /* end n loop */

  *nn = n;

  return (gosa);
}

/**
 * The error of each element is stored in an array, and the for_each_n() 
 * function is used to compute the gradual equation and store the error.
 * Finally, the reduce() function is used to compute the sum of the errors and return a value.
*/
float jacobi_for_each_n(int* nn, int nnmax, Matrix* _a, Matrix* _b, Matrix* _c, Matrix* _p, Matrix* _bnd, Matrix* _wrk1, Matrix* _wrk2) {
  /* -stdpar=gpu to deal with structure pointer members that cannot be changed */
  Matrix a = *_a;
  Matrix b = *_b;
  Matrix c = *_c;
  Matrix p = *_p;
  Matrix bnd = *_bnd;
  Matrix wrk1 = *_wrk1;
  Matrix wrk2 = *_wrk2;

  int n;
  int imax = p.mrows - 2;
  int jmax = p.mcols - 2;
  int kmax = p.mdeps - 2;
  int gridSize = imax * jmax * kmax;
  
  /* Declare array to store error for each step */
  float* gosa = new float[gridSize];
  float sum_gosa = 0.0f;

  for (n = 1; n <= nnmax; n++) {
    std::for_each_n(std::execution::par_unseq, std::ranges::views::iota(0).begin(), gridSize, 
      [=](int idx) {
        int i = idx / (jmax * kmax) + 1,
            j = (idx % (jmax * kmax)) / kmax + 1,
            k = idx % kmax + 1;

        float s0 = MP(a, 0, i, j, k) * MP(p, 0, i + 1, j, k) +
                  MP(a, 1, i, j, k) * MP(p, 0, i, j + 1, k) +
                  MP(a, 2, i, j, k) * MP(p, 0, i, j, k + 1) +
                  MP(b, 0, i, j, k) * (MP(p, 0, i + 1, j + 1, k) - MP(p, 0, i + 1, j - 1, k) -
                                        MP(p, 0, i - 1, j + 1, k) + MP(p, 0, i - 1, j - 1, k)) +
                  MP(b, 1, i, j, k) * (MP(p, 0, i, j + 1, k + 1) - MP(p, 0, i, j - 1, k + 1) -
                                        MP(p, 0, i, j + 1, k - 1) + MP(p, 0, i, j - 1, k - 1)) +
                  MP(b, 2, i, j, k) * (MP(p, 0, i + 1, j, k + 1) - MP(p, 0, i - 1, j, k + 1) -
                                        MP(p, 0, i + 1, j, k - 1) + MP(p, 0, i - 1, j, k - 1)) +
                  MP(c, 0, i, j, k) * MP(p, 0, i - 1, j, k) +
                  MP(c, 1, i, j, k) * MP(p, 0, i, j - 1, k) +
                  MP(c, 2, i, j, k) * MP(p, 0, i, j, k - 1) + MP(wrk1, 0, i, j, k);

        float ss = (s0 * MP(a, 3, i, j, k) - MP(p, 0, i, j, k)) * MP(bnd, 0, i, j, k);
        MP(wrk2, 0, i, j, k) = MP(p, 0, i, j, k) + OMEGA * ss;
        
        gosa[idx] = ss * ss;
      });

    std::for_each_n(std::execution::par_unseq, std::ranges::views::iota(0).begin(), gridSize, 
      [=](int idx) {
        int i = idx / (jmax * kmax) + 1,
            j = (idx % (jmax * kmax)) / kmax + 1,
            k = idx % kmax + 1;
        MP(p, 0, i, j, k) = MP(wrk2, 0, i, j, k);
      });
    
    /* Calculate the sum of errors */
    sum_gosa = std::reduce(std::execution::par_unseq, gosa, gosa + gridSize, 0.0f);

    /* Convergence check */
    if (sum_gosa < EPS) break;

  } /* end n loop */

  *nn = n;
  delete[] gosa;
  return (sum_gosa);
}

double second() {
  struct timeval tm;
  double t;

  static int base_sec = 0, base_usec = 0;

  gettimeofday(&tm, nullptr);

  if (base_sec == 0 && base_usec == 0) {
    base_sec = tm.tv_sec;
    base_usec = tm.tv_usec;
    t = 0.0;
  } else {
    t = (double)(tm.tv_sec - base_sec) +
        ((double)(tm.tv_usec - base_usec)) / 1.0e6;
  }

  return t;
}
