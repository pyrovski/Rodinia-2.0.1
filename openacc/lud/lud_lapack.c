#include <clapack.h>

int lud_lapack(float *matrix, int size, int *ipiv){
  return clapack_sgetrf(CblasRowMajor, size, size, matrix, size, ipiv);
}
