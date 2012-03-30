#include <cula_lapack.h>

int lud_cula(float *matrix, int size, int *ipiv){
  //culaStatus culaDeviceSgetrf(int m, int n, culaDeviceFloat* a, int lda, culaDeviceInt* ipiv);
  return culaDeviceSgetrf(size, size, matrix, size, ipiv);
}
