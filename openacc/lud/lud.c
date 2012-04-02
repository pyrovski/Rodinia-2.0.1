/*
 * =====================================================================================
 *
 *       Filename:  suite.c
 *
 *    Description:  The main wrapper for the suite
 *
 *        Version:  1.0
 *        Created:  10/22/2009 08:40:34 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Liang Wang (lw2aw), lw2aw@virginia.edu
 *        Company:  CS@UVa
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <assert.h>

#include <openacc.h>
#include <cula_lapack.h>

#include "common.h"

static int do_verify = 0;
int omp_num_threads = 1;

static struct option long_options[] = {
  /* name, has_arg, flag, val */
  {"input", 1, NULL, 'i'},
  {"size", 1, NULL, 's'},
  {"verify", 0, NULL, 'v'},
  {0,0,0,0}
};

extern void
lud_openacc(float *matrix, int matrix_dim);

extern void
lud_base(float *matrix, int matrix_dim);

void checkStatus(culaStatus status)
{   
  char buf[256];

  if(!status)
    return;

  culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
  printf("%s\n", buf);

  culaShutdown();
  exit(EXIT_FAILURE);
}

int
main ( int argc, char *argv[] )
{
  int matrix_dim = 32; /* default size */
  int status;
  int opt, option_index=0;
  func_ret_t ret;
  const char *input_file = NULL;
  float *matrix, *matrix_copy;
  stopwatch sw;
  int cpuonly = 0, lapack = 0, cula=0;

  while ((opt = getopt_long(argc, argv, "vs:i:n:cla", 
                            long_options, &option_index)) != -1 ) {
    switch(opt){
    case 'a':
      cula =1;
      break;
    case 'c':
      cpuonly = 1;
      break;
    case 'l':
      lapack = 1;
      break;
    case 'i':
      input_file = optarg;
      break;
    case 'n':
      omp_num_threads = atoi(optarg);
      break;
    case 'v':
      do_verify = 1;
      break;
    case 's':
      matrix_dim = atoi(optarg);
      fprintf(stderr, "Currently not supported, use -i instead\n");
      fprintf(stderr, "Usage: %s [-v] [-n no. of threads] [-s matrix_size|-i input_file]\n", argv[0]);
      exit(EXIT_FAILURE);
    case '?':
      fprintf(stderr, "invalid option\n");
      break;
    case ':':
      fprintf(stderr, "missing argument\n");
      break;
    default:
      fprintf(stderr, "Usage: %s [-v] [-n no. of threads] [-s matrix_size|-i input_file]\n",
	      argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  
  if ( (optind < argc) || (optind == 1)) {
    fprintf(stderr, "Usage: %s [-v] [-n no. of threads] [-s matrix_size|-i input_file]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  if (input_file) {
    printf("Reading matrix from file %s\n", input_file);
    ret = create_matrix_from_file(&matrix, input_file, &matrix_dim);
    if (ret != RET_SUCCESS) {
      matrix = NULL;
      fprintf(stderr, "error create matrix from file %s\n", input_file);
      exit(EXIT_FAILURE);
    }
  } else {
    printf("No input file specified!\n");
    exit(EXIT_FAILURE);
  } 

  if (do_verify){
    printf("Before LUD\n");
    print_matrix(matrix, matrix_dim);
    matrix_duplicate(matrix, &matrix_copy, matrix_dim);
  }

  int *ipiv = (int*)malloc(sizeof(int) * matrix_dim);
  if(cpuonly){
    stopwatch_start(&sw);
    lud_base(matrix, matrix_dim);
  } else if(lapack){
    stopwatch_start(&sw);
    lud_lapack(matrix, matrix_dim, ipiv);
  } else if(cula){
    printf("Initializing CULA\n");
    status = culaInitialize();
    checkStatus(status);
    stopwatch_start(&sw);
    status = lud_cula(matrix, matrix_dim, ipiv);
    checkStatus(status);
  } else {
    acc_init(acc_device_nvidia);
    stopwatch_start(&sw);
    lud_openacc(matrix, matrix_dim);
  }
  stopwatch_stop(&sw);
  printf("Time consumed(ms): %lf\n", 1000*get_interval_by_sec(&sw));

  if (do_verify){
    printf("After LUD\n");
    print_matrix(matrix, matrix_dim);
    printf(">>>Verify<<<<\n");
    lud_verify(matrix_copy, matrix, matrix_dim); 
    free(matrix_copy);
  }

  free(matrix);

  return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
