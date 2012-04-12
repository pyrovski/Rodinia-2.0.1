#ifndef KERNEL_H
#define KERNEL_H

#include <sys/time.h>
#include <stdlib.h>
#include <string.h>

#ifdef OMP
#include <omp.h>
#endif

#include "global.h"

#ifdef __cplusplus
extern "C" {
#endif

/* this structure represents a point */
/* these will be passed around to avoid copying coordinates */
typedef struct {
  float weight;
  float *coord;
  long assign;  /* number of point where this one is assigned */
  float cost;  /* cost of that assignment, weight*distance */
} Point;

/* this is the array of points */
typedef struct {
  long num; /* number of points; may not be N if this is a sample */
  int dim;  /* dimensionality */
  Point *p; /* the array itself */
} Points;

inline double gettime() {
  struct timeval t;
  gettimeofday(&t,NULL);
  return (double)t.tv_sec+t.tv_usec*1e-6;
}

inline int isIdentical(float *i, float *j, int D)
// tells whether two points of D dimensions are identical
{
  int a = 0;
  int equal = 1;

  while (equal && a < D) {
    if (i[a] != j[a]) equal = 0;
    else a++;
  }
  if (equal) return 1;
  else return 0;

}

/* comparator for floating point numbers */
inline int floatcomp(const float *i, const float *j)
{
  float a, b;
  a = *(float *)(i);
  b = *(float *)(j);
  if (a > b) return (1);
  if (a < b) return (-1);
  return(0);
}

/* compute Euclidean distance squared between two points */
inline float dist(Point p1, Point p2, int dim)
{
  int i;
  float result=0.0;
  for (i=0;i<dim;i++)
    result += (p1.coord[i] - p2.coord[i])*(p1.coord[i] - p2.coord[i]);
#ifdef INSERT_WASTE
  double s = waste(result);
  result += s;
  result -= s;
#endif
  return(result);
}

double 
pgain(long x, Points *points, double z, long int *numcenters, int pid, 
      pthread_barrier_t* barrier);

#ifdef __cplusplus
}
#endif

#endif
