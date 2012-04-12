#ifndef GLOBAL_H
#define GLOBAL_H

#define MAXNAMESIZE 1024 // max filename length
#define SEED 1
/* increase this to reduce probability of random error */
/* increasing it also ups running time of "speedy" part of the code */
/* SP = 1 seems to be fine */
#define SP 1 // number of repetitions of speedy must be >=1

/* higher ITER --> more likely to get correct # of centers */
/* higher ITER also scales the running time almost linearly */
#define ITER 3 // iterate ITER* k log k times; ITER >= 1

//#define PRINTINFO //comment this out to disable output
#define PROFILE // comment this out to disable instrumentation code
//#define ENABLE_THREADS  // comment this out to disable threads
//#define INSERT_WASTE //uncomment this to insert waste computation into dist function

#define CACHE_LINE 512 // cache line in byte

#ifdef __cplusplus
extern "C" {
#endif

  extern int nproc;
  extern char *is_center, *switch_membership;
  extern int *center_table;

// instrumentation code
#ifdef PROFILE
extern double time_local_search;
extern double time_speedy;
extern double time_select_feasible;
extern double time_gain;
extern double time_shuffle;
extern double time_gain_dist;
extern double time_gain_init;
#endif 


#ifdef __cplusplus
}
#endif

#endif
