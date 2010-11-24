#ifndef _GLOBALSH_
#define _GLOBALSH_

#define NC 3
#define NSPIN 4
#define NDIM 4

//it could be a good idea to define Nt and Ns at compile time
//usually you stick to one volume for quite a while anyways...
#define NSPACE 32
#define NTIME 12
#define VOLSPACE NSPACE*NSPACE*NSPACE
#define VOL4D VOLSPACE*NTIME

//Seed for Random
const unsigned long long int seed = 500000;

const int su2_entries = 4;

#define PI 	3.14159265358979

#endif

#ifdef _OPENMP
# include <omp.h>
#endif
