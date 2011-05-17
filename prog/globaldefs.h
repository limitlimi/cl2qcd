/** @file
 * Global parameters, specified as macros
 *
 * @todo Some, or even most of these, would be better represented as constants
 */

#ifndef _GLOBALSH_
#define _GLOBALSH_

/** Number of colors */
#define NC 3
#define NSPIN 4

/** Number of dimensions of the lattice */
#define NDIM 4

#ifndef _INKERNEL_

/** Number of lattice sites in time direction (t) */
#define NTIME 8
/** Number of lattice sites in space direction (x,y,z each) */
#define NSPACE 8

/** Spatial volume of the lattice */
#define VOLSPACE NSPACE*NSPACE*NSPACE
/** 4-Dimensional Volume of the lattice */
#define VOL4D VOLSPACE*NTIME

#define SPINORSIZE NSPIN*NC
#define HALFSPINORSIZE NSPIN/2*NC
#define SPINORFIELDSIZE SPINORSIZE*NTIME*VOLSPACE
#define EOPREC_SPINORFIELDSIZE SPINORSIZE*NTIME*VOLSPACE/2
#define GAUGEMOMENTASIZE NDIM*VOL4D*(NC*NC-1)
#define GAUGEFIELDSIZE NC*NC*NDIM*VOL4D
#define SU3ALGEBRASIZE NC*NC-1

//startconditions:
#define START_FROM_SOURCE 2
#define COLD_START 0
#define HOT_START 1

#endif //_INKERNEL_

//EVEN ODD
#define EVEN 0
#define ODD 1

#define TRUE 1
#define FALSE 0

/**
 * PI
 * @todo Rather use PI from stdlib
 */
#define PI 	3.14159265358979

#define su2_entries 4

/** Number of threads to use for OpenCL kernels */
#ifdef _USEGPU_
#define NUMTHREADS 128
#else
#define NUMTHREADS 1
#endif

// Definition of numeric constants for the symmetric structure constants d_ijk of su(3)
/** 1/2 */
#define F_1_2   (static_cast<hmc_float>(0.5))
/** 1/(2*sqrt(3)) */
#define F_1_2S3 (static_cast<hmc_float>(0.288675134594813))
/** 1/sqrt(3) */
#define F_1_S3  (static_cast<hmc_float>(0.577350269189626))

// SL: not sure if those should be here: they define details for the exact exponentiation of su(3) into SU(3)
#define _EXACT_EXPONENTIATION_MAX_POWER_  (50)
#define _EXACT_EXPONENTIATION_ACCURACY_   (2.0E-16)

#endif


