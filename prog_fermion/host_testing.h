#ifndef _TESTINGH_
#define _TESTINGH_

#include <cstdlib>
#include <cstdio>

#include "hmcerrs.h"
#include "globaldefs.h"
#include "types.h"
#include "host_operations.h"
#include "host_geometry.h"
#include "host_fermionsolver.h"
#include "host_simplefermionsolver.h"
#include "host_input.h"

void testing_correlator(hmc_gaugefield* gf, inputparameters* parameters);
void testing_spinor();
void testing_fermionmatrix();
void testing_eoprec_spinor();
void testing_fermsolve(hmc_gaugefield* gf);
void print_su3mat(hmc_su3matrix* A);
void print_staplemat(hmc_staplematrix* A);
void print_fullspinorfield(hmc_spinor* in);
void testing_su3mat();
void testing_gaugefield();
void testing_geometry();
void testing_su3matrix(hmc_gaugefield * in, int spacepos, int timepos);
void testing_adjoin(hmc_gaugefield * in, int spacepos, int timepos);
void testing_det_global(hmc_gaugefield * in);
void testing_matrix_ops(hmc_gaugefield * in);

#endif
