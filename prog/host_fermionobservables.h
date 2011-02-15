#ifndef _FERMIONOBSERVABLESH_
#define _FERMIONOBSERVABLESH_

#include <cstdlib>
#include <cstring>
#include <iostream>

#include "host_operations.h"
#include "host_geometry.h"
#include "globaldefs.h"
#include "hmcerrs.h"
#include "types.h"
#include "host_solver.h"

hmc_error simple_correlator(hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, int cgmax);

#endif