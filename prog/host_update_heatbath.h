/** @file
 * Heatbath update
 */

#ifndef _UPDATE_HEATBATHH_
#define _UPDATE_HEATBATHH_

#include "globaldefs.h"
#include "types.h"
#include "hmcerrs.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "host_operations_spinor.h"
#include "host_operations_spinorfield.h"
#include "host_geometry.h"
#include "host_random.h"

void calc_staple(hmc_gaugefield * field, hmc_staplematrix * dest, const int pos, const int t, const int mu_in);

void heatbath_update (hmc_gaugefield * gaugefield, const hmc_float beta);
void heatbath_overrelax (hmc_gaugefield * gaugefield);
void heatbath_update_checkerboard (hmc_gaugefield * gaugefield, const hmc_float beta);
void heatbath_overrelax_checkerboard (hmc_gaugefield * gaugefield);

#endif /* _UPDATE_HEATBATHH_ */
