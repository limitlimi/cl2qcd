/** @file
 * Provides class Gaugefield_hmc, inherited from Gaugefield_inversion, adds HMC capabilities
 *
 * @todo global- and local work sizes should be moved inside Opencl_fermions (LZ)
 */
#ifndef _GAUGEFIELDHMCH_
#define _GAUGEFIELDHMCH_

#include "gaugefield.h"
#include "gaugefield_inversion.h"
#include "opencl_fermions.h"
#include "opencl_hmc.h"
#include "types_hmc.h"

class Gaugefield_hmc : public Gaugefield_inversion {
  public:

  /**
   * Initializes the devices, to be called by init()
   * @return Error code as defined in hmcerrs.h
   * @param devicetypes array of cl_device_type handles
   */
  virtual hmc_error init_devices(cl_device_type* devicetypes);

  /**
   * Free gaugefield and device allocations.
   */
  virtual hmc_error finalize();

  /**
   * Free device, called by finalize
   */
  virtual hmc_error free_devices();
  
  /**
   * Returns private member opencl_k * devices
   * @return devices of type opencl_k
   * @todo LZ: CHECK IF THIS MAKES SENSE AT ALL!!!!
   */
  Opencl_hmc * get_devices_hmc ();

	hmc_error perform_hmc_step(int dev, inputparameters *parameters, hmc_observables *obs, int iter, hmc_float rnd_number, usetimer* copy_to_from_dev_timer, usetimer* copy_on_dev_timer);
	
	void print_hmcobservables(hmc_observables obs, int iter, std::string filename);
 private:

};



#endif //_GAUGEFIELDHMCH_
