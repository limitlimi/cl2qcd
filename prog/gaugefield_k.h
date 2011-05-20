/** @file
 * Provides a class for gauge fields, including calculation of transport coefficients
 *
 */
#ifndef _GAUGEFIELDKH_
#define _GAUGEFIELDKH_

#include "gaugefield.h"
#include "opencl_k.h"

class Gaugefield_k : public Gaugefield {

  public:
    
	 /**
	 * Initialize gaugefield and devices for measurement of TK kappa
	 *
	 * @param[in] numdevs Number of wanted devices (so far, only 1 makes sense).
	 * @param[in] devicetypes Array of wanted cl_device_types for the devices.
	 * @param[in] input_parameters instance of inputparameters that contains information from input file
	 * @param[in,out] timer Return initialization time.
	 * @return Error code as defined in hmcerrs.h
	 */
	virtual hmc_error init(int numdevs, cl_device_type* devicetypes, inputparameters* input_parameters, usetimer* timer);
	
	hmc_float Q_plaquette();
	
  	/**
	 * Compute the transport coefficient kappa with the energy-momentum-tensor discretized by Karsch&Wyld
	 * @return Error code as defined in hmcerrs.h
	 */	
	hmc_error kappa_karsch ();
	/**
	 * Compute the transport coefficient kappa with the energy-momentum-tensor built by a Clover discretization
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error kappa_clover ();
	
	/**
	 * Returns the transport coefficient kappa computed by Karsch&Wyld's method
	 * @return Result for the transport coefficient kappa
	 */	
	hmc_float get_kappa_karsch ();
	
	/**
	 * Returns the transport coefficient kappa computed by Clover method
	 * @return Result for the transport coefficient kappa
	 */	
	hmc_float get_kappa_clover ();
	
	/**
	 * Set the transport coefficient kappa computed by Karsch&Wyld's method
	 * @param[in] in Result for the transport coefficient kappa
	 * @return Error code as defined in hmcerrs.h
	 */	
	hmc_error set_kappa_karsch (hmc_float in);
	
	/**
	 * Set the transport coefficient kappa computed by Clover method
	 * @param[in] in Result for the transport coefficient kappa
	 * @return Error code as defined in hmcerrs.h
	 */	
	hmc_error set_kappa_clover (hmc_float in);
	
	/**
	 * Compute the transport coefficient kappa with the energy-momentum-tensor discretized by Karsch&Wyld on GPU
	 * @param[in,out] timer time for measurement
	 * @return Error code as defined in hmcerrs.h
	 */	
	 hmc_error kappa_karsch_gpu (usetimer* timer_karsch);
	
	/**
	 * Compute the transport coefficient kappa with the energy-momentum-tensor discretized by Karsch&Wyld  on GPU
	 * @param[in] local_work_size OpenCL local_work_size
	 * @return Error code as defined in hmcerrs.h
	 */	
	hmc_error kappa_clover_gpu (usetimer* timer_clover);
  
	/**
	 * Returns private member opencl_k * devices
	 * @return devices of type opencl_k
	 */
 	Opencl * get_devices_k ();

  private:
// 	 Opencl_k * devices;
  hmc_float kappa_karsch_val;
  hmc_float kappa_clover_val;
  

};

#endif //_GAUGEFIELDKH_
