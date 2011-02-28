#include "hmc.h"

int main(int argc, char* argv[]) {

  char* progname = argv[0];
  print_hello(progname);

  char* inputfile = argv[1];
  inputparameters parameters;
  parameters.readfile(inputfile);
  print_info(&parameters);

  stringstream gaugeout_name;
  gaugeout_name<<"gaugeobservables_beta"<<parameters.get_beta();
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialization
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  sourcefileparameters parameters_source;
  hmc_gaugefield * gaugefield;
  gaugefield = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));
  hmc_rndarray rndarray;

  init_gaugefield(gaugefield,&parameters,&inittime);
  init_random_seeds(rnd, rndarray, &inittime);
   
  opencl gpu(CL_DEVICE_TYPE_CPU, &inittime);

  cout << "initial values of observables:\n\t" ;
  print_gaugeobservables(gaugefield, &polytime, &plaqtime);

  gpu.copy_gaugefield_to_device(gaugefield, &copytime);
  gpu.copy_rndarray_to_device(rndarray, &copytime);

  gpu.testing(gaugefield);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Heatbath
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int nsteps = parameters.get_heatbathsteps();
  cout<<"perform "<<nsteps<<" heatbath steps on OpenCL device..."<<endl;
  for(int i = 0; i<nsteps; i++){
    gpu.run_heatbath(parameters.get_beta(), local_work_size, global_work_size, &updatetime);
    gpu.run_overrelax(parameters.get_beta(), local_work_size, global_work_size, &overrelaxtime);
    if( ( (i+1)%parameters.get_writefrequency() ) == 0 ) {
       gpu.gaugeobservables(local_work_size, global_work_size, &plaq, &tplaq, &splaq, &pol, &plaqtime, &polytime);
       print_gaugeobservables(plaq, tplaq, splaq, pol, i, gaugeout_name.str());
    }
    if( parameters.get_saveconfigs()==TRUE && ( (i+1)%parameters.get_savefrequency() ) == 0 ) {
      gpu.get_gaugefield_from_device(gaugefield, &copytime);
      save_gaugefield(gaugefield, &parameters, i);
      print_gaugeobservables(gaugefield, &plaqtime, &polytime, i, gaugeout_name.str());
    }
  }

  gpu.get_gaugefield_from_device(gaugefield, &copytime);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Final Output
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  totaltime.add();
  save_gaugefield(gaugefield, &parameters, nsteps);  
  time_output(&totaltime, &inittime, &polytime, &plaqtime, &updatetime, &overrelaxtime, &copytime);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // free variables
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  free(gaugefield);
  gpu.finalize();
  
  return HMC_SUCCESS;
}
