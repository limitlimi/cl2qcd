#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"

#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Gaugeforce
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const exec_name = "f_gauge";
std::string const version = "0.1";

class Device : public Opencl_Module_Hmc {
  cl_kernel testKernel;
  meta::Counter counter1, counter2, counter3, counter4;
public:
  Device(cl_command_queue queue, const meta::Inputparameters& params, int maxcomp, std::string double_ext, unsigned int dev_rank) : Opencl_Module_Hmc(params, &counter1, &counter2, &counter3, &counter4) {
    Opencl_Module_Hmc::init(queue, maxcomp, double_ext, dev_rank); /* init in body for proper this-pointer */
  };
  ~Device() {
    finalize();
  };
  void runTestKernel(cl_mem out, cl_mem gf, int gs, int ls);
  void fill_kernels();
  void clear_kernels();
};

const std::string SOURCEFILE = std::string(SOURCEDIR)
#ifdef _USEDOUBLEPREC_
  + "/tests/f_gauge_input_1";
#else
 + "/tests/f_gauge_input_1_single";
#endif

const std::string SOURCEFILE_REC12 = std::string(SOURCEDIR)
#ifdef _USEDOUBLEPREC_
  + "/tests/f_gauge_input_rec12";
#else
 + "/tests/f_gauge_input_rec12_single";
#endif

const char * PARAMS[] = {"foo", SOURCEFILE.c_str()};
const meta::Inputparameters INPUT(2, PARAMS);

const char * PARAMS_REC12[] = {"foo", SOURCEFILE_REC12.c_str()};
const meta::Inputparameters INPUT_REC12(2, PARAMS_REC12);

class Dummyfield : public Gaugefield_hybrid {
public:
  Dummyfield(cl_device_type device_type, meta::Inputparameters inputfile) : Gaugefield_hybrid(inputfile) {
    init(1, device_type);
    meta::print_info_hmc(exec_name.c_str(), inputfile);
  };
  virtual void init_tasks();
  virtual void finalize_opencl();
  
  hmc_float get_squarenorm();
  void runTestKernel();
  
private:
  void fill_buffers();
  void clear_buffers();
  cl_mem out;
  cl_mem sqnorm;
  hmc_float * sf_out;
};

void Dummyfield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new Device(queue[0], get_parameters(), get_max_compute_units(0), get_double_ext(0), 0);

	fill_buffers();
}

void Dummyfield::finalize_opencl()
{
	clear_buffers();
	Gaugefield_hybrid::finalize_opencl();
}

void fill_with_zero(hmc_float * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i] = 0.;
	}
	return;
}

void Dummyfield::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	cl_int err;
	cl_context context = opencl_modules[0]->get_context();

	size_t NUM_ELEMENTS_AE = meta::get_vol4d(get_parameters()) * NDIM * meta::get_su3algebrasize();
	sf_out = new hmc_float[NUM_ELEMENTS_AE];
	fill_with_zero(sf_out, NUM_ELEMENTS_AE);

	Device * device = static_cast<Device*>(opencl_modules[0]);
	out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, device->get_gaugemomentum_buffer_size(), 0, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	device->importGaugemomentumBuffer(out, reinterpret_cast<ae*>(sf_out));

	//create buffer for squarenorm on device
	sqnorm = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
}

void Device::fill_kernels()
{
	//one only needs some kernels up to now. to save time during compiling they are put in here by hand
	Opencl_Module_Hmc::fill_kernels();
	testKernel = createKernel("gauge_force") << basic_fermion_code << "types_hmc.h"  << "operations_gaugemomentum.cl" << "force_gauge.cl";
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	clReleaseMemObject(out);
	clReleaseMemObject(sqnorm);
	delete[] sf_out;
}

void Device::clear_kernels()
{
	clReleaseKernel(testKernel);
	Opencl_Module::clear_kernels();
}

void Device::runTestKernel(cl_mem out, cl_mem gf, int gs, int ls)
{
	cl_int err;
	err = clSetKernelArg(testKernel, 0, sizeof(cl_mem), &gf);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 1, sizeof(cl_mem), &out);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	enqueueKernel(testKernel, gs, ls);
}

hmc_float Dummyfield::get_squarenorm()
{
  //this always call the gaugemomentum sqnorm kernel
	static_cast<Device*>(opencl_modules[0])->set_float_to_gaugemomentum_squarenorm_device(out, sqnorm);
	// get stuff from device
	hmc_float result;
	cl_int err = clEnqueueReadBuffer(*queue, sqnorm, CL_TRUE, 0, sizeof(hmc_float), &result, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	return result;
}

void Dummyfield::runTestKernel()
{
	int gs = 0, ls = 0;
	if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = meta::get_spinorfieldsize(get_parameters());
		ls = 64;
	} else if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_CPU) {
		gs = opencl_modules[0]->get_max_compute_units();
		ls = 1;
	}
	Device * device = static_cast<Device*>(opencl_modules[0]);
	device->runTestKernel(out, device->get_gaugefield(), gs, ls);
}

BOOST_AUTO_TEST_CASE( F_GAUGE )
{
  hmc_float ref_val;
  logger.info() << "Test CPU and GPU version of kernel";
  logger.info() << "\tf_gauge";
  logger.info() << "against reference value";

  logger.info() << "Choosing reference value";
  //CP: I will not check if cpu and gpu have different starting conditions, this should never be the case...
  if(cpu.get_parameters().get_startcondition() == meta::Inputparameters::cold_start) {
    logger.info() << "Use cold config..." ;
    ref_val = 0.;
  } else{
    logger.info() << "Use specific config..";
    logger.warn() << "The reference value has to be adjusted manually if this config is changed!";
    ref_val = 52723.299867438494;
  }
  logger.info() << "reference value:\t" << ref_val;

  hmc_float prec = 1e-8;  
  logger.info() << "acceptance precision: " << prec;

  logger.info() << "Init CPU device";
  Dummyfield cpu(CL_DEVICE_TYPE_CPU, INPUT);
  logger.info() << "gaugeobservables: ";
  cpu.print_gaugeobservables_from_task(0, 0);
  cpu.runTestKernel();
  logger.info() << "|f_gauge|^2:";
  hmc_float cpu_res;
  cpu_res = cpu.get_squarenorm();
  logger.info() << cpu_res;
  BOOST_MESSAGE("Tested CPU");
  
  logger.info() << "Init GPU device";
  Dummyfield dummy(CL_DEVICE_TYPE_GPU, INPUT);
  logger.info() << "gaugeobservables: ";
  dummy.print_gaugeobservables_from_task(0, 0);
  dummy.runTestKernel();
  logger.info() << "|f_gauge|^2:";
  hmc_float gpu_res;
  gpu_res = dummy.get_squarenorm();
  logger.info() << gpu_res;
  BOOST_MESSAGE("Tested GPU");
  
  logger.info() << "Compare CPU result to reference value";
  BOOST_REQUIRE_CLOSE(cpu_res, ref_val, prec);

  logger.info() << "Compare GPU result to reference value";
  BOOST_REQUIRE_CLOSE(gpu_res, ref_val, prec);
  
  logger.info() << "Compare CPU and GPU results";
  BOOST_REQUIRE_CLOSE(cpu_res, gpu_res, prec);
}

BOOST_AUTO_TEST_CASE( F_GAUGE_REC12 )
{
  hmc_float ref_val;
  logger.info() << "Test CPU and GPU version of kernel";
  logger.info() << "\tf_gauge";
  logger.info() << "against reference value";

  logger.info() << "Choosing reference value";
  //CP: I will not check if cpu and gpu have different starting conditions, this should never be the case...
  if(cpu.get_parameters().get_startcondition() == meta::Inputparameters::cold_start) {
    logger.info() << "Use cold config..." ;
    ref_val = 0.;
  } else{
    logger.info() << "Use specific config..";
    logger.warn() << "The reference value has to be adjusted manually if this config is changed!";
    ref_val = 52723.3;
  }
  logger.info() << "reference value:\t" << ref_val;

  hmc_float prec = 1e-8;  
  logger.info() << "acceptance precision: " << prec;

  logger.info() << "Init CPU device";
  Dummyfield cpu(CL_DEVICE_TYPE_CPU, INPUT_REC12);
  logger.info() << "gaugeobservables: ";
  cpu.print_gaugeobservables_from_task(0, 0);
  cpu.runTestKernel();
  logger.info() << "|f_gauge|^2:";
  hmc_float cpu_res;
  cpu_res = cpu.get_squarenorm();
  logger.info() << cpu_res;
  BOOST_MESSAGE("Tested CPU");
  
  logger.info() << "Init GPU device";
  Dummyfield dummy(CL_DEVICE_TYPE_GPU, INPUT_REC12);
  logger.info() << "gaugeobservables: ";
  dummy.print_gaugeobservables_from_task(0, 0);
  dummy.runTestKernel();
  logger.info() << "|f_gauge|^2:";
  hmc_float gpu_res;
  gpu_res = dummy.get_squarenorm();
  logger.info() << gpu_res;
  BOOST_MESSAGE("Tested GPU");
  
  logger.info() << "Compare CPU result to reference value";
  BOOST_REQUIRE_CLOSE(cpu_res, ref_val, prec);

  logger.info() << "Compare GPU result to reference value";
  BOOST_REQUIRE_CLOSE(gpu_res, ref_val, prec);
  
  logger.info() << "Compare CPU and GPU results";
  BOOST_REQUIRE_CLOSE(cpu_res, gpu_res, prec);
}

