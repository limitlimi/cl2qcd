#include "prng.hpp"

#include "../../logger.hpp"
#include "../device.hpp"

using namespace std;

static std::string collect_build_options(hardware::Device * device, const meta::Inputparameters& params);

static std::string collect_build_options(hardware::Device * device, const meta::Inputparameters& params)
{
	std::ostringstream options;
	if(params.get_use_same_rnd_numbers() ) options <<  " -D _SAME_RND_NUMBERS_ ";
#ifdef USE_PRNG_NR3
	options << "-D USE_PRNG_NR3";
#elif defined(USE_PRNG_RANLUX)
	options << "-D USE_PRNG_RANLUX -D RANLUXCL_MAXWORKITEMS=" << hardware::buffers::get_prng_buffer_size(device);
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
	return options.str();
}

#ifdef USE_PRNG_NR3
void hardware::code::PRNG::copy_rndstate_to_device(nr3_state_dev* rndarray) const
{
	prng_buffer.load(rndarray);
}

void hardware::code::PRNG::copy_rndstate_from_device(nr3_state_dev* rndarray) const
{
	prng_buffer.dump(rndarray);
}
#endif // USE_PRNG_NR3

const hardware::buffers::PRNGBuffer& hardware::code::PRNG::get_prng_buffer() const noexcept
{
	return prng_buffer;
}

hardware::code::PRNG::PRNG(const meta::Inputparameters& params, hardware::Device * device)
	: Opencl_Module(params, device), prng_buffer(device)
{
#ifdef USE_PRNG_NR3
	// Prepare random number arrays, for each task and device separately
	const size_t num_rndstates = prng_buffer.get_elements();
	rndarray = new nr3_state_dev[num_rndstates];
	nr3_init_seeds(rndarray, "rand_seeds", num_rndstates);
	prng_buffer.load(rndarray);

	prng_code = ClSourcePackage(collect_build_options(get_device(), get_parameters())) << "random.cl";

	// Prepare random number arrays, for each task and device separately
	const size_t num_rndstates = prng_buffer.get_elements();
	rndarray = new nr3_state_dev[num_rndstates];
	nr3_init_seeds(rndarray, "rand_seeds", num_rndstates);
	prng_buffer.load(rndarray);
#elif defined(USE_PRNG_RANLUX)
	prng_code = ClSourcePackage(collect_build_options(get_device(), get_parameters())) << "ranluxcl/ranluxcl.cl" << "random.cl";
	init_kernel = createKernel("prng_ranlux_init") << get_device()->get_gaugefield_code()->get_sources() << prng_code << "random_ranlux_init.cl";
	initialize(&prng_buffer, get_parameters().get_host_seed() + 1);
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
}

hardware::code::PRNG::~PRNG()
{
#ifdef USE_PRNG_NR3
	delete [] rndarray;
#elif defined(USE_PRNG_RANLUX)
	cl_int clerr = clReleaseKernel(init_kernel);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
}

ClSourcePackage hardware::code::PRNG::get_sources() const noexcept
{
	return prng_code;
}

#ifdef USE_PRNG_RANLUX
void hardware::code::PRNG::initialize(const hardware::buffers::PRNGBuffer * buffer, cl_uint seed)
{
	cl_int clerr;
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes(init_kernel, &ls, &gs, &num_groups);
	if(seed > (10e9 / gs)) { // see ranluxcl source as to why
		/// @todo upgrade to newer ranluxcl to avoid this restcition
		throw Invalid_Parameters("Host seed is too large!", "<< 10e9", get_parameters().get_host_seed());
	}
	clerr = clSetKernelArg(init_kernel, 0, sizeof(cl_uint), &seed);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(init_kernel, 1, sizeof(cl_mem), buffer->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	get_device()->enqueue_kernel(init_kernel, gs, ls);
}
#endif /* USE_PRNG_RANLUX */