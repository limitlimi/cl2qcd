#include "opencl_module_hmc.h"

#include "logger.hpp"
#include "meta/util.hpp"
#include "hardware/device.hpp"

using namespace std;

static std::string collect_build_options(hardware::Device * device, const meta::Inputparameters& params);

static std::string collect_build_options(hardware::Device * device, const meta::Inputparameters& params)
{
	using namespace hardware::buffers;

	std::ostringstream options;
	options <<  "-D BETA=" << params.get_beta() << " -D GAUGEMOMENTASIZE=" << meta::get_vol4d(params) * NDIM;
	//in case of tlsym gauge action
	if(meta::get_use_rectangles(params) == true) {
		options <<  " -D C0=" << meta::get_c0(params) << " -D C1=" << meta::get_c1(params);
	}
	if(check_Gaugemomentum_for_SOA(device)) {
		options << " -D GAUGEMOMENTA_STRIDE=" << get_Gaugemomentum_buffer_stride(meta::get_vol4d(params) * NDIM, device);
	}
	return options.str();
}

void Opencl_Module_Hmc::fill_kernels()
{
	basic_hmc_code = get_device()->get_fermion_code()->get_sources() << ClSourcePackage(collect_build_options(get_device(), get_parameters())) << "types_hmc.h" << "operations_gaugemomentum.cl";
	ClSourcePackage prng_code = get_device()->get_prng_code()->get_sources();

	//init kernels for HMC
	if(get_parameters().get_use_eo() == true) {
		fermion_force_eo = createKernel("fermion_force_eo") << basic_hmc_code << "fermionmatrix.cl" << "force_fermion_eo.cl";
	} 
	fermion_force = createKernel("fermion_force") << basic_hmc_code << "fermionmatrix.cl" << "force_fermion.cl";
	_set_zero_gaugemomentum = createKernel("set_zero_gaugemomentum") << basic_hmc_code <<  "gaugemomentum_zero.cl";
	md_update_gaugefield = createKernel("md_update_gaugefield") << basic_hmc_code << "md_update_gaugefield.cl";
	md_update_gaugemomenta = createKernel("md_update_gaugemomenta") << basic_hmc_code  << "md_update_gaugemomenta.cl";
	gauge_force = createKernel("gauge_force") << basic_hmc_code  << "force_gauge.cl";
	if(meta::get_use_rectangles(get_parameters()) == true) {
		//at the time of writing this kernel, the OpenCL compiler crashed the kernel using optimizations
		gauge_force_tlsym = createKernel("gauge_force_tlsym") << basic_hmc_code << "force_gauge_tlsym.cl";
	}
	if(get_parameters().get_use_smearing() == true) {
		stout_smear_fermion_force = createKernel("stout_smear_fermion_force") << basic_hmc_code << "force_fermion_stout_smear.cl";
	}
	gaugemomentum_squarenorm = createKernel("gaugemomentum_squarenorm") << basic_hmc_code << "gaugemomentum_squarenorm.cl";
	if(get_device()->get_prefers_soa()) {
		gaugemomentum_convert_to_soa = createKernel("gaugemomentum_convert_to_soa") << basic_hmc_code << "gaugemomentum_convert.cl";
		gaugemomentum_convert_from_soa = createKernel("gaugemomentum_convert_from_soa") << basic_hmc_code << "gaugemomentum_convert.cl";
	} else {
		gaugemomentum_convert_to_soa = 0;
		gaugemomentum_convert_from_soa = 0;
	}
}

void Opencl_Module_Hmc::clear_kernels()
{
	cl_uint clerr = CL_SUCCESS;

	logger.debug() << "release HMC-kernels.." ;
	if(get_parameters().get_use_eo() == true) {
		clerr = clReleaseKernel(fermion_force_eo);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	} else {
		clerr = clReleaseKernel(fermion_force);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	clerr = clReleaseKernel(md_update_gaugefield);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(md_update_gaugemomenta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(gauge_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(meta::get_use_rectangles(get_parameters()) == true) {
		clerr = clReleaseKernel(gauge_force_tlsym);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	clerr = clReleaseKernel(_set_zero_gaugemomentum);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(get_parameters().get_use_smearing() == true) {
		clerr = clReleaseKernel(stout_smear_fermion_force);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
}

void Opencl_Module_Hmc::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);
	return;
}

////////////////////////////////////////////////////
//Access to members

const hardware::buffers::Gaugemomentum * Opencl_Module_Hmc::get_clmem_p()
{
	return &clmem_p;
}

const hardware::buffers::Gaugemomentum * Opencl_Module_Hmc::get_clmem_new_p()
{
	return &clmem_new_p;
}

const hardware::buffers::SU3 * Opencl_Module_Hmc::get_new_u()
{
	return &new_u;
}

const hardware::buffers::Plain<spinor> * Opencl_Module_Hmc::get_clmem_phi()
{
	return &clmem_phi;
}

const hardware::buffers::Plain<spinor> * Opencl_Module_Hmc::get_clmem_phi_inv()
{
	return &clmem_phi_inv;
}

const hardware::buffers::Spinor * Opencl_Module_Hmc::get_clmem_phi_eo()
{
	return &clmem_phi_eo;
}

const hardware::buffers::Spinor * Opencl_Module_Hmc::get_clmem_phi_inv_eo()
{
	return &clmem_phi_inv_eo;
}

const hardware::buffers::Plain<spinor> * Opencl_Module_Hmc::get_clmem_phi_mp()
{
	return &clmem_phi_mp;
}

const hardware::buffers::Spinor * Opencl_Module_Hmc::get_clmem_phi_mp_eo()
{
	return &clmem_phi_mp_eo;
}

hardware::buffers::Plain<hmc_float> * Opencl_Module_Hmc::get_clmem_s_fermion_init()
{
	return &clmem_s_fermion_init;
}

hardware::buffers::Plain<hmc_float> * Opencl_Module_Hmc::get_clmem_s_fermion_mp_init()
{
	return &clmem_s_fermion_mp_init;
}

size_t Opencl_Module_Hmc::get_read_write_size(const std::string& in) const
{
//Depending on the compile-options, one has different sizes...
	size_t D = meta::get_float_size(get_parameters());
	//this returns the number of entries in an su3-matrix
	size_t R = meta::get_mat_size(get_parameters());
	//this is the number of spinors in the system (or number of sites)
	size_t S = meta::get_spinorfieldsize(get_parameters());
	size_t Seo = meta::get_eoprec_spinorfieldsize(get_parameters());
	//this is the number of links in the system (and of gaugemomenta)
	size_t G = meta::get_vol4d(get_parameters()) * NDIM;
	//factor for complex numbers
	int C = 2;
	//this is the same as in the function above
	//NOTE: 1 spinor has NC*NDIM = 12 complex entries
	//NOTE: 1 ae has NC*NC-1 = 8 real entries
	int A = meta::get_su3algebrasize();
	if (in == "md_update_gaugefield") {
		//this kernel reads 1 ae and 1 su3 matrix and writes 1 su3 matrix for every link
		return (A + (1 + 1) * R * C) * D * G;
	}
	if (in == "md_update_gaugemomenta") {
		//this kernel reads 2 ae and writes 1 ae per link
		return ((2 + 1) * A) * D * G;
	}
	if (in == "gauge_force") {
		//this kernel reads ingredients for 1 staple plus 1 su3matrix and writes 1 ae for every link
		return G * D * (R * C * ( 6 * (NDIM - 1) + 1 ) + A );
	}
	if (in == "gauge_force_tlsym") {
		//this kernel reads ingredients for 1 rect-staple plus 1 su3matrix and writes 1 ae for every link
		//the rect staple is the same as the normal staple, but with 3 add. contributions and 2 add. matrices in each contribution, so instead of 2 * 3 = 6, one reads 6 * 5 matrices in each direction
		return G * D * (R * C * ( 6 * 5 * (NDIM - 1 ) + 1 ) + A );
	}
	if (in == "fermion_force") {
		//this kernel reads 16 spinors, 8 su3matrices and writes 8 ae per site
		//NOTE: the kernel now runs over all ae instead of all sites, but this must be equivalent!
		return (C * 12 * (16) + C * 8 * R + 8 * A) * D * S;
	}
	if (in == "fermion_force_eo") {
		//this kernel reads 16 spinors, 8 su3matrices and writes 1 ae per site
		return (C * 12 * (16) + C * 8 * R + 8 * A) * D * Seo;
	}
	if (in == "set_zero_gaugemomentum;") {
		//this kernel writes 1 ae per link
		return G * D * A;
	}
	if (in == "gaugemomentum_squarenorm") {
		//this kernel reads 1 ae and writes 1 float per link
		return G * D * ( A + 1);
	}
	if (in == "stout_smear_fermion_force") {
		return 10000000000000000000;
	}
	return 0;
}

uint64_t Opencl_Module_Hmc::get_flop_size(const std::string& in) const
{
	//this is the number of spinors in the system (or number of sites)
	size_t S = meta::get_spinorfieldsize(get_parameters());
	size_t Seo = meta::get_eoprec_spinorfieldsize(get_parameters());
	//this is the number of links in the system (and of gaugemomenta)
	uint64_t G = meta::get_vol4d(get_parameters()) * NDIM;
	//NOTE: 1 ae has NC*NC-1 = 8 real entries
	uint64_t A = meta::get_su3algebrasize();
	//this returns the number of entries in an su3-matrix
	uint64_t R = meta::get_mat_size(get_parameters());
	//this is the same as in the function above
	if (in == "md_update_gaugefield") {
		//this kernel performs one exp(i ae) ( = 327 flops + 1 su3 mult ) and 1 su3 mult per link
		return (meta::get_flop_su3_su3() * ( 1 + 1)  + 327 ) * G;
	}
	if (in == "md_update_gaugemomenta") {
		//this kernel performs 1 real mult and 1 real add per ae
		return (1 + 1) * A * G;
	}
	if (in == "gauge_force") {
		//this kernel calculates 1 staple (= 4*ND-1 su3_su3 + 2*ND-1 su3_add), 1 su3*su3, 1 tr_lambda_u (19 flops) plus 8 add and 8 mult per ae
		return ( 4 * (NDIM - 1) * meta::get_flop_su3_su3() + 2 * (NDIM - 1) * 18 + 1 * meta::get_flop_su3_su3() + 19  + A * ( 1 + 1 )
		       ) * G;
	}
	if (in == "gauge_force_tlsym") {
		//this kernel calculates 1 rect-staple (= 24*ND-1 su3_su3 + 6*ND-1 su3_add), 1 su3*su3, 1 tr_lambda_u (19 flops) plus 8 add and 8 mult per ae
		//24 = 6 contr. per dir, 4 mat_mat per contr.
		return ( 24 * (NDIM - 1) * meta::get_flop_su3_su3() + 6 * (NDIM - 1) * 18 + 1 * meta::get_flop_su3_su3() + 19  + A * ( 1 + 1 )
		       ) * G;
	}
	if (in == "fermion_force") {
		//this kernel performs NDIM * ( 4 * su3vec_acc (6 flops) + tr(v*u) (126 flops) + tr_lambda_u(19 flops) + update_ae(8*2 flops) + su3*su3 + su3*complex (flop_complex_mult * R ) ) per site
		//NOTE: the kernel now runs over all ae instead of all sites, but this must be equivalent!
		return Seo * NDIM * ( 4 * 6 + 126 + 19 + 8 * 2 + meta::get_flop_su3_su3() + meta::get_flop_complex_mult() * R );
	}
	if (in == "fermion_force_eo") {
		//this kernel performs NDIM * ( 4 * su3vec_acc (6 flops) + tr(v*u) (126 flops) + tr_lambda_u(19 flops) + update_ae(8*2 flops) + su3*su3 + su3*complex (flop_complex_mult * R ) ) per site
		return Seo * NDIM * ( 4 * 6 + 126 + 19 + 8 * 2 + meta::get_flop_su3_su3() + meta::get_flop_complex_mult() * R );
	}
	if (in == "set_zero_gaugemomentum;") {
		//this kernel performs 0 mults
		return 0;
	}
	if (in == "gaugemomentum_squarenorm") {
		//this kernel performs 8 real mults and 8-1 real adds per ae
		return (8 + 7) * A * G;
	}
	if (in == "stout_smear_fermion_force") {
		return 10000000000000000000;
	}
	return 0;
}

void Opencl_Module_Hmc::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	Opencl_Module::print_profiling(filename, md_update_gaugefield);
	Opencl_Module::print_profiling(filename, md_update_gaugemomenta);
	Opencl_Module::print_profiling(filename, gauge_force);
	Opencl_Module::print_profiling(filename, gauge_force_tlsym);
	Opencl_Module::print_profiling(filename, fermion_force);
	Opencl_Module::print_profiling(filename, fermion_force_eo);
	Opencl_Module::print_profiling(filename, _set_zero_gaugemomentum);
	Opencl_Module::print_profiling(filename, gaugemomentum_squarenorm);
	Opencl_Module::print_profiling(filename, stout_smear_fermion_force);
}

////////////////////////////////////////////////////
//Methods needed for the HMC-algorithm
void Opencl_Module_Hmc::generate_spinorfield_gaussian(const hardware::buffers::PRNGBuffer * prng)
{
  auto spinor_code = get_device()->get_spinor_code();
	if(get_parameters().get_use_eo() == true) {
	  spinor_code->generate_gaussian_spinorfield_eo_device(get_clmem_phi_inv_eo(), prng);
	} else {
	  spinor_code->generate_gaussian_spinorfield_device(get_clmem_phi_inv(), prng);
	}
	return;
}

void Opencl_Module_Hmc::md_update_spinorfield(const hardware::buffers::SU3 * gaugefield, hmc_float kappa, hmc_float mubar)
{
	auto fermion_code = get_device()->get_fermion_code();

	//suppose the initial gaussian field is saved in clmem_phi_inv (see above).
	//  then the "phi" = Dpsi from the algorithm is stored in clmem_phi
	//  which then has to be the source of the inversion
	if(get_parameters().get_use_eo() == true) {
		fermion_code->Qplus_eo (&clmem_phi_inv_eo, &clmem_phi_eo , gaugefield, kappa, mubar);
		if(logger.beDebug()) fermion_code->print_info_inv_field(&clmem_phi_eo, true, "\tinit field after update ");
	} else {
		fermion_code->Qplus(&clmem_phi_inv, &clmem_phi , gaugefield, kappa, mubar);
		if(logger.beDebug()) fermion_code->print_info_inv_field(&clmem_phi, false, "\tinit field after update ");
	}
}

void Opencl_Module_Hmc::md_update_spinorfield_mp(usetimer * solvertimer, const hardware::buffers::SU3 * gaugefield)
{
	using namespace hardware::buffers;

	auto fermion_code = get_device()->get_fermion_code();
	auto spinor_code = get_device()->get_spinor_code();

	///@todo solvertimer is not used here yet...
	//suppose the initial gaussian field is saved in clmem_phi_inv (see above).
	//  then the "phi" = Dpsi from the algorithm is stored in clmem_phi
	//  which then has to be the source of the inversion
	//in the mass preconditioning case, this is a bit more complicated and involves an inversion
	if(get_parameters().get_use_eo() == true) {
		//CP: Init tmp spinorfield
		Spinor sf_eo_tmp(clmem_phi_inv_eo.get_elements(), get_device());

		//sf_eo_tmp = Qplus_eo(light_mass) phi_inv_eo
		fermion_code->Qplus_eo (&clmem_phi_inv_eo, &sf_eo_tmp , gaugefield);

		//Now one needs ( Qplus_eo )^-1 (heavy_mass) using sf_eo_tmp as source to get phi_mp_eo
		//use always bicgstab here
		logger.debug() << "\t\t\tstart solver";

		/** @todo at the moment, we can only put in a cold spinorfield
		 * or a point-source spinorfield as trial-solution
		 */
		spinor_code->set_zero_spinorfield_eoprec_device(get_clmem_phi_mp_eo());

		fermion_code->gamma5_eo_device(get_clmem_phi_mp_eo());

		int converged = -1;
		if(logger.beDebug()) fermion_code->print_info_inv_field(get_clmem_phi_mp_eo(), true, "\tinv. field before inversion ");
		if(logger.beDebug()) fermion_code->print_info_inv_field(&sf_eo_tmp, true, "\tsource before inversion ");
		converged = fermion_code->bicgstab_eo(::Qplus_eo(fermion_code), this->get_clmem_phi_mp_eo(), &sf_eo_tmp, &new_u, get_parameters().get_solver_prec(), get_parameters().get_kappa_mp(), meta::get_mubar_mp(get_parameters()));
		if (converged < 0) {
			if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
			else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
		} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
		if(logger.beDebug()) fermion_code->print_info_inv_field(get_clmem_phi_mp_eo(), true, "\tinv. field after inversion ");
		if(logger.beDebug()) fermion_code->print_info_inv_field(&clmem_phi_mp_eo, true, "\tinit field after update ");
	} else {
		//CP: Init tmp spinorfield
		Plain<spinor> sf_tmp(clmem_phi_inv.get_elements(), get_device());

		//sf_tmp = Qplus(light_mass) phi_inv
		fermion_code->Qplus (&clmem_phi_inv, &sf_tmp , gaugefield);

		//Now one needs ( Qplus )^-1 (heavy_mass) using sf_tmp as source to get phi_mp
		//use always bicgstab here
		logger.debug() << "\t\t\tstart solver";

		/** @todo at the moment, we can only put in a cold spinorfield
		 * or a point-source spinorfield as trial-solution
		 */
		spinor_code->set_zero_spinorfield_device(get_clmem_phi_mp());
		fermion_code->gamma5_device(get_clmem_phi_mp());

		int converged = -1;
		if(logger.beDebug()) fermion_code->print_info_inv_field(get_clmem_phi_mp(), false, "\tinv. field before inversion ");
		if(logger.beDebug()) fermion_code->print_info_inv_field(&sf_tmp, false, "\tsource before inversion ");
		converged = fermion_code->bicgstab(::Qplus(fermion_code), this->get_clmem_phi_mp(), &sf_tmp, &new_u, get_parameters().get_solver_prec(), get_parameters().get_kappa_mp(), get_mubar_mp(get_parameters()));
		if (converged < 0) {
			if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
			else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
		} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
		if(logger.beDebug()) fermion_code->print_info_inv_field(get_clmem_phi_mp(), false, "\tinv. field after inversion ");
		if(logger.beDebug()) fermion_code->print_info_inv_field(&clmem_phi_mp, false, "\tinit field after update ");
	}
}

//this function takes to args kappa and mubar because one has to use it with different masses when mass-prec is used and when not
void Opencl_Module_Hmc::calc_fermion_force(usetimer * solvertimer, hmc_float kappa, hmc_float mubar)
{
	auto fermion_code = get_device()->get_fermion_code();
	auto spinor_code = get_device()->get_spinor_code();

	int converged = -1;
	if(get_parameters().get_use_eo() == true) {
		//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			/**
			 * The first inversion calculates
			 * X_even = phi = (Qplusminus_eo)^-1 psi
			 * out of
			 * Qplusminus_eo phi_even = psi
			 */
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			/**
			 * Trial solution for the spinorfield
			 */
			spinor_code->set_eoprec_spinorfield_cold_device(fermion_code->get_inout_eo());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\t\t\tinv. field before inversion ");
			converged = fermion_code->cg_eo(::QplusQminus_eo(fermion_code), fermion_code->get_inout_eo(), this->get_clmem_phi_eo(), &new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\t\t\tinv. field after inversion ");

			/**
			 * Y_even is now just
			 *  Y_even = (Qminus_eo) X_even = (Qminus_eo) (Qplusminus_eo)^-1 psi =
			 *    = (Qplus_eo)^-1 psi
			 */
			fermion_code->Qminus_eo(fermion_code->get_inout_eo(), &clmem_phi_inv_eo, &new_u, kappa, mubar);
		} else {
			///@todo if wanted, solvertimer has to be used here..
			//logger.debug() << "\t\tcalc fermion force ingredients using bicgstab with eo.";
			/**
			 * The first inversion calculates
			 * Y_even = phi = (Qplus_eo)^-1 psi
			 * out of
			 * Qplus_eo phi = psi
			 * This is also the energy of the final field!
			*/
			//logger.debug() << "\t\tcalc Y_even...";
			//logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			/**
			 * Trial solution for the spinorfield
			 */
			spinor_code->set_zero_spinorfield_eoprec_device(fermion_code->get_inout_eo());
			fermion_code->gamma5_eo_device(fermion_code->get_inout_eo());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\t\t\tinv. field before inversion ");
			if(logger.beDebug()) fermion_code->print_info_inv_field(get_clmem_phi_eo(), true, "\t\t\tsource before inversion ");
			converged = fermion_code->bicgstab_eo(::Qplus_eo(fermion_code), fermion_code->get_inout_eo(), this->get_clmem_phi_eo(), &new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\t\t\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			hardware::buffers::copyData(&clmem_phi_inv_eo, fermion_code->get_inout_eo());

			/**
			 * Now, one has to calculate
			 * X = (Qminus_eo)^-1 Y = (Qminus_eo)^-1 (Qplus_eo)^-1 psi = (QplusQminus_eo)^-1 psi ??
			 * out of
			 * Qminus_eo clmem_inout_eo = clmem_phi_inv_eo
			 * Therefore, when calculating the final energy of the spinorfield,
			 *  one can just take clmem_phi_inv_eo (see also above)!!
			 */

			//logger.debug() << "\t\tcalc X_even...";
			//copy former solution to clmem_source
			hardware::buffers::copyData(fermion_code->get_source_even(), fermion_code->get_inout_eo());

			//this sets clmem_inout cold as trial-solution
			spinor_code->set_eoprec_spinorfield_cold_device(fermion_code->get_inout_eo());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\t\t\tinv. field before inversion ");
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_source_even(), true, "\t\t\tsource before inversion ");
			converged = fermion_code->bicgstab_eo(::Qminus_eo(fermion_code), fermion_code->get_inout_eo(), fermion_code->get_source_even(), &new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\t\t\tinv. field after inversion ");
		}
		/**
		 * At this point, one has calculated X_odd and Y_odd.
		 * If one has a fermionmatrix
		 *  M = R + D
		 * these are:
		 *  X_odd = -R(-mu)_inv D X_even
		 *  Y_odd = -R(mu)_inv D Y_even
		 */

		/** @fixme below usages of dslash should work, but only because we use bicgstab above
		    proper implementation needs to make sure this is always the case */

		///@NOTE the following calculations could also go in a new function for convenience
		//calculate X_odd
		//therefore, clmem_tmp_eo_1 is used as intermediate state. The result is saved in clmem_inout, since
		//  this is used as a default in the force-function.
		if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
			fermion_code->dslash_eo_device(fermion_code->get_inout_eo(), fermion_code->get_tmp_eo_1(), &new_u, ODD, kappa);
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_1(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			fermion_code->dslash_eo_device(fermion_code->get_inout_eo(), fermion_code->get_tmp_eo_1(), &new_u, ODD, kappa);
			fermion_code->M_tm_inverse_sitediagonal_minus_device(fermion_code->get_tmp_eo_1(), fermion_code->get_tmp_eo_2(), mubar);
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_2(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		}

		/**
		 * If needed, this gives additional debug informations
		 */
// hard debugging ist not SOA safe in the shown vesion!
//		bool debug_hard = false;
//		cl_mem x_tmp, y_tmp;
//		cl_mem x_eo_tmp, x_eo_tmp2, y_eo_tmp, y_eo_tmp2;

// will break on SOA devices
//		if(logger.beDebug() && debug_hard) {
//			int spinorfield_size = sizeof(spinor) * get_spinorfieldsize(get_parameters());
//			x_tmp = create_rw_buffer(spinorfield_size);
//			const hardware::buffers::Plain<spinor>
//			int eo_spinorfield_size = sizeof(spinor) * get_eoprec_spinorfieldsize(get_parameters());
//			x_eo_tmp = create_rw_buffer(eo_spinorfield_size);
//			x_eo_tmp2 = create_rw_buffer(eo_spinorfield_size);
//
//			this->convert_from_eoprec_device(get_inout_eo(), get_tmp_eo_1(), x_tmp);
//			print_info_inv_field(get_inout_eo(), true, "\t\t\t\tX_even ");
//			print_info_inv_field(get_tmp_eo_1(), true, "\t\t\t\tX_odd ");
//			print_info_inv_field(x_tmp, false, "\t\t\t\tX = (X_even, X_odd) ");
//
//			//save x_even and x_odd temporarily
//			copy_buffer_on_device(get_inout_eo(), x_eo_tmp, eo_spinorfield_size);
//			copy_buffer_on_device(get_tmp_eo_1(), x_eo_tmp2, eo_spinorfield_size);
//			print_info_inv_field(x_eo_tmp, true, "\t\t\tx_even:\t");
//			print_info_inv_field(x_eo_tmp2, true, "\t\t\tx_odd:\t");
//		}

		//logger.debug() << "\t\tcalc eo fermion_force F(Y_even, X_odd)...";
		//Calc F(Y_even, X_odd) = F(clmem_phi_inv_eo, clmem_tmp_eo_1)
		fermion_force_eo_device(&clmem_phi_inv_eo,  fermion_code->get_tmp_eo_1(), EVEN, kappa);

		//calculate Y_odd
		//therefore, clmem_tmp_eo_1 is used as intermediate state. The result is saved in clmem_phi_inv, since
		//  this is used as a default in the force-function.
		if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
			fermion_code->dslash_eo_device(&clmem_phi_inv_eo, fermion_code->get_tmp_eo_1(), &new_u, ODD, kappa);
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_1(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			fermion_code->dslash_eo_device(&clmem_phi_inv_eo, fermion_code->get_tmp_eo_1(), &new_u, ODD, kappa);
			fermion_code->M_tm_inverse_sitediagonal_device(fermion_code->get_tmp_eo_1(), fermion_code->get_tmp_eo_2(), mubar);
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_2(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		}

		//logger.debug() << "\t\tcalc eoprec fermion_force F(Y_odd, X_even)...";
		//Calc F(Y_odd, X_even) = F(clmem_tmp_eo_1, clmem_inout_eo)
		fermion_force_eo_device(fermion_code->get_tmp_eo_1(), fermion_code->get_inout_eo(), ODD, kappa);

// will break on SOA devices
//		if(logger.beDebug() && debug_hard) {
//			int spinorfield_size = sizeof(spinor) * get_spinorfieldsize(get_parameters());
//			y_tmp = create_rw_buffer(spinorfield_size);
//			int eo_spinorfield_size = sizeof(spinor) * get_eoprec_spinorfieldsize(get_parameters());
//			y_eo_tmp = create_rw_buffer(eo_spinorfield_size);
//			y_eo_tmp2 = create_rw_buffer(eo_spinorfield_size);
//			//tmp field for differences
//			cl_mem sf_diff = create_rw_buffer(eo_spinorfield_size);
//
//			this->convert_from_eoprec_device(&clmem_phi_inv_eo, get_tmp_eo_1(), y_tmp);
//			print_info_inv_field(&clmem_phi_inv_eo, true, "\t\t\t\tY_even ");
//			print_info_inv_field(get_tmp_eo_1(), true, "\t\t\t\tY_odd ");
//			print_info_inv_field(y_tmp, false, "\t\t\t\tY = (Y_even, Yodd) ");
//
//			//save y_even and y_odd temporarily
//			copy_buffer_on_device(&clmem_phi_inv_eo, y_eo_tmp, eo_spinorfield_size);
//			copy_buffer_on_device(get_tmp_eo_1(), y_eo_tmp2, eo_spinorfield_size);
//			print_info_inv_field(y_eo_tmp, true, "\t\t\ty_even:\t");
//			print_info_inv_field(y_eo_tmp2, true, "\t\t\ty_odd:\t");
//
//			hmc_float compare_to_non_eo_diff;
//			hmc_complex compare_to_non_eo_scal;
//			hmc_complex scalprod;
//			cl_mem scal_tmp = create_rw_buffer(sizeof(hmc_complex));
//
//			logger.debug() << "\t\t\tperform some tests on eo force ingredients:";
//			//calculate scalar product
//			this->set_complex_to_scalar_product_eoprec_device(x_eo_tmp, x_eo_tmp2, scal_tmp);
//			get_buffer_from_device(scal_tmp, &scalprod, sizeof(hmc_complex));
//			logger.debug() << std::scientific << "\t\t\t(x_even, x_odd):\t" << scalprod.re << "\t" << scalprod.im;
//			//calculate difference
//			this->saxpy_eoprec_device(x_eo_tmp2, x_eo_tmp, get_clmem_one(), sf_diff);
//			print_info_inv_field(sf_diff, true, "\t\t\t(x_even - x_odd):\t");
//
//			//calculate scalar product
//			this->set_complex_to_scalar_product_eoprec_device(y_eo_tmp, y_eo_tmp2, scal_tmp);
//			get_buffer_from_device(scal_tmp, &scalprod, sizeof(hmc_complex));
//			logger.debug() << std::scientific << "\t\t\t(y_even, y_odd):\t" << scalprod.re << "\t" << scalprod.im;
//			//calculate difference
//			this->saxpy_eoprec_device(y_eo_tmp2, y_eo_tmp, get_clmem_one(), sf_diff);
//			print_info_inv_field(sf_diff, true, "\t\t\t(y_even - y_odd):\t");
//
//			//calculate scalar product
//			this->set_complex_to_scalar_product_eoprec_device(x_eo_tmp, y_eo_tmp2, scal_tmp);
//			get_buffer_from_device(scal_tmp, &scalprod, sizeof(hmc_complex));
//			logger.debug() << std::scientific << "\t\t\t(x_even, y_odd):\t" << scalprod.re << "\t" << scalprod.im;
//			//calculate difference
//			this->saxpy_eoprec_device(y_eo_tmp2, x_eo_tmp, get_clmem_one(), sf_diff);
//			print_info_inv_field(sf_diff, true, "\t\t\t(x_even - y_odd):\t");
//
//			//calculate scalar product
//			this->set_complex_to_scalar_product_eoprec_device(y_eo_tmp, x_eo_tmp2, scal_tmp);
//			get_buffer_from_device(scal_tmp, &scalprod, sizeof(hmc_complex));
//			logger.debug() << std::scientific << "\t\t\t(y_even, x_odd):\t" << scalprod.re << "\t" << scalprod.im;
//			//calculate difference
//			this->saxpy_eoprec_device(x_eo_tmp2, y_eo_tmp, get_clmem_one(), sf_diff);
//			print_info_inv_field(sf_diff, true, "\t\t\t(y_even - x_odd):\t");
//
//			//calculate scalar product
//			this->set_complex_to_scalar_product_eoprec_device(x_eo_tmp, y_eo_tmp, scal_tmp);
//			get_buffer_from_device(scal_tmp, &scalprod, sizeof(hmc_complex));
//			logger.debug() << std::scientific << "\t\t\t(x_even, y_even):\t" << scalprod.re << "\t" << scalprod.im;
//
//			compare_to_non_eo_scal.re = scalprod.re;
//			compare_to_non_eo_scal.im = scalprod.im;
//
//			//calculate difference
//			this->saxpy_eoprec_device(x_eo_tmp, y_eo_tmp, get_clmem_one(), sf_diff);
//			compare_to_non_eo_diff = print_info_inv_field(sf_diff, true, "\t\t\t(y_even - x_even):\t");
//
//			//calculate scalar product
//			this->set_complex_to_scalar_product_eoprec_device(x_eo_tmp2, y_eo_tmp2, scal_tmp);
//			get_buffer_from_device(scal_tmp, &scalprod, sizeof(hmc_complex));
//			logger.debug() << std::scientific << "\t\t\t(x_even, y_odd):\t" << scalprod.re << "\t" << scalprod.im;
//			//calculate difference
//			this->saxpy_eoprec_device(x_eo_tmp2, y_eo_tmp2, get_clmem_one(), sf_diff);
//			compare_to_non_eo_diff += print_info_inv_field(sf_diff, true, "\t\t\t(x_odd - y_odd):\t");
//
//			compare_to_non_eo_scal.re += scalprod.re;
//			compare_to_non_eo_scal.im += scalprod.im;
//
//			logger.debug() << "\t\t\tnon-eo result should be:";
//			logger.debug() << std::scientific <<  "\t\t\t(x,y):\t" << compare_to_non_eo_scal.re << " " << compare_to_non_eo_scal.im;
//			logger.debug() << std::scientific << "\t\t\t(x-y):\t" << compare_to_non_eo_diff ;
//
//			int clerr = clReleaseMemObject(scal_tmp);
//			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
//			clerr = clReleaseMemObject(sf_diff);
//			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
//
//		}

//		if(logger.beDebug() && debug_hard) {
//			logger.debug() << "\t\t\tperform checks on eo force:";
//			logger.debug() << "\t\t\tF(x_e, x_o):";
//			this->set_zero_clmem_force_device();
//			fermion_force_eo_device(x_eo_tmp, x_eo_tmp2, EVEN, kappa);
//			logger.debug() << "\t\t\tF(x_e, y_o):";
//			this->set_zero_clmem_force_device();
//			fermion_force_eo_device(x_eo_tmp, y_eo_tmp2, EVEN, kappa);
//			logger.debug() << "\t\t\tF(x_e, y_e):";
//			this->set_zero_clmem_force_device();
//			fermion_force_eo_device(x_eo_tmp, y_eo_tmp, ODD, kappa);
//			logger.debug() << "\t\t\tF(y_e, y_o):";
//			this->set_zero_clmem_force_device();
//			fermion_force_eo_device(y_eo_tmp, y_eo_tmp2, EVEN, kappa);
//			logger.debug() << "\t\t\tF(y_e, x_o):";
//			this->set_zero_clmem_force_device();
//			fermion_force_eo_device(y_eo_tmp, x_eo_tmp2, EVEN, kappa);
//			logger.debug() << "\t\t\tF(y_e, x_e):";
//			this->set_zero_clmem_force_device();
//			fermion_force_eo_device(y_eo_tmp, x_eo_tmp, ODD, kappa);
//			logger.debug() << "\t\t\tF(x_o, x_e):";
//			this->set_zero_clmem_force_device();
//			fermion_force_eo_device(x_eo_tmp2, x_eo_tmp, ODD, kappa);
//			logger.debug() << "\t\t\tF(x_o, y_e):";
//			this->set_zero_clmem_force_device();
//			fermion_force_eo_device(x_eo_tmp2, y_eo_tmp, ODD, kappa);
//			logger.debug() << "\t\t\tF(x_o, y_o):";
//			this->set_zero_clmem_force_device();
//			fermion_force_eo_device(x_eo_tmp2, y_eo_tmp2, EVEN, kappa);
//			logger.debug() << "\t\t\tF(y_o, x_e):";
//			this->set_zero_clmem_force_device();
//			fermion_force_eo_device(y_eo_tmp2, x_eo_tmp, ODD, kappa);
//			logger.debug() << "\t\t\tF(y_o, x_o):";
//			this->set_zero_clmem_force_device();
//			fermion_force_eo_device(y_eo_tmp2, x_eo_tmp2, EVEN, kappa);
//			logger.debug() << "\t\t\tF(y_o, y_e):";
//			this->set_zero_clmem_force_device();
//			fermion_force_eo_device(y_eo_tmp2, y_eo_tmp, ODD, kappa);
//			this->set_zero_clmem_force_device();
//		}

//		if(logger.beDebug() && debug_hard) {
//			//free fields from hard debugging
//			int clerr = clReleaseMemObject(y_eo_tmp);
//			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
//			clerr = clReleaseMemObject(y_eo_tmp2);
//			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
//			clerr = clReleaseMemObject(x_eo_tmp);
//			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
//			//clerr = clReleaseMemObject(x_eo_tmp2);
//			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
//		}

//		if(logger.beDebug() && debug_hard) {
//			int clerr = CL_SUCCESS;
//			hmc_complex scalprod;
//			logger.debug() << "\t\t\tperform some tests on non-eo force ingredients:";
//			print_info_inv_field(x_tmp, false, "\t\t\tx:\t");
//			print_info_inv_field(y_tmp, false, "\t\t\ty:\t");
//
//			//calculate scalar product
//			cl_mem scal_tmp = create_rw_buffer(sizeof(hmc_complex));
//			this->set_complex_to_scalar_product_device(x_tmp, y_tmp, scal_tmp);
//			get_buffer_from_device(scal_tmp, &scalprod, sizeof(hmc_complex));
//			logger.debug() << std::scientific << "\t\t\t(x,y):\t" << scalprod.re << "\t" << scalprod.im;
//			clerr = clReleaseMemObject(scal_tmp);
//
//			//calculate difference
//			this->saxpy_device(y_tmp, x_tmp, get_clmem_one(), get_inout());
//			print_info_inv_field(get_inout(), false, "\t\t\t(y-x):\t");
//
//			//calculate non-eo force
//			this->set_zero_clmem_force_device();
//			copy_buffer_on_device(x_tmp, &clmem_phi_inv, meta::get_spinorfieldsize(get_parameters()) * sizeof(spinor));
//			copy_buffer_on_device(y_tmp, get_inout(), meta::get_spinorfieldsize(get_parameters()) * sizeof(spinor));
//
//			print_info_inv_field(&clmem_phi_inv, false, "\t\t\tx before:\t");
//			print_info_inv_field(get_inout(), false, "\t\t\ty before:\t");
//
//			fermion_force_device(&clmem_phi_inv, get_inout(), kappa);
//
//			print_info_inv_field(&clmem_phi_inv, false, "\t\t\tx after:\t");
//			print_info_inv_field(get_inout(), false, "\t\t\ty after:\t");
//
//			//free fields
//			clerr = clReleaseMemObject(x_tmp);
//			clerr = clReleaseMemObject(y_tmp);
//		}
	} else {
		//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			/**
			 * The first inversion calculates
			 * X = phi = (Qplusminus)^-1 psi
			 * out of
			 * Qplusminus phi = psi
			 */
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			/**
			 * Trial solution for the spinorfield
			 */
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
			converged = fermion_code->cg(::QplusQminus(fermion_code), fermion_code->get_inout(), get_clmem_phi(), &new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

			/**
			 * Y is now just
			 *  Y = (Qminus) X = (Qminus) (Qplusminus)^-1 psi =
			 *    = (Qplus)^-1 psi
			 */
			fermion_code->Qminus(fermion_code->get_inout(), &clmem_phi_inv, &new_u, kappa, mubar);

		} else  {
			logger.debug() << "\t\tcalc fermion force ingredients using bicgstab without eo";

			/**
			 * The first inversion calculates
			 * Y = phi = (Qplus)^-1 psi
			 * out of
			 * Qplus phi = psi
			 * This is also the energy of the final field!
			 */
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			/**
			 * Trial solution for the spinorfield
			 */
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
			converged = fermion_code->bicgstab(::Qplus(fermion_code), fermion_code->get_inout(), get_clmem_phi(), &new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			hardware::buffers::copyData(&clmem_phi_inv, fermion_code->get_inout());

			/**
			 * Now, one has to calculate
			 * X = (Qminus)^-1 Y = (Qminus)^-1 (Qplus)^-1 psi = (QplusQminus)^-1 psi ??
			 * out of
			 * Qminus clmem_inout = clmem_phi_inv
			 * Therefore, when calculating the final energy of the spinorfield,
			 *  one can just take clmem_phi_inv (see also above)!!
			 */

			//copy former solution to clmem_source
			hardware::buffers::copyData(fermion_code->get_source(), fermion_code->get_inout());
			logger.debug() << "\t\t\tstart solver";

			//this sets clmem_inout cold as trial-solution
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			converged = fermion_code->bicgstab(::Qminus(fermion_code), fermion_code->get_inout(), fermion_code->get_source(), &new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

		}
		if(logger.beDebug()) {
			fermion_code->print_info_inv_field(&clmem_phi_inv, false, "\tY ");
			fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tX ");
		}
		logger.debug() << "\t\tcalc fermion_force...";
		fermion_force_device(&clmem_phi_inv, fermion_code->get_inout(), kappa);
	}
}

void Opencl_Module_Hmc::calc_fermion_force_detratio(usetimer * solvertimer, const hardware::buffers::SU3 * gaugefield)
{
	auto fermion_code = get_device()->get_fermion_code();
	auto spinor_code = get_device()->get_spinor_code();

	/**
	 * For detratio = det(kappa, mubar) / det(kappa2, mubar2) = det(Q_1^+Q_1^-) / det(Q_2^+Q_2^-)
	 * the force has almost the same ingredients as in the above case:
	 *   F(detratio) = - ( - phi^+ deriv(Q_2) X + Y^+ deriv(Q_1) X  ) + h.c.;
	 * where deriv(Q_i) is the same fct. as above with different parameters,
	 * X = (Q_1^+ Q_1^-)^-1 Q_2^+ phi
	 * Y = (Q_1^+)^-1 Q_2^+ phi
	 * (In the case of Q_2 = 1 = const., one recovers the expression for the "normal" force
	 * The main differences are:
	 *   - invert Q_2^+ phi, not phi for X and Y
	 *   - one additional force term with different mass-parameters and without Y
	 */
	int converged = -1;
	hmc_float kappa = get_parameters().get_kappa();
	hmc_float mubar = get_mubar(get_parameters());
	hmc_float kappa2 = get_parameters().get_kappa_mp();
	hmc_float mubar2 = meta::get_mubar_mp(get_parameters());
	if(get_parameters().get_use_eo() == true) {
		//CP: Init tmp spinorfield
		hardware::buffers::Spinor sf_eo_tmp(clmem_phi_eo.get_elements(), get_device());
		//the source is now Q_2^+ phi = sf_eo_tmp
		fermion_code->Qplus_eo (get_clmem_phi_eo(), &sf_eo_tmp , gaugefield, kappa2, mubar2);
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			/**
			 * The first inversion calculates
			 * X_even = phi = (Qplusminus_eo)^-1 sf_eo_tmp = (Qplusminus_eo)^-1 Q_2^+ phi
			 * out of
			 * Qplusminus_eo phi = sf_eo_tmp
			 */
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			/**
			 * Trial solution for the spinorfield
			 */
			spinor_code->set_eoprec_spinorfield_cold_device(fermion_code->get_inout_eo());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\t\t\tinv. field before inversion ");
			converged = fermion_code->cg_eo(::QplusQminus_eo(fermion_code), fermion_code->get_inout_eo(), &sf_eo_tmp, &new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\t\t\tinv. field after inversion ");

			/**
			 * Y_even is now just
			 *  Y_even = (Qminus_eo) X_even = (Qminus_eo) (Qplusminus_eo)^-1 sf_eo_tmp =
			 *    = (Qplus_eo)^-1 Q_2^+ psi
			 */
			fermion_code->Qminus_eo(fermion_code->get_inout_eo(), &clmem_phi_inv_eo, &new_u, kappa, mubar);
		} else {
			///@todo if wanted, solvertimer has to be used here..
			//logger.debug() << "\t\tcalc fermion force ingredients using bicgstab with eo.";
			/**
			 * The first inversion calculates
			 * Y_even = phi = (Qplus_eo)^-1 psi
			 * out of
			 * Qplus_eo phi = psi
			 */
			//logger.debug() << "\t\tcalc Y_even...";
			//logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			/**
			 * Trial solution for the spinorfield
			 */
			spinor_code->set_zero_spinorfield_eoprec_device(fermion_code->get_inout_eo());
			fermion_code->gamma5_eo_device(fermion_code->get_inout_eo());

			//if(logger.beDebug()) fermion_code->print_info_inv_field(get_inout_eo(), true, "\t\t\tinv. field before inversion ");
			//if(logger.beDebug()) print_info_inv_field(get_clmem_phi_eo(), true, "\t\t\tsource before inversion ");
			converged = fermion_code->bicgstab_eo(::Qplus_eo(fermion_code), fermion_code->get_inout_eo(), &sf_eo_tmp, &new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			//if(logger.beDebug()) print_info_inv_field(get_inout_eo(), true, "\t\t\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			hardware::buffers::copyData(&clmem_phi_inv_eo, fermion_code->get_inout_eo());

			/**
			 * Now, one has to calculate
			 * X = (Qminus_eo)^-1 Y = (Qminus_eo)^-1 (Qplus_eo)^-1 sf_eo_tmp = (QplusQminus_eo)^-1 sf_eo_tmp = (QplusQminus_eo)^-1 Q_2^+ psi
			 * out of
			 * Qminus_eo clmem_inout_eo = clmem_phi_inv_eo
			 */

			//logger.debug() << "\t\tcalc X_even...";
			//copy former solution to clmem_source
			hardware::buffers::copyData(fermion_code->get_source_even(), fermion_code->get_inout_eo());

			//this sets clmem_inout cold as trial-solution
			spinor_code->set_eoprec_spinorfield_cold_device(fermion_code->get_inout_eo());

			//if(logger.beDebug()) print_info_inv_field(get_inout_eo(), true, "\t\t\tinv. field before inversion ");
			//if(logger.beDebug()) print_info_inv_field(get_source_even(), true, "\t\t\tsource before inversion ");
			converged = fermion_code->bicgstab_eo(::Qminus_eo(fermion_code), fermion_code->get_inout_eo(), fermion_code->get_source_even(), &new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			//if(logger.beDebug()) print_info_inv_field(get_inout_eo(), true, "\t\t\tinv. field after inversion ");
		}
		/**
		 * At this point, one has to calculate X_odd and Y_odd.
		 * If one has a fermionmatrix
		 *  M = R + D
		 * these are:
		 *  X_odd = -R(-mu)_inv D X_even
		 *  Y_odd = -R(mu)_inv D Y_even
		 */

		/** @fixme below usages of dslash should work, but only because we use bicgstab above
		proper implementation needs to make sure this is always the case */

		///@NOTE the following calculations could also go in a new function for convenience
		//calculate X_odd
		//therefore, clmem_tmp_eo_1 is used as intermediate state. The result is saved in clmem_inout, since
		//  this is used as a default in the force-function.
		if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
			fermion_code->dslash_eo_device(fermion_code->get_inout_eo(), fermion_code->get_tmp_eo_1(), &new_u, ODD);
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_1(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			fermion_code->dslash_eo_device(fermion_code->get_inout_eo(), fermion_code->get_tmp_eo_1(), &new_u, ODD);
			fermion_code->M_tm_inverse_sitediagonal_minus_device(fermion_code->get_tmp_eo_1(), fermion_code->get_tmp_eo_2());
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_2(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		}

		//logger.debug() << "\t\tcalc eo fermion_force F(Y_even, X_odd)...";
		//Calc F(Y_even, X_odd) = F(clmem_phi_inv_eo, clmem_tmp_eo_1)
		fermion_force_eo_device(&clmem_phi_inv_eo,  fermion_code->get_tmp_eo_1(), EVEN, kappa);

		//calculate Y_odd
		//therefore, clmem_tmp_eo_1 is used as intermediate state. The result is saved in clmem_phi_inv, since
		//  this is used as a default in the force-function.
		if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
			fermion_code->dslash_eo_device(&clmem_phi_inv_eo, fermion_code->get_tmp_eo_1(), &new_u, ODD);
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_1(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			fermion_code->dslash_eo_device(&clmem_phi_inv_eo, fermion_code->get_tmp_eo_1(), &new_u, ODD);
			fermion_code->M_tm_inverse_sitediagonal_device(fermion_code->get_tmp_eo_1(), fermion_code->get_tmp_eo_2());
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_2(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		}

		//logger.debug() << "\t\tcalc eoprec fermion_force F(Y_odd, X_even)...";
		//Calc F(Y_odd, X_even) = F(clmem_tmp_eo_1, clmem_inout_eo)
		fermion_force_eo_device(fermion_code->get_tmp_eo_1(), fermion_code->get_inout_eo(), ODD, kappa);

		/**
		 *Now, one has the additional term - phi^+ deriv(Q_2) X
		 *This works exactly the same as above, except replacing Y -> -phi and different mass parameters
		 */

		//Y is not needed anymore, therefore use clmem_phi_inv_eo to store -phi
		spinor_code->sax_eoprec_device(get_clmem_phi_eo(), fermion_code->get_clmem_minusone(), &clmem_phi_inv_eo);

		//logger.debug() << "\t\tcalc eo fermion_force F(Y_even, X_odd)...";
		//Calc F(Y_even, X_odd) = F(clmem_phi_inv_eo, clmem_tmp_eo_1)
		fermion_force_eo_device(&clmem_phi_inv_eo,  fermion_code->get_tmp_eo_1(), EVEN, kappa2);

		//calculate phi_odd
		//this works in the same way as with Y above, since -phi_even is saved in the same buffer as Y_even
		//therefore, clmem_tmp_eo_1 is used as intermediate state. The result is saved in clmem_phi_inv, since
		//  this is used as a default in the force-function.
		if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
			fermion_code->dslash_eo_device(&clmem_phi_inv_eo, fermion_code->get_tmp_eo_1(), &new_u, ODD);
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_1(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			fermion_code->dslash_eo_device(&clmem_phi_inv_eo, fermion_code->get_tmp_eo_1(), &new_u, ODD);
			fermion_code->M_tm_inverse_sitediagonal_device(fermion_code->get_tmp_eo_1(), fermion_code->get_tmp_eo_2());
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_2(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		}

		//logger.debug() << "\t\tcalc eoprec fermion_force F(Y_odd, X_even)...";
		//Calc F(Y_odd, X_even) = F(clmem_tmp_eo_1, clmem_inout_eo)
		fermion_force_eo_device(fermion_code->get_tmp_eo_1(), fermion_code->get_inout_eo(), ODD, kappa2);
	} else {
		//CP: Init tmp spinorfield
		hardware::buffers::Plain<spinor> sf_tmp(clmem_phi.get_elements(), get_device());
		//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
		//the source is now Q_2^+ phi = sf_tmp
		fermion_code->Qplus (get_clmem_phi(), &sf_tmp , gaugefield, kappa2, mubar2);
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			/**
			 * The first inversion calculates
			 * X = phi = (Qplusminus)^-1 sf_tmp
			 * out of
			 * Qplusminus phi = sf_tmp
			 */
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			/**
			 * Trial solution for the spinorfield
			 */
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
			converged = fermion_code->cg(::QplusQminus(fermion_code), fermion_code->get_inout(), &sf_tmp, &new_u, get_parameters().get_force_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

			/**
			 * Y is now just
			 *  Y = (Qminus) X = (Qminus) (Qplusminus)^-1 sf_tmp =
			 *    = (Qplus)^-1 sf_tmp
			 */
			fermion_code->Qminus(fermion_code->get_inout(), &clmem_phi_inv, &new_u);

		} else  {
			logger.debug() << "\t\tcalc fermion force ingredients using bicgstab without eo";

			/**
			 * The first inversion calculates
			 * Y = phi = (Qplus)^-1 sf_tmp
			 * out of
			 * Qplus phi = sf_tmp
			 * This is also the energy of the final field!
			 */
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			/**
			 * Trial solution for the spinorfield
			 */
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
			converged = fermion_code->bicgstab(::Qplus(fermion_code), fermion_code->get_inout(), &sf_tmp, &new_u, get_parameters().get_force_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			hardware::buffers::copyData(&clmem_phi_inv, fermion_code->get_inout());

			/**
			 * Now, one has to calculate
			 * X = (Qminus)^-1 Y = (Qminus)^-1 (Qplus)^-1 psi = (QplusQminus)^-1 psi ??
			 * out of
			 * Qminus clmem_inout = clmem_phi_inv
			 * Therefore, when calculating the final energy of the spinorfield,
			 *  one can just take clmem_phi_inv (see also above)!!
			 */

			//copy former solution to clmem_source
			hardware::buffers::copyData(fermion_code->get_source(), fermion_code->get_inout());
			logger.debug() << "\t\t\tstart solver";

			//this sets clmem_inout cold as trial-solution
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			converged = fermion_code->bicgstab(::Qminus(fermion_code), fermion_code->get_inout(), fermion_code->get_source(), &new_u, get_parameters().get_force_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

		}
		if(logger.beDebug()) {
			fermion_code->print_info_inv_field(&clmem_phi_inv, false, "\tY ");
			fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tX ");
		}
		logger.debug() << "\t\tcalc fermion_force...";
		fermion_force_device(&clmem_phi_inv, fermion_code->get_inout(), kappa);

		/**
		 *Now, one has the additional term - phi^+ deriv(Q_2) X
		 *This works exactly the same as above, except replacing Y -> -phi and different mass parameters
		 */

		//Y is not needed anymore, therefore use clmem_phi_inv_eo to store -phi
		spinor_code->sax_device(get_clmem_phi(), fermion_code->get_clmem_minusone(), &clmem_phi_inv);

		fermion_force_device(&clmem_phi_inv, fermion_code->get_inout(), kappa2);
	}
}

void Opencl_Module_Hmc::calc_gauge_force()
{
	logger.debug() << "\t\tcalc gauge_force...";
	gauge_force_device();
	if(meta::get_use_rectangles(get_parameters()) == true) {
		logger.debug() << "\t\tcalc rect gauge_force...";
		gauge_force_tlsym_device();
	}
}

hmc_float Opencl_Module_Hmc::calc_s_fermion()
{
	auto fermion_code = get_device()->get_fermion_code();
	auto spinor_code = get_device()->get_spinor_code();

	logger.debug() << "calc final fermion energy...";
	//this function essentially performs the same steps as in the force-calculation, but with higher precision.
	//  therefore, comments are deleted here...
	//  Furthermore, in the bicgstab-case, the second inversions are not needed
	int converged = -1;
	if(get_parameters().get_use_eo() == true) {
		//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			logger.debug() << "\t\t\tstart solver";

			spinor_code->set_eoprec_spinorfield_cold_device(fermion_code->get_inout_eo());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\tinv. field before inversion ");
			converged = fermion_code->cg_eo(::QplusQminus_eo(fermion_code), fermion_code->get_inout_eo(), get_clmem_phi_eo(), &new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\tinv. field after inversion ");

			fermion_code->Qminus_eo(fermion_code->get_inout_eo(), &clmem_phi_inv_eo, &new_u);
		} else {
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			  * or a point-source spinorfield as trial-solution
			  */
			spinor_code->set_zero_spinorfield_eoprec_device(fermion_code->get_inout_eo());
			fermion_code->gamma5_eo_device(fermion_code->get_inout_eo());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\tinv. field before inversion ");
			if(logger.beDebug()) fermion_code->print_info_inv_field(get_clmem_phi_eo(), true, "\tsource before inversion ");
			converged = fermion_code->bicgstab_eo(::Qplus_eo(fermion_code), fermion_code->get_inout_eo(), get_clmem_phi_eo(), &new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			hardware::buffers::copyData(&clmem_phi_inv_eo, fermion_code->get_inout_eo());
		}
	} else {
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			  * or a point-source spinorfield as trial-solution
			  */
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			converged = fermion_code->cg(::QplusQminus(fermion_code), fermion_code->get_inout(), get_clmem_phi(), &new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

			fermion_code->Qminus(fermion_code->get_inout(), &clmem_phi_inv, &new_u);

		} else  {

			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			  * or a point-source spinorfield as trial-solution
			  */
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			converged = fermion_code->bicgstab(::Qplus(fermion_code), fermion_code->get_inout(), get_clmem_phi(), &new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			hardware::buffers::copyData(&clmem_phi_inv, fermion_code->get_inout());
		}
	}
	///@todo: this can be moved in the ifs above!!
	if(get_parameters().get_use_eo() == true) {
		get_device()->get_spinor_code()->set_float_to_global_squarenorm_eoprec_device(&clmem_phi_inv_eo, &clmem_s_fermion);
	} else {
		get_device()->get_spinor_code()->set_float_to_global_squarenorm_device(&clmem_phi_inv, &clmem_s_fermion);
	}
	hmc_float tmp;
	clmem_s_fermion.dump(&tmp);
	return tmp;
}

hmc_float Opencl_Module_Hmc::calc_s_fermion_mp(const hardware::buffers::SU3 * gaugefield)
{
	auto fermion_code = get_device()->get_fermion_code();
	auto spinor_code = get_device()->get_spinor_code();

	logger.debug() << "calc final fermion energy...";
	//this function essentially performs the same steps as in the non mass-prec case, however, one has to apply one more matrix multiplication
	//  therefore, comments are deleted here...
	//  Furthermore, in the bicgstab-case, the second inversions are not needed
	int converged = -1;
	if(get_parameters().get_use_eo() == true) {
		//CP: Init tmp spinorfield
		hardware::buffers::Spinor sf_eo_tmp(clmem_phi_mp_eo.get_elements(), get_device());

		//sf_eo_tmp = Qplus_eo(heavy_mass) phi_mp_eo
		fermion_code->Qplus_eo (get_clmem_phi_mp_eo(), &sf_eo_tmp , gaugefield, get_parameters().get_kappa_mp(), meta::get_mubar_mp(get_parameters()));
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			logger.debug() << "\t\t\tstart solver";
			spinor_code->set_eoprec_spinorfield_cold_device(fermion_code->get_inout_eo());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\tinv. field before inversion ");
			converged = fermion_code->cg_eo(::QplusQminus_eo(fermion_code), fermion_code->get_inout_eo(), &sf_eo_tmp, &new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\tinv. field after inversion ");

			fermion_code->Qminus_eo(fermion_code->get_inout_eo(), &clmem_phi_inv_eo, &new_u);
		} else {
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			spinor_code->set_zero_spinorfield_eoprec_device(fermion_code->get_inout_eo());
			fermion_code->gamma5_eo_device(fermion_code->get_inout_eo());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\tinv. field before inversion ");
			if(logger.beDebug()) fermion_code->print_info_inv_field(&sf_eo_tmp, true, "\tsource before inversion ");
			converged = fermion_code->bicgstab_eo(::Qplus_eo(fermion_code), fermion_code->get_inout_eo(), &sf_eo_tmp, &new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			hardware::buffers::copyData(&clmem_phi_inv_eo, fermion_code->get_inout_eo());
		}
	} else {
		//CP: Init tmp spinorfield
		hardware::buffers::Plain<spinor> sf_tmp(clmem_phi_mp.get_elements(), get_device());

		//sf_eo_tmp = Qplus(light_mass) phi_mp
		fermion_code->Qplus (get_clmem_phi_mp(), &sf_tmp , gaugefield, get_parameters().get_kappa_mp(), meta::get_mubar_mp(get_parameters()));
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			converged = fermion_code->cg(::QplusQminus(fermion_code), fermion_code->get_inout(), &sf_tmp, &new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

			fermion_code->Qminus(fermion_code->get_inout(), &clmem_phi_inv, &new_u);

		} else  {

			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			converged = fermion_code->bicgstab(::Qplus(fermion_code), fermion_code->get_inout(), &sf_tmp, &new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			hardware::buffers::copyData(&clmem_phi_inv, fermion_code->get_inout());
		}
	}
	///@todo: this can be moved in the ifs above!!
	if(get_parameters().get_use_eo() == true) {
		get_device()->get_spinor_code()->set_float_to_global_squarenorm_eoprec_device(&clmem_phi_inv_eo, &clmem_s_fermion);
	} else {
		get_device()->get_spinor_code()->set_float_to_global_squarenorm_device(&clmem_phi_inv, &clmem_s_fermion);
	}
	hmc_float tmp;
	clmem_s_fermion.dump(&tmp);
	return tmp;
}

hmc_observables Opencl_Module_Hmc::metropolis(hmc_float rnd, hmc_float beta, const hardware::buffers::SU3 * gaugefield)
{
	auto gf_code = get_device()->get_gaugefield_code();
	auto fermion_code = get_device()->get_fermion_code();
	auto gm_code = get_device()->get_gaugemomentum_code();

	//Calc Hamiltonian
	logger.debug() << "Calculate Hamiltonian";
	hmc_float deltaH = 0.;
	hmc_float s_old = 0.;
	hmc_float s_new = 0.;

	//Gauge-Part
	hmc_float tplaq, splaq, plaq;
	hmc_float tplaq_new, splaq_new, plaq_new;
	hmc_float rect_new = 0.;
	hmc_float rect = 0.;
	hmc_complex poly;
	hmc_complex poly_new;
	//In this call, the observables are calculated already with appropiate Weighting factor of 2.0/(VOL4D*NDIM*(NDIM-1)*NC)
	gf_code->gaugeobservables(gaugefield, &plaq,  &tplaq, &splaq, &poly);
	gf_code->gaugeobservables(&new_u, &plaq_new,  &tplaq_new, &splaq_new, &poly_new);
	//plaq has to be divided by the norm-factor to get s_gauge
	hmc_float factor = 1. / (meta::get_plaq_norm(get_parameters()));
	if(meta::get_use_rectangles(get_parameters()) == true) {
		gf_code->gaugeobservables_rectangles(gaugefield, &rect);
		gf_code->gaugeobservables_rectangles(&new_u, &rect_new);
		hmc_float c0 = meta::get_c0(get_parameters());
		hmc_float c1 = meta::get_c1(get_parameters());
		deltaH = - beta * ( c0 * (plaq - plaq_new) / factor + c1 * ( rect - rect_new )  );
		s_old = - beta * ( c0 * (plaq) / factor + c1 * ( rect )  );
		s_new = - beta * ( c0 * (plaq_new) / factor + c1 * ( rect_new )  );

	} else {
		/** NOTE: the minus here is introduced to fit tmlqcd!!! */
		deltaH = -(plaq - plaq_new) * beta / factor;
		s_old = -(plaq ) * beta / factor;
		s_new = -(plaq_new) * beta / factor;
	}

	logger.debug() << "\tS_gauge(old field) = " << setprecision(10) << s_old;
	logger.debug() << "\tS_gauge(new field) = " << setprecision(10) << s_new;
	logger.info() << "\tdeltaS_gauge = " << setprecision(10) << deltaH;

	//Gaugemomentum-Part
	hmc_float p2, new_p2;
	gm_code->set_float_to_gaugemomentum_squarenorm_device(&clmem_p, &clmem_p2);
	gm_code->set_float_to_gaugemomentum_squarenorm_device(&clmem_new_p, &clmem_new_p2);
	clmem_p2.dump(&p2);
	clmem_new_p2.dump(&new_p2);
	//the energy is half the squarenorm
	deltaH += 0.5 * (p2 - new_p2);

	logger.debug() << "\tS_gaugemom(old field) = " << setprecision(10) << 0.5 * p2;
	logger.debug() << "\tS_gaugemom(new field) = " << setprecision(10) << 0.5 * new_p2;
	logger.info() << "\tdeltaS_gaugemom = " << setprecision(10) << 0.5 * (p2 - new_p2);

	//Fermion-Part:
	if(! get_parameters().get_use_gauge_only() ) {
		hmc_float spinor_energy_init, s_fermion_final;
		//initial energy has been computed in the beginning...
		clmem_s_fermion_init.dump(&spinor_energy_init);
		// sum_links phi*_i (M^+M)_ij^-1 phi_j
		s_fermion_final = calc_s_fermion();
		deltaH += spinor_energy_init - s_fermion_final;

		logger.debug() << "\tS_ferm(old field) = " << setprecision(10) <<  spinor_energy_init;
		logger.debug() << "\tS_ferm(new field) = " << setprecision(10) << s_fermion_final;
		logger.info() << "\tdeltaS_ferm = " << spinor_energy_init - s_fermion_final;
		if( get_parameters().get_use_mp() ) {
			hmc_float spinor_energy_mp_init, s_fermion_mp_final;
			//initial energy has been computed in the beginning...
			clmem_s_fermion_mp_init.dump(&spinor_energy_mp_init);
			// sum_links phi*_i (M^+M)_ij^-1 phi_j
			s_fermion_mp_final = calc_s_fermion_mp(gaugefield);
			deltaH += spinor_energy_mp_init - s_fermion_mp_final;

			logger.debug() << "\tS_ferm_mp(old field) = " << setprecision(10) <<  spinor_energy_mp_init;
			logger.debug() << "\tS_ferm_mp(new field) = " << setprecision(10) << s_fermion_mp_final;
			logger.info() << "\tdeltaS_ferm_mp = " << spinor_energy_init - s_fermion_mp_final;
		}
	}
	//Metropolis-Part
	hmc_float compare_prob;
	if(deltaH < 0) {
		compare_prob = exp(deltaH);
	} else {
		compare_prob = 1.0;
	}
	logger.info() << "\tdeltaH = " << deltaH << "\tAcc-Prop = " << compare_prob;
	hmc_observables tmp;
	if(rnd <= compare_prob) {
		tmp.accept = 1;
		tmp.plaq = plaq_new;
		tmp.tplaq = tplaq_new;
		tmp.splaq = splaq_new;
		tmp.poly = poly_new;
		tmp.deltaH = deltaH;
		tmp.prob = compare_prob;
		if(meta::get_use_rectangles(get_parameters()) ) tmp.rectangles = rect_new / get_rect_norm(get_parameters());
	} else {
		tmp.accept = 0;
		tmp.plaq = plaq;
		tmp.tplaq = tplaq;
		tmp.splaq = splaq;
		tmp.poly = poly;
		tmp.deltaH = deltaH;
		tmp.prob = compare_prob;
		if(meta::get_use_rectangles(get_parameters()) ) tmp.rectangles = rect / get_rect_norm(get_parameters());
	}

	return tmp;
}

void Opencl_Module_Hmc::calc_spinorfield_init_energy(hardware::buffers::Plain<hmc_float> * dest)
{
	auto spinor_code = get_device()->get_spinor_code();

	//Suppose the initial spinorfield is saved in phi_inv
	//  it is created in generate_gaussian_spinorfield_device
	if(get_parameters().get_use_eo() == true) {
		spinor_code->set_float_to_global_squarenorm_eoprec_device(&clmem_phi_inv_eo, dest);
	} else {
		spinor_code->set_float_to_global_squarenorm_device(&clmem_phi_inv, dest);
	}
}

void Opencl_Module_Hmc::md_update_gaugemomentum_device(hmc_float eps)
{
	md_update_gaugemomentum_device(&clmem_force, &clmem_new_p, eps);

	if(logger.beDebug()) {
	  auto gm_code = get_device()->get_gaugemomentum_code();
		hardware::buffers::Plain<hmc_float> force_tmp(1, get_device());
		hmc_float resid;
		gm_code->set_float_to_gaugemomentum_squarenorm_device(&clmem_new_p, &force_tmp);
		force_tmp.dump(&resid);
		logger.debug() <<  "\tupdated gaugemomenta energy:\t" << resid;
		if(resid != resid) {
			throw Print_Error_Message("calculation of gm gave nan! Aborting...", __FILE__, __LINE__);
		}
	}
}

void Opencl_Module_Hmc::md_update_gaugemomentum_device(const hardware::buffers::Gaugemomentum * in, const hardware::buffers::Gaugemomentum * out, hmc_float eps)
{
	//__kernel void md_update_gaugemomenta(hmc_float eps, __global ae * p_inout, __global ae* force_in){
	hmc_float tmp = eps;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(md_update_gaugemomenta, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(md_update_gaugemomenta, 0, sizeof(hmc_float), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(md_update_gaugemomenta, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(md_update_gaugemomenta, 2, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(md_update_gaugemomenta , gs2, ls2);
}

void Opencl_Module_Hmc::md_update_gaugefield_device(hmc_float eps)
{
	md_update_gaugefield_device(&clmem_new_p, &new_u, eps);
}

void Opencl_Module_Hmc::md_update_gaugefield_device(const hardware::buffers::Gaugemomentum * gm_in, const hardware::buffers::SU3 * gf_out, hmc_float eps)
{
	// __kernel void md_update_gaugefield(hmc_float eps, __global ae * p_in, __global ocl_s_gaugefield * u_inout){
	hmc_float tmp = eps;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(md_update_gaugefield, &ls2, &gs2, &num_groups);
	//set arguments
	//this is always applied to clmem_force
	int clerr = clSetKernelArg(md_update_gaugefield, 0, sizeof(hmc_float), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(md_update_gaugefield, 1, sizeof(cl_mem), gm_in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(md_update_gaugefield, 2, sizeof(cl_mem), gf_out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( md_update_gaugefield , gs2, ls2);
}

void Opencl_Module_Hmc::set_zero_clmem_force_device()
{
  auto gm_code = get_device()->get_gaugemomentum_code();
	gm_code->set_zero_gaugemomentum(&clmem_force);
}

void Opencl_Module_Hmc::set_zero_gaugemomentum(const hardware::buffers::Gaugemomentum * buf)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(_set_zero_gaugemomentum, &ls2, &gs2, &num_groups);
	//set arguments
	//this is always applied to clmem_force
	int clerr = clSetKernelArg(_set_zero_gaugemomentum, 0, sizeof(cl_mem), buf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(_set_zero_gaugemomentum , gs2, ls2);
}

void Opencl_Module_Hmc::gauge_force_device()
{
	gauge_force_device(&new_u, &clmem_force);

	if(logger.beDebug()) {
	  	  auto gm_code = get_device()->get_gaugemomentum_code();
		hardware::buffers::Plain<hmc_float> gauge_force_tmp(1, get_device());
		hmc_float gauge_force_energy = 0.;
		gm_code->set_float_to_gaugemomentum_squarenorm_device(&clmem_force, &gauge_force_tmp);
		gauge_force_tmp.dump(&gauge_force_energy);

		logger.debug() <<  "\t\t\tgauge force:\t" << gauge_force_energy;

		if(gauge_force_energy != gauge_force_energy) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
	}

	//recalculate force with local buffer to get only this contribution to the force vector
	if(logger.beDebug()) {
		const hardware::buffers::Gaugemomentum force2(clmem_force.get_elements(), clmem_force.get_device());
		auto gm_code = get_device()->get_gaugemomentum_code();
		//init new buffer to zero
		gm_code->set_zero_gaugemomentum(&force2);

		//re-calculate force
		gauge_force_device(&new_u, &force2);

		hardware::buffers::Plain<hmc_float> check_force_tmp(1, get_device());
		hmc_float check_force_energy = 0.;
		gm_code->set_float_to_gaugemomentum_squarenorm_device(&force2, &check_force_tmp);
		check_force_tmp.dump(&check_force_energy);
		//logger.debug() <<  "\t\t\t\tforce contribution:\t" << check_force_energy;
		if(check_force_energy != check_force_energy) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
	}

}

void Opencl_Module_Hmc::gauge_force_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Gaugemomentum * out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gauge_force, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(gauge_force, 0, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(gauge_force, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( gauge_force , gs2, ls2);
}


void Opencl_Module_Hmc::gauge_force_tlsym_device()
{
	gauge_force_tlsym_device(&new_u, &clmem_force);

	if(logger.beDebug()) {
		hardware::buffers::Plain<hmc_float> gauge_force_tlsym_tmp(1, get_device());
		auto gm_code = get_device()->get_gaugemomentum_code();
		hmc_float gauge_force_tlsym_energy = 0.;
		gm_code->set_float_to_gaugemomentum_squarenorm_device(&clmem_force, &gauge_force_tlsym_tmp);
		gauge_force_tlsym_tmp.dump(&gauge_force_tlsym_energy);

		logger.debug() <<  "\t\t\tgauge force tlsym:\t" << gauge_force_tlsym_energy;

		if(gauge_force_tlsym_energy != gauge_force_tlsym_energy) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
	}

	//recalculate force with local buffer to get only this contribution to the force vector
	if(logger.beDebug()) {
		const hardware::buffers::Gaugemomentum force2(clmem_force.get_elements(), clmem_force.get_device());
		auto gm_code = get_device()->get_gaugemomentum_code();
		//init new buffer to zero
		gm_code->set_zero_gaugemomentum(&force2);

		gauge_force_tlsym_device(&new_u, &force2);

		hardware::buffers::Plain<hmc_float> check_force_tmp(1, get_device());
		hmc_float check_force_energy = 0.;
		gm_code->set_float_to_gaugemomentum_squarenorm_device(&force2, &check_force_tmp);
		check_force_tmp.dump(&check_force_energy);
		logger.debug() <<  "\t\t\t\tforce contribution:\t" << check_force_energy;
		if(check_force_energy != check_force_energy) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
	}

}

void Opencl_Module_Hmc::gauge_force_tlsym_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Gaugemomentum * out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gauge_force_tlsym, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(gauge_force_tlsym, 0, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(gauge_force_tlsym, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( gauge_force_tlsym , gs2, ls2);
}

void Opencl_Module_Hmc::fermion_force_device(const hardware::buffers::Plain<spinor> * Y, const hardware::buffers::Plain<spinor> * X, hmc_float kappa)
{
	using namespace hardware::buffers;

	fermion_force_device(Y, X, &new_u, &clmem_force, kappa);

	if(logger.beDebug()) {
		Plain<hmc_float> noneo_force_tmp(1, get_device());
		auto gm_code = get_device()->get_gaugemomentum_code();
		hmc_float noneo_force_energy = 0.;
		gm_code->set_float_to_gaugemomentum_squarenorm_device(&clmem_force, &noneo_force_tmp);
		noneo_force_tmp.dump(&noneo_force_energy);
		//logger.debug() <<  "\t\t\tnon-eo force:\t" << noneo_force_energy;
		if(noneo_force_energy != noneo_force_energy) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
	}

	//recalculate force with local buffer to get only this contribution to the force vector
	if(logger.beDebug()) {
		const hardware::buffers::Gaugemomentum force2(clmem_force.get_elements(), clmem_force.get_device());
		auto gm_code = get_device()->get_gaugemomentum_code();
		//init new buffer to zero
		gm_code->set_zero_gaugemomentum(&force2);

		//re-calculate force
		fermion_force_device(Y, X, &new_u, &force2, kappa);

		Plain<hmc_float> check_force_tmp(1, get_device());
		hmc_float check_force_energy = 0.;
		gm_code->set_float_to_gaugemomentum_squarenorm_device(&force2, &check_force_tmp);
		check_force_tmp.dump(&check_force_energy);
		//logger.debug() <<  "\t\t\t\tforce contribution:\t" << check_force_energy;
		if(check_force_energy != check_force_energy) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
	}

}

void Opencl_Module_Hmc::fermion_force_device(const hardware::buffers::Plain<spinor> * Y, const hardware::buffers::Plain<spinor> * X, const hardware::buffers::SU3 * gf, const hardware::buffers::Gaugemomentum * out, hmc_float kappa)
{
	//get kappa
	hmc_float kappa_tmp;
	if(kappa == ARG_DEF) kappa_tmp = get_parameters().get_kappa();
	else kappa_tmp = kappa;

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(fermion_force, &ls2, &gs2, &num_groups);
	//set arguments
	//fermion_force(field, Y, X, out);
	int clerr = clSetKernelArg(fermion_force, 0, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 1, sizeof(cl_mem), Y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 2, sizeof(cl_mem), X->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 3, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 4, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( fermion_force , gs2, ls2);
}

//the argument kappa is set to ARG_DEF as default
void Opencl_Module_Hmc::fermion_force_eo_device(const hardware::buffers::Spinor * Y, const hardware::buffers::Spinor * X, int evenodd, hmc_float kappa)
{
	using namespace hardware::buffers;

	fermion_force_eo_device(Y, X, &new_u, &clmem_force, evenodd, kappa);

	if(logger.beDebug()) {
		Plain<hmc_float> force_tmp(1, get_device());
		auto gm_code = get_device()->get_gaugemomentum_code();
		hmc_float resid;
		gm_code->set_float_to_gaugemomentum_squarenorm_device(&clmem_force, &force_tmp);
		force_tmp.dump(&resid);
		//logger.debug() <<  "\t\t\teoprec force:\t" << resid;

		if(resid != resid) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
	}

	//recalculate force with local buffer, giving only this contribution to the force
	if(logger.beDebug()) {
		const hardware::buffers::Gaugemomentum force2(clmem_force.get_elements(), clmem_force.get_device());
		auto gm_code = get_device()->get_gaugemomentum_code();
		//init new buffer to zero
		gm_code->set_zero_gaugemomentum(&force2);

		//re-calculate force
		fermion_force_eo_device(Y, X, &new_u, &force2, evenodd, kappa);

		Plain<hmc_float> check_force_tmp(1, get_device());
		hmc_float check_force_energy = 0.;
		gm_code->set_float_to_gaugemomentum_squarenorm_device(&force2, &check_force_tmp);
		check_force_tmp.dump(&check_force_energy);
		//logger.debug() <<  "\t\t\t\tforce contribution:\t" << check_force_energy;
		if(check_force_energy != check_force_energy) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
	}
}

//the argument kappa is set to ARG_DEF as default
void Opencl_Module_Hmc::fermion_force_eo_device(const hardware::buffers::Spinor * Y, const hardware::buffers::Spinor * X, const hardware::buffers::SU3 * gf, const hardware::buffers::Gaugemomentum * out, int evenodd, hmc_float kappa)
{
	//get kappa
	hmc_float kappa_tmp;
	if(kappa == ARG_DEF) kappa_tmp = get_parameters().get_kappa();
	else kappa_tmp = kappa;

	//fermion_force(field, Y, X, out);
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(fermion_force_eo, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(fermion_force_eo, 0, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force_eo, 1, sizeof(cl_mem), Y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force_eo, 2, sizeof(cl_mem), X->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force_eo, 3, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force_eo, 4, sizeof(int), &evenodd);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force_eo, 5, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( fermion_force_eo , gs2, ls2);
}

void Opencl_Module_Hmc::stout_smeared_fermion_force_device(std::vector<const hardware::buffers::SU3 *>& gf_intermediate)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(stout_smear_fermion_force, &ls2, &gs2, &num_groups);
	//set arguments
}

void Opencl_Module_Hmc::set_float_to_gaugemomentum_squarenorm_device(const hardware::buffers::Gaugemomentum * clmem_in, const hardware::buffers::Plain<hmc_float> * out)
{
	auto spinor_code = get_device()->get_spinor_code();

	//__kernel void gaugemomentum_squarenorm(__global ae * in, __global hmc_float * out){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gaugemomentum_squarenorm, &ls2, &gs2, &num_groups);

	const hardware::buffers::Plain<hmc_float> clmem_global_squarenorm_buf_glob(num_groups, get_device());

	int clerr = clSetKernelArg(gaugemomentum_squarenorm, 0, sizeof(cl_mem), clmem_in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(gaugemomentum_squarenorm, 1, sizeof(cl_mem), clmem_global_squarenorm_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(gaugemomentum_squarenorm, 2, sizeof(hmc_float) * ls2, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	get_device()->enqueue_kernel(gaugemomentum_squarenorm, gs2, ls2);

	spinor_code->global_squarenorm_reduction(out, &clmem_global_squarenorm_buf_glob);
}

void Opencl_Module_Hmc::importGaugemomentumBuffer(const hardware::buffers::Gaugemomentum * dest, const ae * const data)
{
	cl_int clerr;
	if(dest->is_soa()) {
		hardware::buffers::Plain<ae> tmp(meta::get_vol4d(get_parameters()) * NDIM, dest->get_device());
		tmp.load(data);

		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(gaugemomentum_convert_to_soa, &ls2, &gs2, &num_groups);
		clerr = clSetKernelArg(gaugemomentum_convert_to_soa, 0, sizeof(cl_mem), dest->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		clerr = clSetKernelArg(gaugemomentum_convert_to_soa, 1, sizeof(cl_mem), tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		get_device()->enqueue_kernel(gaugemomentum_convert_to_soa, gs2, ls2);
	} else {
		dest->load(data);
	}
}

void Opencl_Module_Hmc::exportGaugemomentumBuffer(ae * const dest, const hardware::buffers::Gaugemomentum * buf)
{
	cl_int clerr;
	if(buf->is_soa()) {
		hardware::buffers::Plain<ae> tmp(meta::get_vol4d(get_parameters()) * NDIM, buf->get_device());

		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(gaugemomentum_convert_from_soa, &ls2, &gs2, &num_groups);
		clerr = clSetKernelArg(gaugemomentum_convert_from_soa, 0, sizeof(cl_mem), tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		clerr = clSetKernelArg(gaugemomentum_convert_from_soa, 1, sizeof(cl_mem), buf->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		get_device()->enqueue_kernel(gaugemomentum_convert_from_soa, gs2, ls2);

		tmp.dump(dest);
	} else {
		buf->dump(dest);
	}
}

Opencl_Module_Hmc::Opencl_Module_Hmc(const meta::Inputparameters& params, hardware::Device * device)
	: Opencl_Module(params, device),
	  clmem_s_fermion_init(1, device),
	  clmem_s_fermion_mp_init(1, device),
	  clmem_p2(1, device),
	  clmem_new_p2(1, device),
	  clmem_s_fermion(1, device),
	  clmem_p(meta::get_vol4d(params) * NDIM, device),
	  clmem_new_p(meta::get_vol4d(params) * NDIM, device),
	  new_u(device->get_gaugefield_code()->get_gaugefield()->get_elements(), device),
	  clmem_force(meta::get_vol4d(params) * NDIM, device),
	  clmem_phi_inv(meta::get_spinorfieldsize(params), device),
	  clmem_phi_inv_eo(meta::get_eoprec_spinorfieldsize(params), device),
	  clmem_phi(meta::get_spinorfieldsize(params), device),
	  clmem_phi_mp(meta::get_spinorfieldsize(params), device),
	  clmem_phi_eo(meta::get_eoprec_spinorfieldsize(params), device),
	  clmem_phi_mp_eo(meta::get_eoprec_spinorfieldsize(params), device)
{
	fill_kernels();
}

Opencl_Module_Hmc::~Opencl_Module_Hmc()
{
	clear_kernels();
}
