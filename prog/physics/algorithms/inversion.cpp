/** @file
 * Implementation of the inversion algorithms
 */

#include "inversion.hpp"
#include "../meta/util.hpp"
#include "solver.hpp"
#include <cassert>
#include "../lattices/util.hpp"

static void invert_M_nf2_upperflavour(const physics::lattices::Spinorfield* result, const physics::lattices::Gaugefield& gaugefield, const physics::lattices::Spinorfield* source, const hardware::System& system);

void physics::algorithms::perform_inversion(const std::vector<const physics::lattices::Spinorfield*> * result, physics::lattices::Gaugefield* gaugefield, const std::vector<const physics::lattices::Spinorfield*>& sources, const hardware::System& system)
{
	int num_sources = sources.size();
	auto params = system.get_inputparameters();

	//apply stout smearing if wanted
	if(params.get_use_smearing())
		gaugefield->smear();

	for(int k = 0; k < num_sources; k++) {
		logger.debug() << "calling solver..";
		invert_M_nf2_upperflavour(result->at(k), *gaugefield, sources[k], system);
	}

	if(params.get_use_smearing())
		gaugefield->unsmear();
}

static void invert_M_nf2_upperflavour(const physics::lattices::Spinorfield* result, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield* source, const hardware::System& system)
{
	using namespace physics::lattices;
	using namespace physics::algorithms::solvers;
	using namespace physics::fermionmatrix;

	/** This solves the sparse-matrix system
	 *  A x = b
	 *  with  x == result
	 *        A == gaugefield
	 *        b == source
	 * using a Krylov-Solver (BiCGStab or CG)
	 */

	// assert a single GPU
	assert(result->get_buffers().size() == 1);

	auto params = system.get_inputparameters();

	auto result_buf = result->get_buffers().at(0);
	auto source_buf = source->get_buffers().at(0);
	auto gf_buf = gf.get_buffers().at(0);
	auto device = result_buf->get_device();

	int converged = -1;
	hardware::code::Fermions * solver = device->get_fermion_code();
	auto spinor_code = device->get_spinor_code();

	if(!params.get_use_eo()) {
		//noneo case
		//Trial solution
		///@todo this should go into a more general function
		result->cold();
		if(params.get_solver() == meta::Inputparameters::cg) {
			Spinorfield tmp(system);
			//to use cg, one needs an hermitian matrix, which is QplusQminus
			//the source must now be gamma5 b, to obtain the desired solution in the end
			copyData(&tmp, source);
			tmp.gamma5();
			QplusQminus f_neo(params.get_kappa(), meta::get_mubar(params), system);
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(result_buf, false, "\tinv. field before inversion ");
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(source_buf, false, "\tsource before inversion ");
			converged = cg(result, f_neo, gf, tmp, system, params.get_solver_prec());
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(result_buf, false, "\tinv. field after inversion ");
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(source_buf, false, "\tsource after inversion ");
			copyData(&tmp, result);
			//now, calc Qminus result_buf to obtain x = A^⁻1 b
			Qminus qminus(params.get_kappa(), meta::get_mubar(params), system);
			qminus(result, gf, tmp);
		} else {
			M f_neo(params.get_kappa(), meta::get_mubar(params), system);
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(result_buf, false, "\tinv. field before inversion ");
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(source_buf, false, "\tsource before inversion ");
			converged = bicgstab(result, f_neo, gf, *source, system, params.get_solver_prec());
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(result_buf, false, "\tinv. field after inversion ");
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(source_buf, false, "\tsource after inversion ");
		}
	} else {
		/**
		 * If even-odd-preconditioning is used, the inversion is split up
		 * into even and odd parts using Schur decomposition, assigning the
		 * non-trivial inversion to the even sites (see DeGran/DeTar p 174ff).
		 */
		//init some helping buffers
		const hardware::buffers::Spinor clmem_source_even  (meta::get_eoprec_spinorfieldsize(params), device);
		const hardware::buffers::Spinor clmem_source_odd  (meta::get_eoprec_spinorfieldsize(params), device);
		const hardware::buffers::Spinor clmem_tmp_eo_1  (meta::get_eoprec_spinorfieldsize(params), device);
		const hardware::buffers::Spinor clmem_tmp_eo_2  (meta::get_eoprec_spinorfieldsize(params), device);
		const hardware::buffers::Plain<hmc_complex> clmem_one (1, device);
		const hardware::buffers::Spinor result_buf_eo  (meta::get_eoprec_spinorfieldsize(params), device);
		const hardware::buffers::Plain<hmc_complex> clmem_mone (1, solver->get_device());
		hmc_complex one = hmc_complex_one;
		hmc_complex mone = {-1.,0.};
		clmem_one.load(&one);
		clmem_mone.load(&mone);

		//convert source and input-vector to eoprec-format
		/**
		 * This currently is a workaround connceted to issue #387
		 * the roles of even/odd vectors are interchanged!
		 * @todo: fix
		 * original code:
		 * spinor_code->convert_to_eoprec_device(&clmem_source_even, &clmem_source_odd, source_buf);
		 * workaround:
		 */
		spinor_code->convert_to_eoprec_device(&clmem_source_odd, &clmem_source_even, source_buf);

		if(logger.beDebug()) solver->print_info_inv_field(source_buf, false, "\tsource before inversion ");
		if(logger.beDebug()) solver->print_info_inv_field(&clmem_source_even, true, "\teven source before inversion ");
		if(logger.beDebug()) solver->print_info_inv_field(&clmem_source_odd, true, "\todd source before inversion ");

		//prepare sources
		/**
		 * This changes the even source according to (with A = M + D):
		 *  b_e = b_e - D_eo M_inv b_o
		 */
		if(params.get_fermact() == meta::Inputparameters::wilson) {
			//in this case, the diagonal matrix is just 1 and falls away.
			solver->dslash_eo_device(&clmem_source_odd, &clmem_tmp_eo_1, gf_buf, EVEN);
			spinor_code->saxpy_eoprec_device(&clmem_source_even, &clmem_tmp_eo_1, &clmem_one, &clmem_source_even);
		} else if(params.get_fermact() == meta::Inputparameters::twistedmass) {
			solver->M_tm_inverse_sitediagonal_device(&clmem_source_odd, &clmem_tmp_eo_1);
			solver->dslash_eo_device(&clmem_tmp_eo_1, &clmem_tmp_eo_2, gf_buf, EVEN);
			spinor_code->saxpy_eoprec_device(&clmem_source_even, &clmem_tmp_eo_2, &clmem_one, &clmem_source_even);
		}

		//Trial solution
		///@todo this should go into a more general function
		spinor_code->set_eoprec_spinorfield_cold_device(&result_buf_eo);
		logger.debug() << "start eoprec-inversion";
		//even solution
		if(params.get_solver() == meta::Inputparameters::cg) {
			//to use cg, one needs an hermitian matrix, which is QplusQminus
			//the source must now be gamma5 b, to obtain the desired solution in the end
			solver->gamma5_eo_device(&clmem_source_even);
			hardware::code::QplusQminus_eo f_eo(solver);
			if(logger.beDebug()) solver->print_info_inv_field(&result_buf_eo, true, "\tinv field before inversion ");
			if(logger.beDebug()) solver->print_info_inv_field(&clmem_source_even, true, "\tsource before inversion ");
			converged = solver->cg_eo(f_eo, &result_buf_eo, &clmem_source_even, gf_buf, params.get_solver_prec());
			if(logger.beDebug()) solver->print_info_inv_field(&result_buf_eo, true, "\tinv field after inversion ");
			if(logger.beDebug()) solver->print_info_inv_field(&clmem_source_even, true, "\tsource after inversion ");
			//now, calc Qminus result_buf_eo to obtain x = A^⁻1 b
			//therefore, use source as an intermediate buffer
			solver->Qminus_eo(&result_buf_eo, &clmem_source_even, gf_buf, params.get_kappa(), meta::get_mubar(params));
			//save the result to result_buf
			hardware::buffers::copyData(&result_buf_eo, &clmem_source_even);
		} else {
			hardware::code::Aee f_eo(solver);
			if(logger.beDebug()) solver->print_info_inv_field(&result_buf_eo, true, "\tinv field before inversion ");
			if(logger.beDebug()) solver->print_info_inv_field(&clmem_source_even, true, "\tsource before inversion ");
			converged = solver->bicgstab_eo(f_eo, &result_buf_eo, &clmem_source_even, gf_buf, params.get_solver_prec());
			if(logger.beDebug()) solver->print_info_inv_field(&result_buf_eo, true, "\tinv field after inversion ");
			if(logger.beDebug()) solver->print_info_inv_field(&clmem_source_even, true, "\tsource after inversion ");
		}

		//odd solution
		/** The odd solution is obtained from the even one according to:
		 *  x_o = M_inv D x_e - M_inv b_o
		 * @todo: find out why it must be
		 *  x_o = -(M_inv D x_e + M_inv b_o)
		 */
		if(params.get_fermact() == meta::Inputparameters::wilson) {
			//in this case, the diagonal matrix is just 1 and falls away.
			solver->dslash_eo_device(&result_buf_eo, &clmem_tmp_eo_1, gf_buf, ODD);
			spinor_code->saxpy_eoprec_device(&clmem_tmp_eo_1, &clmem_source_odd, &clmem_mone, &clmem_tmp_eo_1);
			spinor_code->sax_eoprec_device(&clmem_tmp_eo_1, &clmem_mone, &clmem_tmp_eo_1);
		} else if(params.get_fermact() == meta::Inputparameters::twistedmass) {
			solver->dslash_eo_device(&result_buf_eo, &clmem_tmp_eo_2, gf_buf, ODD);
			solver->M_tm_inverse_sitediagonal_device(&clmem_tmp_eo_2, &clmem_tmp_eo_1);
			solver->M_tm_inverse_sitediagonal_device(&clmem_source_odd, &clmem_tmp_eo_2);
			spinor_code->saxpy_eoprec_device(&clmem_tmp_eo_1, &clmem_tmp_eo_2, &clmem_mone, &clmem_tmp_eo_1);
			spinor_code->sax_eoprec_device(&clmem_tmp_eo_1, &clmem_mone, &clmem_tmp_eo_1);
		}

		 ///CP: whole solution
		//convert source and input-vector to eoprec-format
		/**
		 * This currently is a workaround connceted to issue #387
		 * the roles of even/odd vectors are interchanged!
		 * @todo: fix
		 * original code:
		 * //CP: suppose the even sol is saved in inout_eoprec, the odd one in clmem_tmp_eo_1
		 * spinor_code->convert_from_eoprec_device(&result_buf_eo, &clmem_tmp_eo_1, result_buf);
		 * workaround:
		 */
		//CP: suppose the odd sol is saved in inout_eoprec, the even one in clmem_tmp_eo_1
		spinor_code->convert_from_eoprec_device( &clmem_tmp_eo_1, &result_buf_eo, result_buf);
	}

	if(logger.beDebug()) solver->print_info_inv_field(result_buf, false, "\tsolution ");

	if (converged < 0) {
		if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
		else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
	} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
}

