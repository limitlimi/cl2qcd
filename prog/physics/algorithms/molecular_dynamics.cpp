/** @file
 * Implementation of the molecular dynamics algorithms
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * (c) 2012-2013 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
 */

#include "molecular_dynamics.hpp"

#include <stdexcept>
#include "../fermionmatrix/fermionmatrix.hpp"
#include "../../meta/util.hpp"
#include "solver.hpp"
#include "solver_shifted.hpp"
#include "forces.hpp"
#include "../../hardware/code/molecular_dynamics.hpp"
#include "../lattices/util.hpp"


void physics::algorithms::md_update_gaugefield(const physics::lattices::Gaugefield * const gf, const physics::lattices::Gaugemomenta& gm, const hmc_float eps)
{
	auto gf_bufs = gf->get_buffers();
	auto gm_bufs = gm.get_buffers();
	size_t num_bufs = gf_bufs.size();
	if(num_bufs != gm_bufs.size()) {
		throw std::invalid_argument("The gaugefield and the gaugemomenta use different devices.");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto gf_buf = gf_bufs[i];
		auto gm_buf = gm_bufs[i];
		auto code = gf_buf->get_device()->get_molecular_dynamics_code();
		code->md_update_gaugefield_device(gm_buf, gf_buf, eps);
	}

	logger.debug() << "\tHMC [UP]:\tupdate GF [" << eps << "]";
}

void physics::algorithms::md_update_spinorfield(const physics::lattices::Spinorfield * const out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& orig, const hardware::System& system, const hmc_float kappa, const hmc_float mubar)
{
	logger.debug() << "\tHMC [UP]:\tupdate SF";
	physics::fermionmatrix::Qplus qplus(kappa, mubar, system);
	qplus(out, gf, orig);
	log_squarenorm("Spinorfield after update", *out);
}

void physics::algorithms::md_update_spinorfield(const physics::lattices::Spinorfield_eo * const out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& orig, const hardware::System& system, const hmc_float kappa, const hmc_float mubar)
{
	logger.debug() << "\tHMC [UP]:\tupdate SF";
	physics::fermionmatrix::Qplus_eo qplus(kappa, mubar, system);
	qplus(out, gf, orig);
	log_squarenorm("Spinorfield after update", *out);
}

/**
 * This is the update of the Staggeredfield_eo in the RHMC. Actually, it is a tool to obtain
 * the new Staggeredfield_eo at the beginning of each iteration, given a gaussian drawn field.
 * In fact, this function makes the operator (Mdag*M) to some rational power act on the
 * Staggeredfield_eo "orig". To be clearer, here the following operation is implemented:
 * @code
 *   out = (a_0 + \sum_{i=1}^k a_i * (Mdag*M + b_i)^{-1} ) * orig 
 *       = a_0 * orig + \sum_{i=1}^k a_i * [(Mdag*M + b_i)^{-1} * orig]
 * @endcode
 * where k is the order of the rational approximation, a_0, a_i and b_i are the coefficients.
 */
void physics::algorithms::md_update_spinorfield(const physics::lattices::Staggeredfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Staggeredfield_eo& orig, physics::algorithms::Rational_Coefficients& coeff, const hardware::System& system, hmc_float mass)
{
	logger.debug() << "\tRHMC [UP]:\tupdate SF";
	auto params = system.get_inputparameters();
	const physics::fermionmatrix::MdagM_eo fm(system, mass);
	
	//Temporary fields for shifted inverter
	logger.trace() << "\t\tstart solver...";
	std::vector<physics::lattices::Staggeredfield_eo *> X;
	for(int i=0; i<coeff.Get_order(); i++)
		X.push_back(new physics::lattices::Staggeredfield_eo(system));
	//Here the inversion must be performed with high precision, because it'll be used for Metropolis test
	const int iterations = physics::algorithms::solvers::cg_m(X, coeff.Get_b(), fm, gf, orig, system, params.get_solver_prec());
	logger.trace() << "\t\t...end solver in " << iterations << " iterations";
	
	physics::lattices::sax(out, {coeff.Get_a0(), 0.}, orig);
	for(int i=0; i<coeff.Get_order(); i++)
		physics::lattices::saxpy(out, {(coeff.Get_a())[i], 0.}, orig, *out);
	
	log_squarenorm("Staggeredfield_eo after update", *out);
}

/**
 * template for md_update_spinorfield_mp
 * this needs 3 fermionmatrices in order to use cg as default solver (because for the cg one needs a hermitian matrix)
 */
template<class FERMIONMATRIX, class FERMIONMATRIX_CONJ, class FERMIONMATRIX_HERM, class SPINORFIELD> static void md_update_spinorfield_mp(const SPINORFIELD * const out, const physics::lattices::Gaugefield& gf, const SPINORFIELD& orig, const hardware::System& system, const hmc_float kappa, const hmc_float mubar)
{
        SPINORFIELD tmp(system);
	FERMIONMATRIX qplus(kappa, mubar, system);
	auto params = system.get_inputparameters();

	log_squarenorm("Spinorfield before update: ", orig);

	qplus(&tmp, gf, orig);

	/**
	 * Now one needs ( Qplus )^-1 (heavy_mass) using tmp as source to get phi_mp
	 */
	try {
	  /**
	   * Use BiCGStab as default here
	   * an exception will be thrown if the solver cannot solve
	   */
	  logger.debug() << "\t\t\tstart solver";
	  
	  /** 
	   * @todo at the moment, we can only put in a cold spinorfield
	   * or a point-source spinorfield as trial-solution
	   */
	  out->zero();
	  out->gamma5();
	  
	  FERMIONMATRIX qplus_mp(params.get_kappa_mp(), meta::get_mubar_mp(params), system);
	  const int iterations = physics::algorithms::solvers::bicgstab(out, qplus_mp, gf, tmp, system, params.get_solver_prec());
	} //try
	catch (physics::algorithms::solvers::SolverException& e ) {
	  logger.fatal() << e.what();
	  logger.info() << "Retry with CG...";
	  SPINORFIELD tmp2(system);
	  /** 
	   * @todo at the moment, we can only put in a cold spinorfield
	   * or a point-source spinorfield as trial-solution
	   */
	  tmp2.zero();
	  tmp2.gamma5();

	  FERMIONMATRIX_HERM fm_herm(params.get_kappa_mp(), meta::get_mubar_mp(params), system);
	  const int iterations_cg = physics::algorithms::solvers::cg(&tmp2, fm_herm, gf, tmp, system, params.get_solver_prec());
	  FERMIONMATRIX_CONJ fm_conj(params.get_kappa_mp(), meta::get_mubar_mp(params), system);
	  fm_conj(out, gf, tmp2);
	}
}

void physics::algorithms::md_update_spinorfield_mp(const physics::lattices::Spinorfield * const out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& orig, const hardware::System& system, const hmc_float kappa, const hmc_float mubar)
{
  logger.debug() << "\tHMC [UP]:\tupdate SF_MP";
	  using physics::fermionmatrix::Qplus;
	  using physics::fermionmatrix::Qminus;
	  using physics::fermionmatrix::QplusQminus;

	  ::md_update_spinorfield_mp<Qplus, Qminus, QplusQminus>(out, gf, orig, system, kappa, mubar);
}

void physics::algorithms::md_update_spinorfield_mp(const physics::lattices::Spinorfield_eo * const out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& orig, const hardware::System& system, const hmc_float kappa, const hmc_float mubar)
{
	logger.debug() << "\tHMC [UP]:\tupdate SF_MP";
	  using physics::fermionmatrix::Qplus_eo;
	  using physics::fermionmatrix::Qminus_eo;
	  using physics::fermionmatrix::QplusQminus_eo;

	  ::md_update_spinorfield_mp<Qplus_eo, Qminus_eo, QplusQminus_eo>(out, gf, orig, system, kappa, mubar);
}

template<class SPINORFIELD> static void md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const SPINORFIELD& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar)
{
	using namespace physics::algorithms;

	physics::lattices::Gaugemomenta delta_p(system);
	delta_p.zero();
	calc_total_force(&delta_p, gf, phi, system, kappa, mubar);

	logger.debug() << "\tHMC [UP]:\tupdate GM [" << eps << "]";
	physics::lattices::saxpy(inout, -1.*eps, delta_p);
}
void physics::algorithms::md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar)
{
	::md_update_gaugemomentum(inout, eps, gf, phi, system, kappa, mubar);
}
void physics::algorithms::md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar)
{
	::md_update_gaugemomentum(inout, eps, gf, phi, system, kappa, mubar);
}

void physics::algorithms::md_update_gaugemomentum_gauge(const physics::lattices::Gaugemomenta * const gm, const hmc_float eps, const physics::lattices::Gaugefield& gf, const hardware::System& system)
{
	const physics::lattices::Gaugemomenta force(system);
	force.zero();
	calc_gauge_force(&force, gf, system);
	log_squarenorm("\tHMC [UP]:\tFORCE [GAUGE]:\t", force);

	logger.debug() << "\tHMC [UP]:\tupdate GM [" << eps << "]";
	physics::lattices::saxpy(gm, -1.*eps, force);
}

template<class SPINORFIELD> static void md_update_gaugemomentum_fermion(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const SPINORFIELD& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar)
{
	using namespace physics::algorithms;

	const physics::lattices::Gaugemomenta force(system);
	force.zero();
	calc_fermion_force(&force, gf, phi, system, kappa, mubar);
	log_squarenorm("\tHMC [UP]:\tFORCE [DET]:\t", force);

	logger.debug() << "\tHMC [UP]:\tupdate GM [" << eps << "]";
	physics::lattices::saxpy(inout, -1.*eps, force);
}
void physics::algorithms::md_update_gaugemomentum_fermion(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar)
{
	::md_update_gaugemomentum_fermion(inout, eps, gf, phi, system, kappa, mubar);
}
void physics::algorithms::md_update_gaugemomentum_fermion(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar)
{
	::md_update_gaugemomentum_fermion(inout, eps, gf, phi, system, kappa, mubar);
}

template<class SPINORFIELD> static void md_update_gaugemomentum_detratio(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const SPINORFIELD& phi_mp, const hardware::System& system)
{
	using namespace physics::algorithms;

	const physics::lattices::Gaugemomenta force(system);
	force.zero();
	calc_detratio_forces(&force, gf, phi_mp, system);
	log_squarenorm("\tHMC [UP]:\tFORCE [DETRAT]:\t", force);

	logger.debug() << "\tHMC [UP]:\tupdate GM [" << eps << "]";
	physics::lattices::saxpy(inout, -1.*eps, force);
}
void physics::algorithms::md_update_gaugemomentum_detratio(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system)
{
	::md_update_gaugemomentum_detratio(inout, eps, gf, phi, system);
}
void physics::algorithms::md_update_gaugemomentum_detratio(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system)
{
	::md_update_gaugemomentum_detratio(inout, eps, gf, phi, system);
}
