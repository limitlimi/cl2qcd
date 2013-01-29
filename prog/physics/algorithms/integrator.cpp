/** @file
 * Implementation of the integrator algorithms
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * (c) 2012-2013 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
 */

#include "integrator.hpp"
#include "molecular_dynamics.hpp"
#include "../../meta/util.hpp"

template<class SPINORFIELD> static void leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);
template<class SPINORFIELD> static void leapfrog_1ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);
template<class SPINORFIELD> static void leapfrog_2ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);
template<class SPINORFIELD> static void leapfrog_3ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);

template<class SPINORFIELD> static void twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);
template<class SPINORFIELD> static void twomn_1ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);
template<class SPINORFIELD> static void twomn_2ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);
template<class SPINORFIELD> static void twomn_3ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);

void physics::algorithms::leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield& phi, const hardware::System& system)
{
	::leapfrog(gm, gf, phi, system);
}
void physics::algorithms::leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system)
{
	::leapfrog(gm, gf, phi, system);
}

template<class SPINORFIELD> void leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system)
{
	auto params = system.get_inputparameters();

	logger.trace() << "\tHMC [INT]:\tstart leapfrog...";
	//it is assumed that the new gaugefield and gaugemomentum have been set to the old ones already when this function is called the first time
	if(params.get_num_timescales() == 1) {
		leapfrog_1ts(gm, gf, phi, system);
	} else if (params.get_num_timescales() == 2) {
		leapfrog_2ts(gm, gf, phi, system);
	} else if (params.get_num_timescales() == 3) {
		leapfrog_3ts(gm, gf, phi, system);
	} else {
		Print_Error_Message("\tHMC [INT]:\tMore than 3 timescales is not implemented yet. Aborting...");
	}
	logger.trace() << "\tHMC [INT]:\t...finished leapfrog";
}

template<class SPINORFIELD> static void leapfrog_1ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system)
{
	using namespace physics::algorithms;

	auto params = system.get_inputparameters();
	const int n0 = params.get_integrationsteps(0);
	const hmc_float deltaTau0 = params.get_tau() / ((hmc_float) n0);
	const hmc_float deltaTau0_half = 0.5 * deltaTau0;

	md_update_gaugemomentum(gm, deltaTau0_half, *gf, phi, system);
	for(int k = 1; k < n0; k++) {
		md_update_gaugefield(gf, *gm, deltaTau0);
		md_update_gaugemomentum(gm, deltaTau0, *gf, phi, system);
	}
	md_update_gaugefield(gf, *gm, deltaTau0);
	md_update_gaugemomentum(gm, deltaTau0_half, *gf, phi, system);
}

template<class SPINORFIELD> static void leapfrog_2ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system)
{
	using namespace physics::algorithms;

	//this uses 2 timescales (more is not implemented yet): timescale0 for the gauge-part, timescale1 for the fermion part
	//this is done after hep-lat/0209037. See also hep-lat/0506011v2 for a more advanced version
	auto params = system.get_inputparameters();
	const int n0 = params.get_integrationsteps(0);
	const int n1 = params.get_integrationsteps(1);
	const hmc_float deltaTau1 = params.get_tau() / ((hmc_float) n1);
	const hmc_float deltaTau0 = deltaTau1 / ( (hmc_float) n0 );
	const hmc_float deltaTau0_half = 0.5 * deltaTau0;
	const hmc_float deltaTau1_half = 0.5 * deltaTau1;

	//this corresponds to V_s2(deltaTau/2)
	md_update_gaugemomentum_fermion(gm, deltaTau1_half, *gf, phi, system);
	//now, m steps "more" are performed for the gauge-part
	//this corresponds to [V_s1(deltaTau/2/m) V_t(deltaTau/m) V_s1(deltaTau/2/m) ]^m
	for(int l = 0; l < n0; l++) {
		if(l == 0) md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system);
		md_update_gaugefield(gf, *gm, deltaTau0);
		//one has to include the case of n1=1 here
		if(l == n0 - 1 && n1 == 1) md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system);
		else md_update_gaugemomentum_gauge(gm, deltaTau0, *gf, system);
	}
	for(int k = 1; k < n1; k++) {
		//this corresponds to V_s2(deltaTau)
		md_update_gaugemomentum_fermion(gm, deltaTau1, *gf, phi, system);
		for(int l = 0; l < n0; l++) {
			//this corresponds to [V_s1(deltaTau/2/m) V_t(deltaTau/m) V_s1(deltaTau/2/m) ]^m
			// where the first half_step has been carried out above already
			md_update_gaugefield(gf, *gm, deltaTau0);
			//md_update_gaugemomentum_gauge(deltaTau0_half);
			if(l == n0 - 1 && k == n1 - 1) md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system);
			else md_update_gaugemomentum_gauge(gm, deltaTau0, *gf, system);
		}
	}
	//this corresponds to the missing V_s2(deltaTau/2)
	md_update_gaugemomentum_fermion(gm, deltaTau1_half, *gf, phi, system);
}

template<class SPINORFIELD> static void leapfrog_3ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system)
{
	using namespace physics::algorithms;

	// just like with 2 timescales...
	auto params = system.get_inputparameters();
	const int n0 = params.get_integrationsteps(0);
	const int n1 = params.get_integrationsteps(1);
	const int n2 = params.get_integrationsteps(2);
	const hmc_float deltaTau2 = params.get_tau() / ((hmc_float) n2);
	const hmc_float deltaTau1 = deltaTau2 / ( (hmc_float) n1 );
	const hmc_float deltaTau0 = deltaTau1 / ( (hmc_float) n0 );
	const hmc_float deltaTau0_half = 0.5 * deltaTau0;
	const hmc_float deltaTau1_half = 0.5 * deltaTau1;
	const hmc_float deltaTau2_half = 0.5 * deltaTau2;

	//In this case one has to call the "normal" md_update_gaugemomentum_fermion with the heavier mass
	const hmc_float kappa_tmp = params.get_kappa_mp();
	const hmc_float mubar_tmp = meta::get_mubar_mp(params);

	md_update_gaugemomentum_detratio(gm, deltaTau2_half, *gf, phi, system);
	//now, n1 steps "more" are performed for the fermion-part
	for(int l = 0; l < n1; l++) {
		if(l == 0) md_update_gaugemomentum_fermion(gm, deltaTau1_half, *gf, phi, system, kappa_tmp, mubar_tmp);
		//now, n0 steps "more" are performed for the gauge-part
		for(int j = 0; j < n0; j++) {
			if(l == 0 && j == 0) md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system);
			md_update_gaugefield(gf, *gm, deltaTau0);
			if(j == n0 - 1 && l == n1 - 1 && n2 == 1) md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system);
			else md_update_gaugemomentum_gauge(gm, deltaTau0, *gf, system);
		}
		if(l == n1 - 1 && n2 == 1) md_update_gaugemomentum_fermion(gm, deltaTau1_half, *gf, phi, system, kappa_tmp, mubar_tmp);
		else md_update_gaugemomentum_fermion(gm, deltaTau1, *gf, phi, system, kappa_tmp, mubar_tmp);
	}
	//perform n2 - 1 intermediate steps
	for(int k = 1; k < n2; k++) {
		md_update_gaugemomentum_detratio(gm, deltaTau2, *gf, phi, system);
		for(int l = 0; l < n1; l++) {
			for(int j = 0; j < n0; j++) {
				md_update_gaugefield(gf, *gm, deltaTau0);
				if(j == n0 - 1 && l == n1 - 1 &&  k == n2 - 1) md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system);
				else md_update_gaugemomentum_gauge(gm, deltaTau0, *gf, system);
			}
			if(l == n1 - 1 && k == n2 - 1) md_update_gaugemomentum_fermion(gm, deltaTau1_half, *gf, phi, system, kappa_tmp, mubar_tmp);
			else md_update_gaugemomentum_fermion(gm, deltaTau1, *gf, phi, system, kappa_tmp, mubar_tmp);
		}
	}
	md_update_gaugemomentum_detratio(gm, deltaTau2_half, *gf, phi, system);
}

void physics::algorithms::twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield& phi, const hardware::System& system)
{
	::twomn(gm, gf, phi, system);
}
void physics::algorithms::twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system)
{
	::twomn(gm, gf, phi, system);
}

template<class SPINORFIELD> void twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system)
{
	auto params = system.get_inputparameters();

	logger.trace() << "\tHMC [INT]\tstarting 2MN...";
	//it is assumed that the new gaugefield and gaugemomentum have been set to the old ones already when this function is called the first time
	if(params.get_num_timescales() == 1) {
		twomn_1ts(gm, gf, phi, system);
	} else if (params.get_num_timescales() == 2) {
		twomn_2ts(gm, gf, phi, system);
	} else if (params.get_num_timescales() == 3) {
		twomn_3ts(gm, gf, phi, system);
	} else {
		Print_Error_Message("\tHMC [INT]:\tMore than 3 timescales is not implemented yet. Aborting...");
	}
	logger.debug() << "\tHMC [INT]:\tfinished 2MN";
}

template<class SPINORFIELD> void twomn_1ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system)
{
	using namespace physics::algorithms;

	auto params = system.get_inputparameters();
	const int n0 = params.get_integrationsteps(0);
	const hmc_float deltaTau0 = params.get_tau() / ((hmc_float) n0);
	const hmc_float deltaTau0_half = 0.5 * deltaTau0;
	const hmc_float lambda_times_deltaTau0 = deltaTau0 * params.get_lambda(0);
	const hmc_float one_minus_2_lambda = 1. - 2.*params.get_lambda(0);
	const hmc_float one_minus_2_lambda_times_deltaTau0 = one_minus_2_lambda * deltaTau0;

	md_update_gaugemomentum(gm, lambda_times_deltaTau0, *gf, phi, system);
	md_update_gaugefield(gf, *gm, deltaTau0_half);
	md_update_gaugemomentum(gm, one_minus_2_lambda_times_deltaTau0, *gf, phi, system);
	md_update_gaugefield(gf, *gm, deltaTau0_half);

	for(int k = 1; k < n0; k++) {
		md_update_gaugemomentum(gm, 2.*lambda_times_deltaTau0, *gf, phi, system);
		md_update_gaugefield(gf, *gm, deltaTau0_half);
		md_update_gaugemomentum(gm, one_minus_2_lambda_times_deltaTau0, *gf, phi, system);
		md_update_gaugefield(gf, *gm, deltaTau0_half);
	}
	md_update_gaugemomentum(gm, lambda_times_deltaTau0, *gf, phi, system);
}

template<class SPINORFIELD> void twomn_2ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system)
{
	using namespace physics::algorithms;

	//this is done after hep-lat/0209037. See also hep-lat/0506011v2 for a more advanced version
	auto params = system.get_inputparameters();
	const int n0 = params.get_integrationsteps(0);
	const int n1 = params.get_integrationsteps(1);

	//this uses 2 timescales (more is not implemented yet): timescale1 for the gauge-part, timescale2 for the fermion part
	const hmc_float deltaTau1 = params.get_tau() / ((hmc_float) n1);
	//NOTE: With 2MN, the stepsize for the lower integration step is deltaTau1/(2 N0)!!
	const hmc_float deltaTau0 = deltaTau1 / ( 2.* (hmc_float) n0 );
	const hmc_float deltaTau0_half = 0.5 * deltaTau0;
	//hmc_float deltaTau1_half = 0.5 * deltaTau1;
	const hmc_float lambda0_times_deltaTau0 = deltaTau0 * params.get_lambda(0);
	const hmc_float lambda1_times_deltaTau1 = deltaTau1 * params.get_lambda(1);
	const hmc_float one_minus_2_lambda0 = 1. - 2.*params.get_lambda(0);
	const hmc_float one_minus_2_lambda1 = 1. - 2.*params.get_lambda(1);
	const hmc_float one_minus_2_lambda0_times_deltaTau0 = one_minus_2_lambda0 * deltaTau0;
	const hmc_float one_minus_2_lambda1_times_deltaTau1 = one_minus_2_lambda1 * deltaTau1;

	md_update_gaugemomentum_fermion(gm, lambda1_times_deltaTau1, *gf, phi, system);
	//now, n0 steps "more" are performed for the gauge-part
	//this corresponds to [exp(lambda*eps T(V_gauge) ) exp( eps/2 V ) exp( (1 - 2lamdba) *eps T(V_gauge) ) exp( eps/2 V ) exp( lamdba*eps T(V_ga    uge) ) ]^m
	for(int l = 0; l < n0; l++) {
		if(l == 0) md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system);
		md_update_gaugefield(gf, *gm, deltaTau0_half);
		md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
		md_update_gaugefield(gf, *gm, deltaTau0_half);
		md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
	}
	//this corresponds to V_s2( ( 1 - 2lambda) *deltaTau)
	md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system);
	//now, m steps "more" are performed for the gauge-part (again)
	for(int l = 0; l < n0; l++) {
		//the first half step has been carried out above already
		md_update_gaugefield(gf, *gm, deltaTau0_half);
		md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
		md_update_gaugefield(gf, *gm, deltaTau0_half);
		//in case one does not perform intermediate steps after this, one must perform a half_step only!
		if(l == n0 - 1 && n1 == 1) md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system);
		else md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
	}
	//the last V_s2(lambda*deltaTau) can be pulled into the intermediate steps
	for(int k = 1; k < n1; k++) {
		//this corresponds to V_s2(deltaTau)
		md_update_gaugemomentum_fermion(gm, 2.*lambda1_times_deltaTau1, *gf, phi, system);
		for(int l = 0; l < n0; l++) {
			//the first half step has been carried out above already
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
		}
		md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system);
		for(int l = 0; l < n0; l++) {
			//the first half step has been carried out above already
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			if(l == n0 - 1 && k == n1 - 1) md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system);
			else md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
		}
	}
	//this corresponds to the missing V_s2(lambda*deltaTau)
	md_update_gaugemomentum_fermion(gm, lambda1_times_deltaTau1, *gf, phi, system);
}

template<class SPINORFIELD> void twomn_3ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system)
{
	using namespace physics::algorithms;

	//just like with 2 timescales...
	auto params = system.get_inputparameters();
	const int n0 = params.get_integrationsteps(0);
	const int n1 = params.get_integrationsteps(1);
	const int n2 = params.get_integrationsteps(2);

	const hmc_float deltaTau2 = params.get_tau() / ((hmc_float) n2);
	//NOTE: With 2MN, the stepsize for the lower integration step is deltaTau1/(2 Ni-1)!!
	const hmc_float deltaTau1 = deltaTau2 / ( 2.* (hmc_float) n1 );
	const hmc_float deltaTau0 = deltaTau1 / ( 2.* (hmc_float) n0 );
	const hmc_float deltaTau1_half = 0.5 * deltaTau1;
	const hmc_float deltaTau0_half = 0.5 * deltaTau0;
	//const hmc_float deltaTau1_half = 0.5 * deltaTau1;
	const hmc_float lambda0_times_deltaTau0 = deltaTau0 * params.get_lambda(0);
	const hmc_float lambda1_times_deltaTau1 = deltaTau1 * params.get_lambda(1);
	const hmc_float lambda2_times_deltaTau2 = deltaTau2 * params.get_lambda(2);
	const hmc_float one_minus_2_lambda0 = 1. - 2.*params.get_lambda(0);
	const hmc_float one_minus_2_lambda1 = 1. - 2.*params.get_lambda(1);
	const hmc_float one_minus_2_lambda2 = 1. - 2.*params.get_lambda(2);
	const hmc_float one_minus_2_lambda0_times_deltaTau0 = one_minus_2_lambda0 * deltaTau0;
	const hmc_float one_minus_2_lambda1_times_deltaTau1 = one_minus_2_lambda1 * deltaTau1;
	const hmc_float one_minus_2_lambda2_times_deltaTau2 = one_minus_2_lambda2 * deltaTau2;

	//In this case one has to call the "normal" md_update_gaugemomentum_fermion with the heavier mass
	const hmc_float kappa = params.get_kappa_mp();
	const hmc_float mubar = meta::get_mubar_mp(params);

	md_update_gaugemomentum_detratio(gm, lambda2_times_deltaTau2, *gf, phi, system);
	for(int l = 0; l < n1; l++) {
		if(l == 0) md_update_gaugemomentum_fermion(gm, lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
		for(int j = 0; j < n0; j++) {
			if(j == 0 && l == 0) md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system);
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
		}
		md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
		for(int j = 0; j < n0; j++) {
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
		}
		md_update_gaugemomentum_fermion(gm, 2.*lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
	}
	md_update_gaugemomentum_detratio(gm, one_minus_2_lambda2_times_deltaTau2, *gf, phi, system);
	for(int l = 0; l < n1; l++) {
		for(int j = 0; j < n0; j++) {
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
		}
		md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
		for(int j = 0; j < n0; j++) {
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			if(j == n0 - 1 && l == n1 - 1 && n2 == 1) md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system);
			else md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
		}
		if(l == n1 - 1 && n2 == 1) md_update_gaugemomentum_fermion(gm, lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
		else md_update_gaugemomentum_fermion(gm, 2.*lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
	}
	for(int k = 1; k < n2; k++) {
		//this corresponds to V_s2(deltaTau)
		md_update_gaugemomentum_detratio(gm, 2.*lambda2_times_deltaTau2, *gf, phi, system);
		for(int l = 0; l < n1; l++) {
			for(int j = 0; j < n0; j++) {
				md_update_gaugefield(gf, *gm, deltaTau0_half);
				md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
				md_update_gaugefield(gf, *gm, deltaTau0_half);
				md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
			}
			md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
			for(int j = 0; j < n0; j++) {
				md_update_gaugefield(gf, *gm, deltaTau0_half);
				md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
				md_update_gaugefield(gf, *gm, deltaTau0_half);
				md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
			}
			md_update_gaugemomentum_fermion(gm, 2.*lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
		}
		md_update_gaugemomentum_detratio(gm, one_minus_2_lambda2_times_deltaTau2, *gf, phi, system);
		for(int l = 0; l < n1; l++) {
			for(int j = 0; j < n0; j++) {
				md_update_gaugefield(gf, *gm, deltaTau0_half);
				md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
				md_update_gaugefield(gf, *gm, deltaTau0_half);
				md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
			}
			md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
			for(int j = 0; j < n0; j++) {
				md_update_gaugefield(gf, *gm, deltaTau0_half);
				md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
				md_update_gaugefield(gf, *gm, deltaTau0_half);
				if(j == n0 - 1 && l == n1 - 1 && k == n2 - 1) md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system);
				else md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
			}
			if(l == n1 - 1 && k == n2 - 1) md_update_gaugemomentum_fermion(gm, lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
			else md_update_gaugemomentum_fermion(gm, 2.*lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
		}
	}
	md_update_gaugemomentum_detratio(gm, lambda2_times_deltaTau2, *gf, phi, system);
}
