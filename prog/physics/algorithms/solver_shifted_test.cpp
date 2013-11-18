/** @file
 * Tests of the multi-shifted inverter algorithm
 * 
 * (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "solver_shifted.hpp"
#include "rational_approximation.hpp"
#include "../../logger.hpp"
#include "../lattices/util.hpp"
#include "../lattices/scalar_complex.hpp"
#include "../lattices/staggeredfield_eo.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::algorithms::solvers::solver_shifted
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(cgm_1)
{
	using namespace physics::lattices;
	using namespace physics::algorithms::solvers;
	
	const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg"};
	meta::Inputparameters params(3, _params);
	hardware::System system(params);
	physics::PRNG prng(system);
	
	//This are some possible values of sigma
	hmc_float pol[5] = {0.0002065381736724, 0.00302707751065980, 0.0200732678058145,
	                                        0.12517586269872370, 1.0029328743375700};
	std::vector<hmc_float> sigma(pol, pol + sizeof(pol)/sizeof(hmc_float));
	physics::fermionmatrix::MdagM_eo matrix(system, 0.567);
	
	//This configuration for the Ref.Code is the same as for example dks_input_5
	Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
	Staggeredfield_eo b(system);
	std::vector<Staggeredfield_eo*> out;
	for(uint i=0; i<sigma.size(); i++)
		out.push_back(new Staggeredfield_eo(system));
	//This field is NOT that of the test explicit_stagg (D_KS_eo) because here the lattice is 4^4
	pseudo_randomize<Staggeredfield_eo, su3vec>(&b, 13);
  
	int iter = cg_m(out, sigma, matrix, gf, b, system, 1.e-23);
	logger.info() << "CG-M algorithm converged in " << iter << " iterations.";
	
	//Once obtained the solution we apply each operator (matrix + sigma) onto the field
	//out and we should obtain the field b. Hence we compare the squarenorms of b and of 
	//(matrix + sigma) * out.
	std::vector<Staggeredfield_eo*> aux;
	std::vector<hmc_float> sqnorm_out;
	hmc_float sqnorm_b = squarenorm(b);
	logger.info() << "                           sqnorm(b)=" << std::setprecision(16) << sqnorm_b;
	for(uint i=0; i<sigma.size(); i++){
		aux.push_back(new Staggeredfield_eo(system));
		matrix(aux[i], gf, *out[i]);
		saxpy(aux[i], {sigma[i],0.}, *out[i], *aux[i]);
		sqnorm_out.push_back(squarenorm(*aux[i]));
		logger.info() << "sqnorm((matrix + sigma[" << i << "]) * out[" << i << "])=" << std::setprecision(16) << sqnorm_out[i];
		BOOST_CHECK_CLOSE(sqnorm_b, sqnorm_out[i], 1.e-8);
	}
}

BOOST_AUTO_TEST_CASE(cgm_2)
{
	using namespace physics::lattices;
	using namespace physics::algorithms::solvers;
	using namespace physics::algorithms;
	
	const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg"};
	meta::Inputparameters params(3, _params);
	hardware::System system(params);
	physics::PRNG prng(system);
	
	//These are some possible values of sigma
	Rational_Approximation approx(8, 1,2, 1.e-5,1);
	std::vector<hmc_float> sigma = approx.Get_b();
	physics::fermionmatrix::MdagM_eo matrix(system, 1.01335);
	
	//This configuration for the Ref.Code is the same as for example dks_input_5
	Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
	Staggeredfield_eo b(system);
	std::vector<Staggeredfield_eo*> out;
	for(uint i=0; i<sigma.size(); i++)
		out.push_back(new Staggeredfield_eo(system));
	//This field is that of the test explicit_stagg, part 2 (D_KS_eo)
	pseudo_randomize<Staggeredfield_eo, su3vec>(&b, 123);
	
	//These are the sqnorms of the output of the CG-M algorithm from the reference code
	std::vector<hmc_float> sqnorms_ref;
	sqnorms_ref.push_back(59.877118728264179026);
	sqnorms_ref.push_back(59.875307563978779513);
	sqnorms_ref.push_back(59.865246895125224569);
	sqnorms_ref.push_back(59.811193256571264953);
	sqnorms_ref.push_back(59.521693061471694364);
	sqnorms_ref.push_back(57.988143162101970063);
	sqnorms_ref.push_back(50.327004662008008040);
	sqnorms_ref.push_back(22.236536652925686042);
	//Now I calculate the fields out
	int iter = cg_m(out, sigma, matrix, gf, b, system, 1.e-12);
	logger.info() << "CG-M algorithm converged in " << iter << " iterations.";
	
	std::vector<hmc_float> sqnorm_out;
	for(uint i=0; i<sigma.size(); i++){
		sqnorm_out.push_back(squarenorm(*out[i]));
		logger.info() << "sqnorm(out[" << i << "])=" << std::setprecision(16) << sqnorm_out[i];
		BOOST_CHECK_CLOSE(sqnorms_ref[i], sqnorm_out[i], 1.e-8);
	}
}


