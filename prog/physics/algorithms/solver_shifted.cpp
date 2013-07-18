/** @file
 * Implementation of the shifted solver algorithms
 *
 * (c) 2013 Alessandro Sciarra <sciarra@th.uni-frankfurt.de>
 */

#include "solver.hpp"
#include "solver_shifted.hpp"

#include "../../logger.hpp"
//#include "../../operations_complex.h"
//#include "../../meta/type_ops.hpp"
//#include "../lattices/util.hpp"
#include "../lattices/scalar_complex.hpp"
//#include <cmath>
#include <sstream>
#include <vector>

static std::string create_log_prefix_cgm(int number) noexcept;
void log_squarenorm_aux(const std::string& msg, const std::vector<physics::lattices::Staggeredfield_eo *> x, const int n) noexcept;


int physics::algorithms::solvers::cg_m(const std::vector<physics::lattices::Staggeredfield_eo *> x, const std::vector<hmc_float> sigma, const physics::fermionmatrix::Fermionmatrix_stagg_eo& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Staggeredfield_eo& b, const hardware::System& system, hmc_float prec)
{
	if(sigma.size() != x.size())
		throw std::invalid_argument("Wrong size of multi-shifted inverter parameters!");
	
	using namespace physics::lattices;
	using physics::algorithms::solvers::SolverStuck;
	using physics::algorithms::solvers::SolverDidNotSolve;

	auto params = system.get_inputparameters();
	
	/// @todo start timer synchronized with device(s)
	klepsydra::Monotonic timer;
	/// @todo make configurable from outside
	const int RESID_CHECK_FREQUENCY = params.get_cg_iteration_block_size();
	const bool USE_ASYNC_COPY = params.get_cg_use_async_copy();
	if(USE_ASYNC_COPY) {
		logger.warn() << "Asynchroneous copying in the CG-M is currently unimplemented!";
	}

	//The number of values of constants sigma as well as the number of fields x
	const int Neqs = sigma.size();
	
	//Auxiliary staggered fields
	const Staggeredfield_eo* r = new Staggeredfield_eo(system);
	const Staggeredfield_eo* p = new Staggeredfield_eo(system);
	std::vector<const Staggeredfield_eo*> ps;
	
	//Auxiliary scalar vectors
	std::vector<const Scalar<hmc_complex>*> alpha;
	std::vector<const Scalar<hmc_complex>*> beta;
	std::vector<const Scalar<hmc_complex>*> zeta_i;   //This is zeta at the step iter-1
	std::vector<const Scalar<hmc_complex>*> zeta_ii;  //This is zeta at the step iter
	std::vector<const Scalar<hmc_complex>*> zeta_iii; //This is zeta at the step iter+1
	//Auxiliary scalars
	const Scalar<hmc_complex> alpha_scalar(system);
	const Scalar<hmc_complex> beta_scalar_prev(system);    //This is beta_scalar at the step iter-1
	const Scalar<hmc_complex> beta_scalar(system);         //This is beta_scalar at the step iter
	
	//Auxiliary containers for temporary saving
	const Staggeredfield_eo* v = new Staggeredfield_eo(system); //this is to store A.p
	const Scalar<hmc_complex> tmp1(system);                     //this is to store (r,r) before updating r
	const Scalar<hmc_complex> tmp2(system);                     //this is to store (r,r) after updating r
	const Scalar<hmc_complex> tmp3(system);                     //this is to store (p,v) as Scalar
	const Scalar<hmc_complex> num(system);                      //this is to store constants numerators
	const Scalar<hmc_complex> den(system);                      //this is to store constants denumerators
	//Auxiliary constants as Scalar
	const Scalar<hmc_complex> minus_one(system);
	minus_one.store(hmc_complex_minusone);
	const Scalar<hmc_complex> one(system);
	one.store(hmc_complex_one);
	
	
	hmc_float resid;
	int iter = 0;
	
	//Initialization auxilary and output quantities
	for(int i=0; i<Neqs; i++){
		x[i]->set_zero();                                     // x[i] = 0 
		ps.push_back(new Staggeredfield_eo(system));    
		copyData(ps[i], b);                                    // ps[i] = b
		alpha.push_back(new Scalar<hmc_complex>(system));
		alpha[i]->store(hmc_complex_zero);                    // alpha[i] = 0
		zeta_i.push_back(new Scalar<hmc_complex>(system));
		zeta_ii.push_back(new Scalar<hmc_complex>(system));
		zeta_iii.push_back(new Scalar<hmc_complex>(system));
		zeta_i[i]->store(hmc_complex_one);                    // zeta_i[i] = 1
		zeta_ii[i]->store(hmc_complex_one);                   // zeta_ii[i] = 1
	}
	copyData(r, b);                          // r = b
	copyData(p, b);                          // p = b
	scalar_product(&tmp1, *r, *r);           // set tmp1 = (r, r) for the first iteration
	beta_scalar.store(hmc_complex_one);      // beta_scalar = 1, here I should set beta_scalar_prev
						  // but in this way I can set beta_scalar_prev at the begin
						  // of the loop over iter recursively.
	
	//At first, in the log string I will report only the data about the first 3 system.
	//This is to avoid to print to shell tens of lines per time, since Neqs can be 
	//also 20 or something like that. In case, modify the int passed to log_squarenorm
	const int report_num = 3;
	log_squarenorm(create_log_prefix_cgm(iter) + "b (initial): ", b);
	log_squarenorm_aux(create_log_prefix_cgm(iter) + "x (initial)", x, report_num);
	log_squarenorm(create_log_prefix_cgm(iter) + "r (initial): ", *r);
	
	//NOTE: Here, most of the complex numbers may also be just hmc_floats.
	//      However, for this one would need some additional functions...
	for(iter = 0; iter < params.get_cgmax(); iter ++) {
		//Update beta_scalar: v=A.p and tmp1=(r,r) and tmp3=(p,v) ---> beta_scalar=(-1)*tmp1/tmp3
		copyData(&beta_scalar_prev,beta_scalar);  //before updating beta_scalar its value is saved
		A(v, gf, *p);
		log_squarenorm(create_log_prefix_cgm(iter) + "v: ", *v);
		scalar_product(&tmp3, *p, *v); 
		divide(&beta_scalar, tmp1, tmp3);  //tmp1 is set from previous iteration
		multiply(&beta_scalar, minus_one, beta_scalar);
		//Update field r: r+=beta_scalar*A.p ---> r = r + beta_scalar*v
		saxpy(r, beta_scalar, *v, *r);
		log_squarenorm(create_log_prefix_cgm(iter) + "r: ", *r);
		
		//We store in tmp2 the quantity (r,r) that we use later. When we check
		//the residuum, then it is already calculated.
		scalar_product(&tmp2, *r, *r);
		//Check whether the algorithm converged
		if(iter % RESID_CHECK_FREQUENCY == 0) {
			//Calculate resid: 
			resid = tmp2.get().re;
			logger.debug() << create_log_prefix_cgm(iter) << "resid: " << resid;
			//test if resid is NAN
			if(resid != resid) {
				logger.fatal() << create_log_prefix_cgm(iter) << "NAN occured!";
				throw SolverStuck(iter, __FILE__, __LINE__);
			}
			if(resid < prec) {
				logger.debug() << create_log_prefix_cgm(iter) << "Solver converged in " << iter << " iterations! resid:\t" << resid;
				// report on performance
				if(logger.beInfo()) {
					
					//TODO
					
				}
				// report on solution
				//log_squarenorm_aux(create_log_prefix_cgm(iter) + "x (final): ", x, report_num);
				return iter;
			}
		}
		//Update alpha_scalar: alpha_scalar = tmp2/tmp1
		divide(&alpha_scalar, tmp2, tmp1);
		//Update field p: p = r + alpha_scalar*p
		saxpy(p, alpha_scalar, *p, *r);
		//Loop over the system equations, namely over the set of sigma values
		for(int k=0; k<Neqs; k++){
			//Update zeta_iii[k]: num = zeta_i[k]*zeta_ii[k]*beta_scalar_prev
			//                    den = beta_scalar*alpha_scalar*(zeta_i[k]-zeta_ii[k])+
			//                        + zeta_i[k]*beta_scalar_prev*(1-sigma[k]*beta_scalar)
			// ---> zeta_iii[k] = num/den
			//TODO
			//Update beta[k]: beta[k] = beta_scalar*zeta_iii[k]/zeta_ii[k]
			//TODO
			//Update x[k]: x[k] = x[k] - beta[k]*ps[k] ---> use num to store (- beta[k]) for saxpy
			//TODO
			//Update alpha[k]: num = alpha_scalar*zeta_iii[k]*beta[k]
			//                 den = zeta_ii[k]*beta_scalar
			// ---> alpha[k] = num/den
			//TODO
			//Update ps[k]: ps[k] = zeta_iii[k]*r + alpha[k]*ps[k]
			//TODO
		}
		
		
		//Set tmp1 from tmp2 for following iteration
		copyData(&tmp1, tmp2);
	}
	
	logger.fatal() << create_log_prefix_cgm(iter) << "Solver did not solve in " << params.get_cgmax() << " iterations. Last resid: " << resid;
	throw SolverDidNotSolve(iter, __FILE__, __LINE__);
}



static std::string create_log_prefix_solver(std::string name, int number) noexcept
{
  using namespace std;
  string separator_big = "\t";
  string separator_small = " ";
  string label = "SOLVER";

  stringstream strnumber;
  strnumber.fill('0');
  /// @todo this should be length(cgmax)
  strnumber.width(6);
  strnumber << right << number;
  stringstream outfilename;
  outfilename << separator_big << label << separator_small << "[" << name << "]" << separator_small << "[" << strnumber.str() << "]:" << separator_big;
  string outputfile = outfilename.str();
  return outputfile;
}

static std::string create_log_prefix_cgm(int number) noexcept
{
  return create_log_prefix_solver("CG-M", number);
}


void log_squarenorm_aux(const std::string& msg, const std::vector<physics::lattices::Staggeredfield_eo *> x, const int n) noexcept
{
	if(logger.beDebug()) {
		hmc_float tmp;
		for(int i=0; i<n; i++){
			tmp = squarenorm(*x[i]);
			logger.debug() << msg << "[field_" << i << "]:" << std::scientific << std::setprecision(10) << tmp;
		}
	}
}



