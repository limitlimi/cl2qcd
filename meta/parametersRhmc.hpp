/** @file
 *
 * Copyright (c) 2014 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
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

#ifndef _META_PARAMETERS_RHMC_HPP_
#define _META_PARAMETERS_RHMC_HPP_

#include <vector>
#include <string>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

namespace meta{
class ParametersRhmc
{
public:
	int get_md_approx_ord() const noexcept;
	int get_metro_approx_ord() const noexcept;
	int get_findminmax_iteration_block_size() const noexcept;
	int get_findminmax_max() const noexcept;
	double get_findminmax_prec() const noexcept;
	bool get_conservative() const noexcept;
	int get_num_tastes() const noexcept;
	double get_approx_lower() const noexcept;
	double get_approx_upper() const noexcept;
	int get_rhmcsteps() const noexcept;
	std::string get_approx_heatbath_file() const noexcept;
	std::string get_approx_md_file() const noexcept;
	std::string get_approx_metropolis_file() const noexcept;
	bool get_read_rational_approximations_from_file() const noexcept;


protected:
	/** @TODO If the rational approximation is read from file than its parameters could differ
	 *        from the following! This means, for example, that one could use get_md_approx_ord()
	 *        to get a value that is not that loaded from the file!
	 *  @TODO If read_rational_approximations_from_file is false it makes no sense to have the
	 *        approx_*_file variables, but this is similar to the gauge configuration.
	 */
	int md_approx_ord;
	int metro_approx_ord;
	int findminmax_iteration_block_size;
	int findminmax_max;
	double findminmax_prec;
	bool conservative; //this is for the strategy in findminmax_eigenvalues
	int num_tastes; //the numerator of the power of the determinant and then of the Rational Approx.
	double approx_lower;
	double approx_upper; //range of validity of the Rational Approximation
	int rhmcsteps;
	bool read_rational_approximations_from_file;
	std::string approx_heatbath_file;
	std::string approx_md_file;
	std::string approx_metropolis_file;

	po::options_description getOptions();
};

}

#endif