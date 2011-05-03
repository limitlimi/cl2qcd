/** @file
 * Input file handling
 */

#ifndef _INPUTH_
#define _INPUTH_

#include "hmcerrs.h"
#include "types.h"
#include "globaldefs.h"

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <cstring>
using namespace std;

/**
 * Parser for the input file.
 */
class inputparameters {
public:
	/**
	 * Default constructor loads default values.
	 */
	inputparameters() {
		set_defaults();
	};
	/**
	 * Parse the given file, overwriting current
	 * parameter values.
	 *
	 * \param ifn Name of the input file to parse
	 * \return Error as defined in hmcerrs.h:
	 *         \li HMC_FILEERROR if file cannot be opened
	 *         \li HMC_SUCCESS otherwise
	 */
	hmc_error readfile(char* ifn);
	/**
	 * Reset all parameters to default values.
	 *
	 * \return Should always return HMC_SUCCESS.
	 */
	hmc_error set_defaults();
	hmc_float get_kappa();
	hmc_float get_beta();
	hmc_float get_tau_fermion();
	hmc_float get_tau_gauge();
	hmc_float get_theta_fermion();
	hmc_float get_theta_gaugefield();
	hmc_float get_mu();
	hmc_float get_csw();
	hmc_float get_chem_pot_re();
	hmc_float get_chem_pot_im();
	int get_cgmax();
	int get_prec();
	int get_startcondition();
	int get_thermalizationsteps();
	int get_heatbathsteps();
	int get_overrelaxsteps();
	int get_hmcsteps();
	int get_saveconfigs();
	int get_savefrequency();
	int get_writefrequency();
	void display_sourcefile();
	void display_sourcefilenumber();
	//CP
	//this is out of laziness
	std::string sourcefile;
	std::string sourcefilenumber;
private:
	hmc_float kappa;
	hmc_float beta;
	hmc_float mu;
	hmc_float csw;
	hmc_float theta_fermion;
	hmc_float theta_gaugefield;
	hmc_float chem_pot_re;
	hmc_float chem_pot_im;
	hmc_float tau_fermion;
	hmc_float tau_gauge;
	int cgmax;
	int prec;
	int startcondition;
	int thermalizationsteps;
	int heatbathsteps;
	int overrelaxsteps;
	int hmcsteps;
	int savefrequency;
	int saveconfigs;
	int writefrequency;
	void val_assign(hmc_float* out, std::string line);
	void val_assign(int * out, std::string line);
	void sourcefilenumber_assign(std::string * out);
	void cond_assign(int * out, std::string line);
	void savecond_assign(int * out, std::string line);
	void val_assign(std::string * out, std::string line);
};

#endif /* _INPUTH_ */
