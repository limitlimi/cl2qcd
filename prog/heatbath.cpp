#include "general_header.h"

#include "meta/util.hpp"
#include "physics/prng.hpp"
#include "physics/lattices/gaugefield.hpp"
#include "physics/algorithms/heatbath.hpp"

int main(int argc, const char* argv[])
{
	try {
		using namespace physics;
		using namespace physics::lattices;
		using physics::algorithms::heatbath;

		meta::Inputparameters parameters(argc, argv);
		switchLogLevel(parameters.get_log_level());

		fstream logfile;
		logfile.open("heatbath.log", std::ios::out | std::ios::app);
		if(logfile.is_open()) {
			meta::print_info_heatbath(argv[0], &logfile, parameters);
			logfile.close();
		} else {
			throw File_Exception("heatbath.log");
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		init_timer.reset();

		hardware::System system(parameters);

		PRNG prng(system);
		Gaugefield gaugefield(system, prng);

		print_gaugeobservables(gaugefield, 0);

		init_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Heatbath
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();

		logger.trace() << "Start thermalization" ;
		int ntherm = parameters.get_thermalizationsteps();
		for(int iter = 0; iter < ntherm; iter++) heatbath(gaugefield, prng);

		int nsteps = parameters.get_heatbathsteps();
		int overrelaxsteps = parameters.get_overrelaxsteps();
		int writefreq = parameters.get_writefrequency();
		int savefreq = parameters.get_savefrequency();

		logger.info() << "Start heatbath";

		for(int i = 0; i < nsteps; i++) {
			heatbath(gaugefield, prng, overrelaxsteps);
			if( ( (i + 1) % writefreq ) == 0 ) {
				//name of file to store gauge observables
				std::string gaugeout_name = meta::get_gauge_obs_file_name(parameters, "");
				print_gaugeobservables(gaugefield, i, gaugeout_name);
			}
			if( savefreq != 0 && ( (i + 1) % savefreq ) == 0 ) {
				gaugefield.save(i+1);
			}
		}

		gaugefield.save(nsteps);
		logger.trace() << "heatbath done";
		perform_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Final Output
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		total_timer.add();
		general_time_output(&total_timer, &init_timer, &perform_timer, &plaq_timer, &poly_timer);

	} //try
	//exceptions from Opencl classes
	catch (Opencl_Error& e) {
		logger.fatal() << e.what();
		exit(1);
	} catch (File_Exception& fe) {
		logger.fatal() << "Could not open file: " << fe.get_filename();
		logger.fatal() << "Aborting.";
		exit(1);
	} catch (Print_Error_Message& em) {
		logger.fatal() << em.what();
		exit(1);
	} catch (Invalid_Parameters& es) {
		logger.fatal() << es.what();
		exit(1);
	}

	return 0;

}
