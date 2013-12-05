#include "generalExecutable.h"

#include "meta/util.hpp"
#include "physics/prng.hpp"
#include "physics/lattices/gaugefield.hpp"
#include "physics/algorithms/heatbath.hpp"

class heatbathExecutable: public generalExecutable
{
public:
	heatbathExecutable(int argc, const char* argv[]) : generalExecutable(argc, argv)
	{
		initializationTimer.reset();
		prng = new physics::PRNG(*system);
		gaugefield = new physics::lattices::Gaugefield(*system, *prng);
		meta::print_info_heatbath(ownName, parameters);
		writeHeatbathLogfile();
		initializationTimer.add();
	}

	void writeHeatbathLogfile() {
		outputToFile.open(filenameForInverterLogfile, std::ios::out | std::ios::app);
		if (outputToFile.is_open()) {
			meta::print_info_heatbath(ownName, &outputToFile, parameters);
			outputToFile.close();
		} else {
			throw File_Exception(filenameForInverterLogfile);
		}
	}

	void invoke(int argc, const char* argv[]) {
		using namespace physics;
		using namespace physics::lattices;
		using physics::algorithms::heatbath;
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Heatbath
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		performanceTimer.reset();
		print_gaugeobservables(*gaugefield, 0);
		logger.trace() << "Start thermalization";
		int ntherm = parameters.get_thermalizationsteps();
		for (int iter = 0; iter < ntherm; iter++)
			heatbath(*gaugefield, *prng);
		int nsteps = parameters.get_heatbathsteps();
		int overrelaxsteps = parameters.get_overrelaxsteps();
		int writefreq = parameters.get_writefrequency();
		int savefreq = parameters.get_savefrequency();
		logger.info() << "Start heatbath";
		for (int i = 0; i < nsteps; i++) {
			heatbath(*gaugefield, *prng, overrelaxsteps);
			if (((i + 1) % writefreq) == 0) {
				//name of file to store gauge observables
				std::string gaugeout_name = meta::get_gauge_obs_file_name(
						parameters, "");
				print_gaugeobservables(*gaugefield, i, gaugeout_name);
			}
			if (savefreq != 0 && ((i + 1) % savefreq) == 0) {
				gaugefield->save(i + 1);
			}
		}

		gaugefield->save(nsteps);
		logger.trace() << "heatbath done";
		performanceTimer.add();
	}
private:
	physics::PRNG * prng;
	physics::lattices::Gaugefield * gaugefield;
	const std::string 	filenameForInverterLogfile 		= "heatbath.log";
};

int main(int argc, const char* argv[])
{
	try {
		heatbathExecutable heatbathInstance(argc, argv);
		heatbathInstance.invoke(argc, argv);
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
