/*
 * @file
 * Declaration of the generalExecutable class.
 * This class provides general features needed by the executables.
 */

#ifndef GENERALEXECUTABLE_H_
#define GENERALEXECUTABLE_H_

#include <fstream>

#include "types.h"
#include "host_use_timer.h"
#include "exceptions.h"
#include "logger.hpp"
#include "meta/util.hpp"
#include "./hardware/system.hpp"
#include "physics/lattices/gaugefield.hpp"

class generalExecutable
{

public:
	/**
	 * Initialize meta::Inputparametes and Shardware::System objects
	 */
	generalExecutable(int argc, const char* argv[]);

	/**
	 * Destructor. Prints runtime information to screen and file.
	 */
	~generalExecutable();

protected:
	const char* ownName;
	const std::string 	filenameForCurrentPrngState 	= "prng.inverter.save";
	usetimer totalRuntimeOfExecutable;
	usetimer initializationTimer;
	usetimer performanceTimer;
	meta::Inputparameters parameters;
	hardware::System * system;
	physics::PRNG * prng;
	physics::lattices::Gaugefield * gaugefield;
	std::ofstream outputToFile;
	const char* generalTimeOutputFilename = "general_time_output";

	void printRuntimeInformationToScreenAndFile();

	void printGeneralTimesToScreen();

	void printGeneralTimesToFile();

	void saveCurrentPrngStateToFile();
};


class multipleConfigurationExecutable : public generalExecutable
{
public:
	multipleConfigurationExecutable(int argc, const char* argv[]) : generalExecutable (argc, argv)
	{
		initializationTimer.reset();
		setIterationVariables();
		initializationTimer.add();
	}

	void performMeasurements();

protected:
	std::string currentConfigurationName;
	int iterationStart;
	int iterationEnd;
	int iterationIncrement;

	void setIterationVariables();

	void initializeGaugefieldAccordingToIterationVariable(int interation);

	void initializeGaugefieldAccordingToConfigurationGivenInSourcefileParameter();

	void initializeGaugefield(int interation);

	void performMeasurementsForSpecificIteration(int interation);

	virtual void performApplicationSpecificMeasurements() {};
};

#endif /* GENERALEXECUTABLE_H_ */
