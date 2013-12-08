/**
 * @file
 * Declaration of the generationExecutable class.
 * This class provides features for the generation of gauge configurations
 * according to certain algorithms.
 **/

#include "generalExecutable.h"

class generationExecutable : public generalExecutable
{
public:
	generationExecutable(int argc, const char* argv[]);

	void generateConfigurations();

protected:
	int writeFrequency;
	int saveFrequency;
	int thermalizationSteps;
	int generationSteps;
	int iteration;
	std::string filenameForGaugeobservables;

	/**
	 * Sets member variables that control the iterations during
	 * the generation of gaugefield configurations.
	 **/
	void setIterationParameters();

	/**
	 * Saves current gaugefield configuration to disk if current iteration
	 * is a multiple of the inputparameter save_frequency or if it is the
	 * last iteration.
	 **/
	void saveGaugefield();

	/**
	 * Saves current prng configuration to disk if current iteration
	 * is a multiple of the inputparameter save_frequency or if it is the
	 * last iteration.
	 **/
	void savePrng();

	/**
	 * Performs thermalization of the physical system according to the algorithm
	 * specified in the "thermalizeAccordingToSpecificAlgorithm" function.
	 **/
	void thermalize();

	/**
	 * Generates a new gauge configuration according to the algorithm
	 * specified in the "generateAccordingToSpecificAlgorithm" function.
	 * Performs measurements on this configuration according to the
	 * function "performOnlineMeasurements".
	 **/
	void generate();

	void virtual thermalizeAccordingToSpecificAlgorithm() {};

	void virtual generateAccordingToSpecificAlgorithm() {};

	/**
	 * Measurements to be performed after each step of configuration generation.
	 * By default this measures the gauge observables.
	 * This function can be replaced by a more specific one for each specific algorithm.
	 */
	void virtual performOnlineMeasurements();
};


