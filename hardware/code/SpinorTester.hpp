#include "../hardware/code/kernelTester.hpp"

#include "../../meta/util.hpp"
#include "../../host_functionality/host_random.h"
#include "../../physics/prng.hpp"
#include "spinors.hpp"
#include "complex.hpp"

class SpinorTester : public KernelTester {
public:
	SpinorTester(std::string kernelName, std::string inputfileIn, int numberOfValues = 1);
	~SpinorTester();
	
protected:
	std::string getSpecificInputfile(std::string inputfileIn);

	spinor * createSpinorfield(size_t numberOfElements, int seed = 123456);
	void fill_with_one(spinor * in, int size);
	void fill_with_random(spinor * in, int size, int seed);
	spinor * createSpinorfieldWithOnesAndZerosDependingOnSiteParity();	
	void fill_with_one_eo(spinor * in, int size, bool eo);
	hmc_float count_sf(spinor * in, int size);
	hmc_float calc_var(hmc_float in, hmc_float mean);
	hmc_float calc_var_sf(spinor * in, int size, hmc_float sum);
	void calcSquarenormAndStoreAsKernelResult(const hardware::buffers::Plain<spinor> * in);
	void calcSquarenormEvenOddAndStoreAsKernelResult(const hardware::buffers::Spinor * in);
	
	const hardware::code::Spinors * code;
	physics::PRNG * prng;

	hardware::buffers::Plain<double> * doubleBuffer;
	
	size_t spinorfieldElements;
	size_t spinorfieldEvenOddElements;
	bool useRandom;
	bool evenOrOdd;
	bool calcVariance;
	hmc_complex alpha_host;
	hmc_complex beta_host;
	int iterations;
};