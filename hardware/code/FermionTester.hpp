#ifndef FERMIONTESTER_HPP
#define FERMIONTESTER_HPP

#include "SpinorTester.hpp"

class FermionTester : public SpinorTester
{
public:
	FermionTester(std::string kernelName, std::string inputfileIn, int numberOfValues = 1):
	SpinorTester(kernelName, getSpecificInputfile(inputfileIn), numberOfValues)
	{
			code = device->get_fermion_code();
			gaugefield = new physics::lattices::Gaugefield(*system, *prng);
	}
	~FermionTester()
	{
		delete gaugefield;
	}
	
protected:
	const hardware::code::Fermions * code;
	physics::lattices::Gaugefield * gaugefield;
	
	std::string getSpecificInputfile(std::string inputfileIn)
	{
		//todo: this is ugly, find a better solution.
		// The problem is that the parent class calls a similar fct.
		return "../fermions/" + inputfileIn;
	}
	
	const hardware::buffers::SU3* getGaugefieldBuffer() {
		return gaugefield->get_buffers()[0];
	}
};

#endif