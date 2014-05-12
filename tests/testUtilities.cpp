#include "testUtilities.hpp"

#include "../host_functionality/logger.hpp"
#include <boost/test/unit_test.hpp>

std::string defaultGpuOption = "--use_gpu=false";
std::string defaultRec12Option = "--use_rec12=false";
std::string defaultSourceDirectory = "../../tests/inputfiles";

static void setArguments(std::string & inputfile_location, std::string & gpu_opt, std::string & rec12_opt, int & num_par, const int param_expect)
{
  logger.info() << "expect command line parameters:";
  logger.info() << "\t<exec_name>\t<source-dir>\t<gpu_usage>\t<rec12_usage>";

	switch(num_par){
		case 0:
			logger.fatal() << "Something went terribly wrong! Did not even get executable name! Aborting...";
			exit(-1);
			break;
		case 1:
			logger.fatal() << "Got only " << num_par << " command line parameters, expected " << param_expect << "! Use default values instead:";
			logger.fatal() << "\t<exec_name>\t" << defaultSourceDirectory << "\t" << defaultGpuOption << "\t" << defaultRec12Option;
			break;
		case 2:
			logger.fatal() << "Got only " << num_par << " command line parameters, expected " << param_expect << "! Use default values instead:";
			logger.fatal() << "\t<exec_name>\t<source-dir>\t" << defaultGpuOption << "\t" << defaultRec12Option;
			
			inputfile_location = boost::unit_test::framework::master_test_suite().argv[1];
			break;
		case 3:
			logger.fatal() << "Got only " << num_par << " command line parameters, expected " << param_expect << "! Use default values instead:";
			logger.fatal() << "\t<exec_name>\t<source-dir>\t<gpu-usage>\t" << defaultRec12Option;
			
			inputfile_location = boost::unit_test::framework::master_test_suite().argv[1];
			gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
			break;

		default:
			if(num_par > param_expect)
			{
			  logger.warn() << "Got " << num_par << " command line parameters, expected only " << param_expect << "! Use only the first " << param_expect << " values";
			  num_par = 4;
			}
			
			inputfile_location = boost::unit_test::framework::master_test_suite().argv[1];
			gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
			rec12_opt =  boost::unit_test::framework::master_test_suite().argv[3];
			break;
	}
}

meta::Inputparameters createParameters(std::string inputfile)
{
	std::string inputfile_location = defaultSourceDirectory;
	std::string gpu_opt = defaultGpuOption;
	std::string rec12_opt = defaultRec12Option;
	int num_par = 0;
	const int param_expect = 4;
  
	num_par = boost::unit_test::framework::master_test_suite().argc;
	setArguments(inputfile_location, gpu_opt, rec12_opt, num_par, param_expect);
	inputfile_location += '/' + inputfile;
	logger.info() << "inputfile used: " << inputfile_location;
	
	const char* _params_cpu[] = {"foo", inputfile_location.c_str(), gpu_opt.c_str() , rec12_opt.c_str(), "--device=0"};
	meta::Inputparameters params(num_par + 1, _params_cpu);
	return params;
}

void printKernelInformation(std::string name)
{
  logger.info() << "Test kernel\t\"" << name << "\"\tagainst reference value";
}
