/*	
RANLUXCLTest, by Ivar Nikolaisen (2011).

This program show how to use the OpenCL implementation of RANLUX
(named RANLUXCL), and checks the correctness of the generated numbers 
against values generated by the original fortran 77 code. It does this 
by generating 1 000 000 values per work-item and checking the final number 
of each work-item agains those generated by the fortran implementation.
The reference values are loaded from the file ranluxclref.binary.
*/

//OpenCL header
#include "CL/cl.hpp"

//Standard headers
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <sys/stat.h> // to get size of file
#include <cstdlib>
#include <cmath>

template <class T>
inline std::string to_string(T x){
	std::ostringstream o;
	o << std::setprecision(16) << x;
	return o.str();
}

//iTimer is a simple class to make timing events easy. It should provide microsecond accuracy
//on both Windows and *nix systems.
#ifndef WIN32
#include <sys/time.h>
	//Implementation for *nix (does it work on BSD systems?)
	class iTimer{
	private:
		timeval StartTime, EndTime;
		double AccumulatedTime;
	public:
		iTimer(){
			gettimeofday(&StartTime, NULL);
			AccumulatedTime = 0;
		} //Constructor
		void Start(){gettimeofday(&StartTime, NULL);}

		double End(){
			gettimeofday(&EndTime, NULL); 
			AccumulatedTime += ((double)EndTime.tv_sec + EndTime.tv_usec/1E6 - ((double)StartTime.tv_sec + StartTime.tv_usec/1E6));
			return AccumulatedTime;
		}

		double Time(){return AccumulatedTime;}
	};

#elif defined WIN32
	#include <Windows.h>
	//Implementation for Windows
	class iTimer{
	private:
		LARGE_INTEGER StartTime, EndTime, Frequency, Temp;
		double AccumulatedTime;
	public:
		iTimer(){ //Constructor
			QueryPerformanceCounter(&StartTime);
			EndTime.QuadPart = 0;
			AccumulatedTime = 0;
			QueryPerformanceFrequency(&Frequency);
		}

		void Start(){QueryPerformanceCounter(&StartTime);}

		double End(){
			QueryPerformanceCounter(&EndTime);
			Temp.QuadPart = EndTime.QuadPart - StartTime.QuadPart;
			AccumulatedTime += (double)(Temp.QuadPart) / Frequency.QuadPart;
			EndTime.QuadPart = StartTime.QuadPart = 0;
			return AccumulatedTime;
		}

		double Time(){return AccumulatedTime;}
	};
#endif //WIN32
std::string Get_CL_Errstring(cl_int err) {
	//Get_CL_Errstring returns the error string associated with a given OpenCL error code.
	switch (err) {
		case CL_SUCCESS:                          return "Success!";
		case CL_DEVICE_NOT_FOUND:                 return "Device not found.";
		case CL_DEVICE_NOT_AVAILABLE:             return "Device not available";
		case CL_COMPILER_NOT_AVAILABLE:           return "Compiler not available";
		case CL_MEM_OBJECT_ALLOCATION_FAILURE:    return "Memory object allocation failure";
		case CL_OUT_OF_RESOURCES:                 return "Out of resources";
		case CL_OUT_OF_HOST_MEMORY:               return "Out of host memory";
		case CL_PROFILING_INFO_NOT_AVAILABLE:     return "Profiling information not available";
		case CL_MEM_COPY_OVERLAP:                 return "Memory copy overlap";
		case CL_IMAGE_FORMAT_MISMATCH:            return "Image format mismatch";
		case CL_IMAGE_FORMAT_NOT_SUPPORTED:       return "Image format not supported";
		case CL_BUILD_PROGRAM_FAILURE:            return "Program build failure";
		case CL_MAP_FAILURE:                      return "Map failure";
		case CL_INVALID_VALUE:                    return "Invalid value";
		case CL_INVALID_DEVICE_TYPE:              return "Invalid device type";
		case CL_INVALID_PLATFORM:                 return "Invalid platform";
		case CL_INVALID_DEVICE:                   return "Invalid device";
		case CL_INVALID_CONTEXT:                  return "Invalid context";
		case CL_INVALID_QUEUE_PROPERTIES:         return "Invalid queue properties";
		case CL_INVALID_COMMAND_QUEUE:            return "Invalid command queue";
		case CL_INVALID_HOST_PTR:                 return "Invalid host pointer";
		case CL_INVALID_MEM_OBJECT:               return "Invalid memory object";
		case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:  return "Invalid image format descriptor";
		case CL_INVALID_IMAGE_SIZE:               return "Invalid image size";
		case CL_INVALID_SAMPLER:                  return "Invalid sampler";
		case CL_INVALID_BINARY:                   return "Invalid binary";
		case CL_INVALID_BUILD_OPTIONS:            return "Invalid build options";
		case CL_INVALID_PROGRAM:                  return "Invalid program";
		case CL_INVALID_PROGRAM_EXECUTABLE:       return "Invalid program executable";
		case CL_INVALID_KERNEL_NAME:              return "Invalid kernel name";
		case CL_INVALID_KERNEL_DEFINITION:        return "Invalid kernel definition";
		case CL_INVALID_KERNEL:                   return "Invalid kernel";
		case CL_INVALID_ARG_INDEX:                return "Invalid argument index";
		case CL_INVALID_ARG_VALUE:                return "Invalid argument value";
		case CL_INVALID_ARG_SIZE:                 return "Invalid argument size";
		case CL_INVALID_KERNEL_ARGS:              return "Invalid kernel arguments";
		case CL_INVALID_WORK_DIMENSION:           return "Invalid work dimension";
		case CL_INVALID_WORK_GROUP_SIZE:          return "Invalid work group size";
		case CL_INVALID_WORK_ITEM_SIZE:           return "Invalid work item size";
		case CL_INVALID_GLOBAL_OFFSET:            return "Invalid global offset";
		case CL_INVALID_EVENT_WAIT_LIST:          return "Invalid event wait list";
		case CL_INVALID_EVENT:                    return "Invalid event";
		case CL_INVALID_OPERATION:                return "Invalid operation";
		case CL_INVALID_GL_OBJECT:                return "Invalid OpenGL object";
		case CL_INVALID_BUFFER_SIZE:              return "Invalid buffer size";
		case CL_INVALID_MIP_LEVEL:                return "Invalid mip-map level";
		default:                                  return "Unknown";
	}
}
void checkErr(cl_int err, const char *Message){
	//checkErr is a helper function to check for errors in OpenCL API calls,
	//and print a provided message along with the error code.
	if(err != CL_SUCCESS){
		std::cerr << "ERROR (" << err << "): " << Message << ": " << Get_CL_Errstring(err) << "\n";
		exit(EXIT_FAILURE);
	}
}
void Print_Help(){
	std::cout << "Usage:\n"
	          << "ranluxcltest.exe <lux> <UseGPU> <Buildlog> <OpenCLInfo> <deviceNr> <platformNr>\n\n"
	          << "Where: *lux is the luxury setting for ranluxcl. Set to 0-4, with\n"
	          << "        4 the best but slowest. Set to above 24 to manually set\n" 
	          << "        the p-value.\n\n"
	          << "       *UseGPU decides whether a GPU is to be used. If set to\n"
	          << "        true (1), a GPU device will be requested.\n\n"
	          << "       *If Buildlog is set true, the buildlog for the OpenCL code\n"
	          << "        will be printed\n\n"
	          << "       *If OpenCLInfo is set true, information about the OpenCL\n"
	          << "        platform and device will be displayed\n\n"
	          << "       *deviceNr selects the device (CPU or GPU depending on UseGPU\n"
	          << "        option), with 0 as the default (indexing starts at zero).\n\n"
	          << "       *platformNr selects the platform, with 0 as the default.\n"
	          << "        (indexing starts at zero)\n\n"
	          << "\nIf only the first argument is given, the rest default to false.\n"
	          << "\nFor example this will run with luxury setting 4 (slowest) on a GPU,\n"
	          << "displaying the build log and information about the OpenCL device and\n"
	          << "platform:\n"
	          << "ranluxcltest.exe 4 1 1 1\n"
	          << "\nCorrectness of final numbers will be checked if lux is 0-4.\n";
}

void Compile_OpenCL_Code(
	std::string BuildOptions,
	std::string FileName,
	cl::Context &context,
	cl::Program &program,
	cl_uint verbosity)
{
	//This function builds the OpenCL kernels and functions using the provided build options.
	//If verbosity is >0 the build log is printed.
	cl_int err, err2;
	cl::STRING_CLASS buildLog;

	//Getting source code
	std::ifstream file(FileName.c_str());
	checkErr(file.is_open() ? CL_SUCCESS:-1, FileName.c_str());
	std::string prog(std::istreambuf_iterator<char>(file), (std::istreambuf_iterator<char>()));
	cl::Program::Sources source(1, std::make_pair(prog.c_str(), prog.length()+1));
	program = cl::Program(context, source);
	file.close();

	//Get vector of devices in context
	VECTOR_CLASS <cl::Device> deviceList;
	err = context.getInfo(CL_CONTEXT_DEVICES, &deviceList);
	checkErr(err, "cl::Context::getInfo()");

	//Building source code
	err2 = program.build(deviceList, BuildOptions.c_str());
		
	//Fetch build log for first device
	err = program.getBuildInfo(deviceList[0], CL_PROGRAM_BUILD_LOG, &buildLog);
	checkErr(err, "getBuildInfo");

	if(verbosity)
		std::cout << "Build log for " << FileName << ":\n" << buildLog << "\n";

	checkErr(err2, "cl::Program::build()");
}

void OpenCL_Initializations(cl::Context &context, 
	                        std::vector <cl::Device> &devices,
	                        cl::CommandQueue &queue,
	                        bool PrintOpenCLInfo,
	                        bool UseOpenCLProfiling,
	                        bool UseGPU,
	                        cl_uint deviceNr,
	                        cl_uint platformNr)
{
	//Performs initializations of OpenCL

	std::vector <cl::Platform> platformList;
	cl_int err;

	//Getting list of platforms
	cl::Platform::get(&platformList);
	checkErr(platformList.size()!=0 ? CL_SUCCESS : -1, "cl::Platform::get");

	if(platformNr >= platformList.size()){
		std::cerr << "Platform nr. " << platformNr << " (indexing starts at zero) cannot be selected as there are only "
		          << platformList.size() << " platforms available. Exiting\n";
		exit(1);
	}

	cl_context_properties cprops[3] = {CL_CONTEXT_PLATFORM, (cl_context_properties)(platformList[platformNr])(), 0};

	//Creating context
	if(UseGPU) context = cl::Context(CL_DEVICE_TYPE_GPU, cprops, NULL, NULL, &err);
	else context = cl::Context(CL_DEVICE_TYPE_CPU, cprops, NULL, NULL, &err);
	checkErr(err, "Context::Context()");

	//Getting devices from context.
	devices = context.getInfo<CL_CONTEXT_DEVICES>();
	checkErr(devices.size() > 0 ? CL_SUCCESS : -1, "devices.size() > 0");

	if(deviceNr >= devices.size()){
		std::cerr << "Device nr. " << deviceNr << " (indexing starts at zero) cannot be selected as there are only "
		          << devices.size() << " devices available. Exiting\n";
		exit(1);
	}

	//Creating command queue
	if(UseOpenCLProfiling == 1) queue = cl::CommandQueue(context, devices[deviceNr], CL_QUEUE_PROFILING_ENABLE, &err);
	else queue = cl::CommandQueue(context, devices[deviceNr], 0, &err);
	checkErr(err, "CommandQueue::CommandQueue()");

	//Print some OpenCL information if requested
	if(PrintOpenCLInfo == 1){
		cl_uint uint;
		size_t sizet;
		cl_ulong ulong;
		std::string str;

		std::cout << "Num platforms:       " << platformList.size() << std::endl;

		std::string PlatformVendor;
		platformList[platformNr].getInfo(CL_PLATFORM_VENDOR, &PlatformVendor);
		std::cout << "Platform Vendor:     " << PlatformVendor << std::endl;

		platformList[platformNr].getInfo(CL_PLATFORM_PROFILE, &str);
		std::cout << "Platform profile:    " << str << "\n";

		platformList[platformNr].getInfo(CL_PLATFORM_VERSION, &str);
		std::cout << "Platform version:    " << str << "\n";

		platformList[platformNr].getInfo(CL_PLATFORM_NAME, &str);
		std::cout << "Platform name:       " << str << "\n";

		platformList[platformNr].getInfo(CL_PLATFORM_EXTENSIONS, &str);
		std::cout << "Platform extensions: " << str << "\n";

		devices[deviceNr].getInfo(CL_DEVICE_NAME, &str);
		std::cout << "Device name:         " << str << "\n";

		cl_device_type devicetype;
		devices[deviceNr].getInfo(CL_DEVICE_TYPE, &devicetype);
		std::cout << "Device type:         ";
		if(devicetype == CL_DEVICE_TYPE_CPU) std::cout << "CPU ";
		if(devicetype == CL_DEVICE_TYPE_GPU) std::cout << "GPU ";
		if(devicetype == CL_DEVICE_TYPE_ACCELERATOR) std::cout << "Accelerator ";
		if(devicetype == CL_DEVICE_TYPE_DEFAULT) std::cout << "Default ";
		std::cout << "\n";

		devices[deviceNr].getInfo(CL_DEVICE_VENDOR_ID, &uint);
		std::cout << "Device vendor ID:    " << uint << "\n";

		devices[deviceNr].getInfo(CL_DEVICE_MAX_COMPUTE_UNITS, &uint);
		std::cout << "Max compute units:   " << uint << "\n";

		devices[deviceNr].getInfo(CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, &uint);
		std::cout << "Max work item dim:   " << uint << "\n";

		devices[deviceNr].getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &sizet);
		std::cout << "Max work group size: " << sizet << "\n";

		devices[deviceNr].getInfo(CL_DEVICE_MAX_MEM_ALLOC_SIZE, &ulong);
		std::cout << "Max mem alloc size:  " << ulong/(1024*1024.0) << " MiB\n";

		devices[deviceNr].getInfo(CL_DEVICE_GLOBAL_MEM_SIZE, &ulong);
		std::cout << "Global memory size:  " << ulong/(1024*1024.0) << " MiB\n";

		devices[deviceNr].getInfo(CL_DEVICE_GLOBAL_MEM_CACHE_TYPE, &ulong);
		std::cout << "Global memory cache: ";
		if(ulong == CL_NONE) std::cout << "None\n";
		if(ulong == CL_READ_ONLY_CACHE) std::cout << "Read only\n";
		if(ulong == CL_READ_WRITE_CACHE) std::cout << "Read/Write cache\n";

		devices[deviceNr].getInfo(CL_DEVICE_LOCAL_MEM_TYPE, &ulong);
		std::cout << "Local memory type:   ";
		if(ulong == CL_LOCAL) std::cout << "Dedicated\n";
		if(ulong == CL_GLOBAL) std::cout << "Global memory\n";

		devices[deviceNr].getInfo(CL_DEVICE_PROFILING_TIMER_RESOLUTION, &sizet);
		std::cout << "Prof Timer Res:      " << sizet << " ns\n";

		devices[deviceNr].getInfo(CL_DRIVER_VERSION, &str);
		std::cout << "OpenCL Driver ver:   " << str << "\n";
	}
}
cl_float *ReadReferenceFile(size_t &NumVals, bool Verbose){
	//Reads reference values from binary file.

	std::string FileName;
	FileName = SOURCEDIR "/ranluxcl/ranluxclref.binary";

	// get file size using stat
	struct stat result;

	if(stat(FileName.c_str(), &result) == 0)
		NumVals = result.st_size/sizeof(cl_float);
	else
		return NULL;
		

	cl_float *refvals = new cl_float [(size_t)NumVals];

	//Get values
	std::ostringstream ost;
	if(Verbose) std::cout << "Reading file: " << FileName << "\n";

	stat(FileName.c_str(), &result);

	ost << FileName.c_str();
	std::ifstream infile;
	infile.open(ost.str().c_str(), std::ios::in | std::ios::binary);
	infile.read((char*)refvals,result.st_size);
	infile.close();

	return refvals;
}
void CheckCorrectness(size_t NumWorkitems,
                      int lux,
                      cl_int KernelCycles,
                      cl_int NumIterations,
                      cl_float *PRNs)
{
	//Checking correctness. A binary file containing reference values generated by the original fortran code
	//is loaded. The file contains value number 1 million for seeds 1 through 51200 for p = 24, 48, 100, 224 and 404
	//(i.e. luxury 0, 1, 2, 3 and 4 of ranluxcl_initialization).

	bool doCorrectnessCheck = 1;
	size_t NumVals, numRefFileValsPerpval=51200, NumValsToCheck=NumWorkitems, numBackupValsPerpval=4, offset=0, numErrors;
	double errorLimit=1E-6, diffSum, avgRef, avgGen;

	std::cout << "\n";

	//The number of values to check is either the number of work-items, or numRefFileValsPerpval if
	//there are more work-items than there are values.
	if(NumWorkitems > numRefFileValsPerpval)
		NumValsToCheck = numRefFileValsPerpval;

	//The values of each of the five p-values are stored one after the other. This sets the offset we need to index from.
	if(lux == 0) offset = 0;
	else if(lux == 1) offset = 1;
	else if(lux == 2) offset = 2;
	else if(lux == 3) offset = 3;
	else if(lux == 4) offset = 4;
	else{
		std::cout << "luxury value not checkable, correctness check will not be performed\n\n";
		doCorrectnessCheck = 0;
	}

	//Correctness check is only valid if exactly 1 million values have been generated per work-item
	if(KernelCycles * NumIterations != 1000000){
		std::cout << "Correctness check only possible if 1000000 values are generated per work-item, but\n" 
		          << KernelCycles * NumIterations << " have been generated. Correctness check will not be performed\n\n";
		doCorrectnessCheck = 0;
	}

	//Performing correctness check
	if(doCorrectnessCheck){
		std::cout << "Performing correctness check. Done in several different ways to avoid errors.\n";

		cl_float *refVals = ReadReferenceFile(NumVals, 1);

		if(refVals == NULL){
			std::cout << "Reference file ranluxclref.binary not found, using backup dataset.\n";
			NumValsToCheck = NumWorkitems;
			if(NumWorkitems > numBackupValsPerpval)
				NumValsToCheck = numBackupValsPerpval;

			refVals = new cl_float [numBackupValsPerpval * 5];

			offset *= numBackupValsPerpval;

			refVals[0+0*numBackupValsPerpval]=(cl_float)0.82777202129364; refVals[1+0*numBackupValsPerpval]=(cl_float)0.14299058914185; 
			refVals[2+0*numBackupValsPerpval]=(cl_float)0.98017895221710; refVals[3+0*numBackupValsPerpval]=(cl_float)0.96829599142075;

			refVals[0+1*numBackupValsPerpval]=(cl_float)0.32378226518631; refVals[1+1*numBackupValsPerpval]=(cl_float)0.68649446964264; 
			refVals[2+1*numBackupValsPerpval]=(cl_float)0.89647895097733; refVals[3+1*numBackupValsPerpval]=(cl_float)0.42210841178894;

			refVals[0+2*numBackupValsPerpval]=(cl_float)0.74214434623718; refVals[1+2*numBackupValsPerpval]=(cl_float)0.58056205511093; 
			refVals[2+2*numBackupValsPerpval]=(cl_float)0.77871483564377; refVals[3+2*numBackupValsPerpval]=(cl_float)0.72025460004807;

			refVals[0+3*numBackupValsPerpval]=(cl_float)0.71860218048096; refVals[1+3*numBackupValsPerpval]=(cl_float)0.76306450366974; 
			refVals[2+3*numBackupValsPerpval]=(cl_float)0.56076055765152; refVals[3+3*numBackupValsPerpval]=(cl_float)0.45423328876495;

			refVals[0+4*numBackupValsPerpval]=(cl_float)0.42760473489761; refVals[1+4*numBackupValsPerpval]=(cl_float)0.10520881414414; 
			refVals[2+4*numBackupValsPerpval]=(cl_float)0.30440622568131; refVals[3+4*numBackupValsPerpval]=(cl_float)0.50632297992706;
		}
		else offset *= numRefFileValsPerpval;

		//Check each value
		numErrors = 0;
		diffSum = avgRef = avgGen = 0;
		for(size_t i=0; i<NumValsToCheck; i++){
			avgRef += refVals[offset + i];
			avgGen += PRNs[i];
			diffSum += std::abs((double)(PRNs[i] - refVals[offset + i]));
			if(std::abs((double)(PRNs[i] - refVals[offset + i])) > errorLimit)
				numErrors++;
		}
		avgRef = avgRef / (double)NumValsToCheck;
		avgGen = avgGen / (double)NumValsToCheck;

		if(numErrors == 0)
			std::cout << "\nCorrectness check PASSED!\n" << "Please also verify that the below statements are correct:\n";
		else
			std::cout << "\nCorrectness check FAILED!\n" << "Presumably the below statements will also not be true:\n";

		std::cout << std::setprecision(6) << std::fixed
		          << "Number of errors found in " << NumValsToCheck << " values (should be zero): " << numErrors << "\n"
		          << "Total summed up difference of all values (should be zero): " << diffSum << "\n"
		          << "The averages of the reference and generated numbers (should be equal): " << avgRef << " and " << avgGen << "\n\n";
	}
}

int main(int argc, char *argv[])
{
	iTimer totalTimer;
	double genTime;
	cl_int err;
	cl_ulong TotalNumbersGenerated;

	//OpenCL variables
	cl::Context context;
	std::vector <cl::Device> devices;
	cl::CommandQueue queue;
	cl::Program program;

	/*Below are some program settings.
	KernelCycles is the number of cycles done per kernel 
	(must be divisible by 4 because each call to ranluxcl produces a float4).

	NumWorkitems is the number of work-items (kernels) to launch per iteration. 
	For correctness check to be applicable KernelCycles * NumIterations must
	equal 1 000 000. */

	//Default values for GPU. Feel free to change anything as long as 
	//KernelCycles * NumIterations == 1 000 000 so that correctness check
	//can be performed. If running out of resources on nvidia hardware
	//try a work-group size of 32. If execution is taking a long time,
	//reduce NumWorkitems.
	cl_int KernelCycles = 1000; //Should be divisible by 4
	size_t WorkgroupSize = 64;
	size_t NumWorkitems = 800 * WorkgroupSize;
	cl_int NumIterations = 1000;
	cl_int deviceNr = 0;
	cl_int platformNr = 0;
	bool ShowBuildLog = 0;
	bool PrintOpenCLInfo = 0;
	bool UseGPU = 0;

	//Reduce amount of work a bit if CPU is selected
	size_t WorkgroupSizeCPU = 1;
	size_t NumWorkitemsCPU = 512;

	if(KernelCycles % 4 != 0){
		std::cout << "KernelCycles not divisible by 4. Exiting\n";
		return 1;
	}

	//Reading parameters
	if(argc < 2){
		Print_Help();
		exit(0);
	}

	//Input values handled strangely to avoid warnings in Visual Studio 10 when casting int to bool.
	int lux = atoi(argv[1]);
	if(argc > 2) UseGPU = atoi(argv[2]) ? 1 : 0;
	if(argc > 3) ShowBuildLog = atoi(argv[3]) ? 1 : 0;
	if(argc > 4) PrintOpenCLInfo = atoi(argv[4]) ? 1 : 0;
	if(argc > 5) deviceNr = atoi(argv[5]);
	if(argc > 6) platformNr = atoi(argv[6]);

	//Different settings when using CPU
	if(!UseGPU){
		WorkgroupSize = WorkgroupSizeCPU;
		NumWorkitems = NumWorkitemsCPU;
	}

	TotalNumbersGenerated = cl_ulong(KernelCycles) * cl_ulong(NumWorkitems) * cl_ulong(NumIterations);

	//Size of ranluxcltab buffer is seven float4 per work-item.
	size_t RANLUXCLTabSize = NumWorkitems * 7 * sizeof(cl_float4);

	if(UseGPU) std::cout << "Using GPU\n";
	else std::cout << "Using CPU\n";

	cl_float* PRNs = new cl_float [NumWorkitems];

	OpenCL_Initializations(context, devices, queue, PrintOpenCLInfo, 0, UseGPU, deviceNr, platformNr);

	//Building OpenCL program
	std::cout << "Building OpenCL program: ";
	std::string BuildOptions;
	std::string BuildLog;
	std::string FileName = SOURCEDIR "/ranluxcl/ranluxcltest_kernels.cl";

	BuildOptions += "-I " SOURCEDIR "/ranluxcl"; //Search for include files in current directory

	//Used to get exact same sequence as original implementation (for correctness check). Should generally
	//NOT be set in other programs, to ensure parallel generators aren't in initially correlated states.
	BuildOptions += " -D RANLUXCL_NO_WARMUP";

	//Set luxury value. If this is not defined the highest (4) is used by default.
	BuildOptions += " -D RANLUXCL_LUX=" + to_string(lux);

	Compile_OpenCL_Code(BuildOptions, FileName, context, program, ShowBuildLog);

	//Setting kernels
	cl::Kernel Kernel_RanluxclInit = cl::Kernel(program, "Kernel_Ranluxcl_Init", &err); checkErr(err, "cl::Kernel Kernel_RanluxclInit");
	cl::Kernel Kernel_PRN = cl::Kernel(program, "Kernel_PRN", &err); checkErr(err, "cl::Kernel kernel_RanluxclInit");

	//Creating OpenCL buffers
	cl::Buffer Buffer_RANLUXCLTab(context, CL_MEM_READ_WRITE, RANLUXCLTabSize, NULL, &err);
	checkErr(err, "Buffer:RANLUXCLTab");
	cl::Buffer Buffer_PRNs(context, CL_MEM_READ_WRITE, NumWorkitems * sizeof(cl_float), NULL, &err);
	checkErr(err, "Buffer:PRNs");

	//Setting arguments to kernels
	err = Kernel_PRN.setArg(0, KernelCycles); checkErr(err, "Kernel::setArg()");
	err = Kernel_PRN.setArg(1, Buffer_RANLUXCLTab); checkErr(err, "Kernel::setArg()");
	err = Kernel_PRN.setArg(2, Buffer_PRNs); checkErr(err, "Kernel::setArg()");

	err = Kernel_RanluxclInit.setArg(0, 0); checkErr(err, "Kernel::setArg()");
	err = Kernel_RanluxclInit.setArg(1, Buffer_RANLUXCLTab); checkErr(err, "Kernel::setArg()");

	//Initialize the generator
	err = queue.enqueueNDRangeKernel(Kernel_RanluxclInit, cl::NullRange, cl::NDRange(NumWorkitems), cl::NDRange(WorkgroupSize));
	checkErr(err, "CommandQueue::enqueueNDRangeKernel()");

	iTimer genTimer;

	for(int Iteration=0; Iteration<NumIterations; Iteration++){
		err = queue.enqueueNDRangeKernel(Kernel_PRN, 
		                                 cl::NullRange, 
		                                 cl::NDRange(NumWorkitems), 
		                                 cl::NDRange(WorkgroupSize));
		checkErr(err, "CommandQueue::enqueueNDRangeKernel()");
	}

	//Just call finish after all commands have been entered into queue
	err = queue.finish(); checkErr(err, "queue.finish()");

	genTimer.End();
	genTime = genTimer.Time();
	
	//Reading final numbers generated
	err = queue.enqueueReadBuffer(Buffer_PRNs, 
	                              CL_TRUE, 
	                              0, 
	                              NumWorkitems * sizeof(cl_float), 
	                              PRNs);
	checkErr(err, "Read PRNs Buffer");

	//Check correctness
	CheckCorrectness(NumWorkitems, lux, KernelCycles, NumIterations, PRNs);

	totalTimer.End();

	std::cout << std::setprecision(2) << std::fixed
	          << "Total program time:      " << totalTimer.Time() << " seconds\n"
	          << "Generation time:         " << genTime << " seconds\n"
	          << "Total numbers generated: " << TotalNumbersGenerated/1E6 << " million\n" 
	          << "Generation speed:        " << TotalNumbersGenerated / (1E6 * genTime) << " million numbers per second" << "\n";

	return 0;
}
