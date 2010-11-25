#include "opencl.h"

using namespace std;

hmc_error init_opencl(){
  cout<<"OpenCL being initialized..."<<endl;
  vector<cl::Platform> platforms;
  vector<cl::Device> devices;
  cl_int clerr = CL_SUCCESS;
  clerr = cl::Platform::get(&platforms);
  cout<<"Number of platforms:    "<<platforms.size()<<endl;
  //LZ: for now, stick to platform 0 without any further checks...
  cout<<"Choose platform number: 0"<<endl;
  cl::Platform platform = platforms[0];
  string info;
  platform.getInfo(CL_PLATFORM_NAME,&info);
  cout<<"\tCL_PLATFORM_NAME:     "<<info<<endl;
  platform.getInfo(CL_PLATFORM_VENDOR,&info);
  cout<<"\tCL_PLATFORM_VENDOR:   "<<info<<endl;
  platform.getInfo(CL_PLATFORM_VERSION,&info);
  cout<<"\tCL_PLATFORM_VERSION:  "<<info<<endl;
  //    platform.getInfo(CL_PLATFORM_PROFILE,&info);
  //    cout<<"\tCL_PLATFORM_PROFILE:  "<<info<<endl;
  cout<<endl;
  platform.getDevices(CL_DEVICE_TYPE_ALL,&devices);
  cout<<"\tAvailable devices: "<<devices.size()<<endl;
  int ndev = 0;
  for(unsigned int m=0; m<devices.size(); m++) {
    string info;
    devices[m].getInfo(CL_DEVICE_NAME,&info);
    cout<<"\t\tCL_DEVICE_NAME:    "<<info<<endl;
    devices[m].getInfo(CL_DEVICE_VENDOR,&info);
    cout<<"\t\tCL_DEVICE_VENDOR:  "<<info<<endl;
    cl_device_type type;
    devices[m].getInfo(CL_DEVICE_TYPE,&type);
    if(type == CL_DEVICE_TYPE_CPU) cout<<"\t\tCL_DEVICE_TYPE:    CPU"<<endl;
    if(type == CL_DEVICE_TYPE_GPU) cout<<"\t\tCL_DEVICE_TYPE:    GPU"<<endl;
    if(type == CL_DEVICE_TYPE_ACCELERATOR) cout<<"\t\tCL_DEVICE_TYPE:    ACCELERATOR"<<endl;
    if(type != CL_DEVICE_TYPE_CPU && type != CL_DEVICE_TYPE_GPU && type != CL_DEVICE_TYPE_ACCELERATOR) {
      cout<<"unexpected CL_DEVICE_TYPE..."<<endl;
      exit(HMC_OCLERROR);
    }
    if(type == wanted_device) ndev = m;
    devices[m].getInfo(CL_DEVICE_VERSION,&info);
    cout<<"\t\tCL_DEVICE_VERSION: "<<info<<endl;
  }
  cout<<"Choose device number:   "<<ndev<<endl;
  cout<<"Create context..."<<endl;
  cl_context_properties contextproperties[3] = {CL_CONTEXT_PLATFORM,(cl_context_properties)(platform)(),0};
  cl::Context context(wanted_device,contextproperties,NULL,NULL,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  cout<<"Create command queue..."<<endl;
  cl::CommandQueue cmdqueue(context,devices[ndev],0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  cout<<"Read kernel source from file: "<<cl_kernels_file<<endl;
  fstream kernelsfile;
  kernelsfile.open(cl_kernels_file.c_str());


  string kernelssource;

  //  kernelsfile.read(kernelssource;
  cout<<kernelssource<<endl;


  return HMC_SUCCESS;
}

hmc_error finalize_opencl(){
  return HMC_SUCCESS;
}
