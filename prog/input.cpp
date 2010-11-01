#include "input.h"

hmc_error inputparameters::set_defaults(){
  kappa = 0.125;
  beta = 4.0;
  mu = 0.006;
  cgmax = 1000;
  return HMC_SUCCESS;
}

hmc_error inputparameters::readfile(char* ifn){
  std::ifstream infile;
  infile.open(ifn);
  if(!infile.is_open()) {
    printf("Could not open input file: %s\n",ifn);
    exit(HMC_FILEERROR);
  }
  while (infile.good()) {
    std::string line;
    infile>>line;
    if(line.find("#")!=std::string::npos) break; //allow comments
    if(line.find("kappa")!=std::string::npos) val_assign(&kappa,line);
    if(line.find("Kappa")!=std::string::npos) val_assign(&kappa,line);
    if(line.find("mu")!=std::string::npos) val_assign(&mu,line);
    if(line.find("Mu")!=std::string::npos) val_assign(&mu,line);
    if(line.find("beta")!=std::string::npos) val_assign(&beta,line);
    if(line.find("Beta")!=std::string::npos) val_assign(&beta,line);
    if(line.find("cgmax")!=std::string::npos) val_assign(&cgmax,line);
    if(line.find("CGmax")!=std::string::npos) val_assign(&cgmax,line);
    if(line.find("Cgmax")!=std::string::npos) val_assign(&cgmax,line);
  }
  return HMC_SUCCESS;
}

void inputparameters::val_assign(hmc_float * out, std::string line) {
  size_t pos = line.find("=");
  std::string value=line.substr(pos+1);
  (*out) = atof(value.c_str());
  return;
}

void inputparameters::val_assign(int * out, std::string line) {
  size_t pos = line.find("=");
  std::string value=line.substr(pos+1);
  (*out) = atoi(value.c_str());
  return;
}


hmc_float inputparameters::get_kappa(){
  return kappa;
}

hmc_float inputparameters::get_beta(){
  return beta;
}

hmc_float inputparameters::get_mu(){
  return mu;
}

int inputparameters::get_cgmax(){
  return cgmax;
}
