#ifndef _USETIMERH_
#define _USETIMERH_

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <sstream>
//#include <boost/lexical_cast.hpp>

#include "globaldefs.h"
#include "input.h"
#include "hmcerrs.h"
#include "timer.h" 
#include "gaugeobservables.h"

class usetimer {
 public:
  usetimer(inputparameters* parameters_) : parameters(parameters_) {};
  void start();
  void reset();
  void measurements(int n);
  void getTimeAndReset(int n);
  void add(int n);
  void zero(int n);
  void output();
 private:
  int time_measurements[10];
  Timer timer;
  inputparameters *parameters;
};

#endif
