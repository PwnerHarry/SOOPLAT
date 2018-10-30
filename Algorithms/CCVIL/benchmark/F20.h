#ifndef _F20_H
#define _F20_H

#include "Benchmarks.h"

class F20:public Benchmarks{
protected:
public:
	F20(RunParameter* runParam);
	F20();
	double compute(double* x) ;
	double compute(vector<double> x) ;
	~F20();
};

#endif

