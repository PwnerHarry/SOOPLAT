#ifndef _F11_H
#define _F11_H

#include "Benchmarks.h"

class F11:public Benchmarks{
protected:
public:
	F11(RunParameter* runParam);
	F11();
	double compute(double* x) ;
	double compute(vector<double> x) ;
	~F11();
};

#endif
