#ifndef _F10_H
#define _F10_H

#include "Benchmarks.h"

class F10:public Benchmarks{
protected:
public:
	F10(RunParameter* runParam);
	F10();
	double compute(double* x) ;
	double compute(vector<double> x) ;
	~F10();
};

#endif
