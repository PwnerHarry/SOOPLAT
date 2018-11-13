#ifndef _F6_H
#define _F6_H

#include "Benchmarks.h"

class F6:public Benchmarks{
protected:
public:
	F6(RunParameter* runParam);
	F6();
	double compute(double* x) ;
	double compute(vector<double> x) ;
	~F6();
};

#endif
