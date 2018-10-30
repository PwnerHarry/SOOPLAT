#ifndef _F19_H
#define _F19_H

#include "Benchmarks.h"

class F19:public Benchmarks{
protected:
public:
	F19(RunParameter* runParam);
	F19();
	double compute(double* x) ;
	double compute(vector<double> x) ;
	~F19();
};

#endif


