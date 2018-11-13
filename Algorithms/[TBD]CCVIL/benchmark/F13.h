#ifndef _F13_H
#define _F13_H

#include "Benchmarks.h"

class F13:public Benchmarks{
protected:
public:
	F13(RunParameter* runParam);
	F13();
	double compute(double* x) ;
	double compute(vector<double> x) ;
	~F13();
};

#endif

