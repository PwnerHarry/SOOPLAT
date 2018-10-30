#ifndef _F9_H
#define _F9_H

#include "Benchmarks.h"

class F9:public Benchmarks{
protected:
public:
	F9(RunParameter* runParam);
	F9();
	double compute(double* x) ;
	double compute(vector<double> x) ;
	~F9();
};

#endif
