#ifndef _F8_H
#define _F8_H

#include "Benchmarks.h"

class F8:public Benchmarks{
protected:
public:
	F8(RunParameter* runParam);
	F8();
	double compute(double* x) ;
	double compute(vector<double> x) ;
	~F8();
};

#endif
