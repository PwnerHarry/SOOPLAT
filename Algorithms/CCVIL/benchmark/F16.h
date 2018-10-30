#ifndef _F16_H
#define _F16_H

#include "Benchmarks.h"

class F16:public Benchmarks{
protected:
	void generateInterArray ( );
public:
	F16(RunParameter* runParam);
	F16();
	double compute(double* x) ;
	double compute(vector<double> x) ;
	~F16();
};

#endif

