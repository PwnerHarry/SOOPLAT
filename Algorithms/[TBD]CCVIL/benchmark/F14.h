#ifndef _F14_H
#define _F14_H

#include "Benchmarks.h"

class F14:public Benchmarks{
protected:
	void generateInterArray ( );
public:
	F14(RunParameter* runParam);
	F14();
	double compute(double* x) ;
	double compute(vector<double> x) ;
	~F14();
};

#endif

