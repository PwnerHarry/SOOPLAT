#ifndef _F17_H
#define _F17_H

#include "Benchmarks.h"

class F17:public Benchmarks{
protected:
	void generateInterArray ( );
public:
	F17(RunParameter* runParam);
	F17();
	double compute(double* x) ;
	double compute(vector<double> x) ;
	~F17();
};

#endif

