#ifndef _F15_H
#define _F15_H

#include "Benchmarks.h"

class F15:public Benchmarks{
protected:
	void generateInterArray ( );
public:
	F15(RunParameter* runParam);
	F15();
	double compute(double* x) ;
	double compute(vector<double> x) ;
	~F15();
};

#endif
