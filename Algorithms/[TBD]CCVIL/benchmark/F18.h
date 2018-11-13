#ifndef _F18_H
#define _F18_H

#include "Benchmarks.h"

class F18:public Benchmarks{
protected:
	void generateInterArray ( );
public:
	F18(RunParameter* runParam);
	F18();
	double compute(double* x) ;
	double compute(vector<double> x) ;
	~F18();
};

#endif

