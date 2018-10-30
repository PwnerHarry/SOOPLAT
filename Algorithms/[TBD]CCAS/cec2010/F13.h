#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F13 :public Benchmarks{
	protected:
	public:
		F13();
		double compute(double* x);
		double compute(vector<double> x);
		~F13();
	};
}


