#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F4 :public Benchmarks{
	protected:
	public:
		F4();
		double compute(double* x);
		double compute(vector<double> x);
		~F4();
	};
}

