#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F7 :public Benchmarks{
	protected:
	public:
		F7();
		double compute(double* x);
		double compute(vector<double> x);
		~F7();
	};
}

