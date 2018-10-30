#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F9 :public Benchmarks{
	protected:
	public:
		F9();
		double compute(double* x);
		double compute(vector<double> x);
		~F9();
	};
}

