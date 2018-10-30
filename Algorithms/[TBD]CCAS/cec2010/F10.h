#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F10 :public Benchmarks{
	protected:
	public:
		F10();
		double compute(double* x);
		double compute(vector<double> x);
		~F10();
	};
}

