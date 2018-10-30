#pragma once

#include "Benchmarks2010.h"

namespace CEC2010 {

	class F1 :public Benchmarks{
	public:
		F1();
		double compute(double* x);
		double compute(vector<double> x);
		~F1();
	};
}

