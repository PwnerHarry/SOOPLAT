#pragma once

#include "Benchmarks2010.h"

namespace CEC2010 {
	class F3 :public Benchmarks{
	protected:

	public:
		F3();
		double compute(double* x);
		double compute(vector<double> x);
		~F3();
	};
}

